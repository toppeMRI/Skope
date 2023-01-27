% Create B0 mapping sequence for GE (for TOPPE interpeter v5)
% Modified from github/HarmonizedMRI/Calibration/b0/GE/b04ge.m

% We need the following Matlab tool(s)
addpath ~/github/toppeMRI/toppe/                     % +toppe toolbox

% GE system parameters.
% If desired (e.g., for waveform design or to reduce PNS), 
% maxSlew and maxGrad can be < system limits.
sys = toppe.systemspecs('maxSlew', 15, 'slewUnit', 'Gauss/cm/ms', ...
    'maxGrad', 5, 'gradUnit', 'Gauss/cm', ...
    'psd_rf_wait', 90, ...       % microseconds
    'psd_grd_wait', 100, ...     % microseconds
    'gradient', 'xrm', ...       % xrm: MR750 (default); hrmw: Premier; hrmb: UHP
    'timessi', 100);    % us

% Acquisition parameters
res = 0.4;  % iso voxel size (cm)
N = [60 60 60];
FOV = N*res;          % cm
flip = 5;             % degrees
deltaTE = [0 1000/440];   % (ms) acquire 2 images with TE extended by 0 and 2.3 ms

% Options
arg.entryFile = 'toppeN.entry';
arg.scanFilePath = '/usr/g/research/pulseq/cal/b0/';
arg.tbw = 8;                     % RF pulse time-bandwidth product
arg.rfDur = 2;                   % RF pulse duration (ms)
arg.ftype = 'min';               % 'min': minimum-phase SLR pulse; 'ls': linear phase
arg.slabThick = 0.8*FOV(3);      % excited slab thickness
arg.rfSpoilSeed = 117;           % RF spoiling phase increment factor (degrees)
arg.exMod         = 'tipdown.mod';
arg.readoutMod    = 'readout.mod';
arg.nCyclesSpoil = 2;   % number of cycles of phase across voxel (along x and z)
arg.fatsat       = false;         % add fat saturation pulse?
arg.fatFreqSign = -1;            % sign of fatsat pulse frequency offset
arg.autoChop = true;             % If true, only acquire data on gradient plateau. See makegre.m

nScans = numel(deltaTE);

% Since we are using the helper function 'makegre' below,
% the in-plane FOV and matrix size must be square.
if N(1) ~= N(2) | FOV(1) ~= FOV(2)
    error('In-plane FOV and matrix be square.');
end

nx = N(1);
ny = N(2);
nz = N(3);

% Write modules.txt
nModules = 2 + arg.fatsat;
fid = fopen('modules.txt', 'wt');
fprintf(fid, 'Total number of unique cores\n');
fprintf(fid, '%d\n', nModules);
fprintf(fid, 'fname  duration(us)    hasRF?  hasDAQ?\n');
if arg.fatsat
    fprintf(fid, '%s\t0\t1\t0\n', 'fatsat.mod');
end
fprintf(fid, '%s\t0\t1\t0\n', arg.exMod);
fprintf(fid, '%s\t0\t0\t1\n', arg.readoutMod);
fclose(fid);

% Write TOPPE entry file (to be placed in /usr/g/research/pulseq/ on scanner).
% This can be edited by hand as needed after copying to scanner.
toppe.writeentryfile(arg.entryFile, ...
    'filePath', arg.scanFilePath, ...
    'b1ScalingFile', arg.exMod, ...
    'readoutFile', arg.readoutMod);

% Create .mod files
if arg.fatsat
    % fat sat module (fatsat.mod)
    fatsat.flip = 90;
    fatsat.slThick = 1e5;    % dummy value (determines slice-select gradient, but we won't use it). Just needs to be large to reduce dead time before+after rf pulse
    fatsat.tbw = 2;
    fatsat.bw = 440;   % Hz
    fatsat.dur = fatsat.tbw*1e3/fatsat.bw;
    b1 = toppe.utils.rf.makeslr(fatsat.flip, fatsat.slThick, fatsat.tbw, ...
        fatsat.dur, 1e-8, sys, ...
        'ftype', 'min', 'type', 'ex', 'writeModFile', false);
    b1 = toppe.makeGElength(b1);
    toppe.writemod(sys, 'rf', b1, 'ofname', 'fatsat.mod', 'desc', 'fat sat pulse');

    fatChemShift = 3.5;  % fat/water chemical shift (ppm)
    fatFreq = arg.fatFreqSign*sys.gamma*1e4*sys.B0*fatChemShift*1e-6;  % Hz
end

% excitation module (tipdown.mod)
[ex.rf, ex.g] = toppe.utils.rf.makeslr(flip, arg.slabThick, ...
    arg.tbw, arg.rfDur, nz*arg.nCyclesSpoil, sys, ...
    'ftype', arg.ftype, ...
    'spoilDerate', 0.5, ...
    'ofname', arg.exMod);

% Data acquisition module (readout.mod).
% Here we use the helper function 'makegre' to do that, but that's not a requirement.
zres = FOV(3)/nz;
toppe.utils.makegre(FOV(1), nx, zres, sys, ... 
    'ofname', arg.readoutMod, ...
    'autoChop', arg.autoChop, ...   
    'ncycles', arg.nCyclesSpoil); 

% Write scanloop.txt
rfphs = 0;              % radians
rfSpoilSeed_cnt = 0;

toppe.write2loop('setup', sys, 'version', 4);  % initialize the file

for iz = -1:nz     % We use iz<1 for approach to steady-state
    fprintf('\b\b\b\b\b\b\b\b\b\b\b\b%d of %d', max(1,iz), nz);
    for iy = 1:ny
        for ite = 1:length(deltaTE)
            % y/z phase encode amplitudes. Turn off during approach to steady-state.
            % My convention is to start at (-kymax, -kzmax)
            a_gy = -((iy-1+0.5)-ny/2)/(ny/2) * (iz>0);  
            a_gz = -((iz-1+0.5)-nz/2)/(nz/2) * (iz>0);

            % fat saturation module
            if(arg.fatsat)
                toppe.write2loop('fatsat.mod', sys, ...
                    'RFoffset', round(fatFreq), ...   % Hz
                    'RFphase', rfphs);         % radians

                rfphs = rfphs + (arg.rfSpoilSeed/180*pi)*rfSpoilSeed_cnt ;  % radians
                rfSpoilSeed_cnt = rfSpoilSeed_cnt + 1;
            end

            % RF excitation module
            toppe.write2loop(arg.exMod, sys, ...
                'RFamplitude', 1.0, ...
                'textra', deltaTE(ite), ...
                'RFphase', rfphs);

            % data acquisition module
            toppe.write2loop(arg.readoutMod, sys, ...
                'Gamplitude', [1.0 a_gy a_gz]', ...
                'DAQphase', rfphs, ...
                'textra', max(deltaTE) - deltaTE(ite), ... % to keep TR constant
                'slice', max(iz,1), 'echo', ite, 'view', iy);

            % Update rf phase (RF spoiling)
            rfphs = rfphs + (arg.rfSpoilSeed/180*pi)*rfSpoilSeed_cnt ;  % radians
            rfSpoilSeed_cnt = rfSpoilSeed_cnt + 1;
        end
    end
end
fprintf('\n');
toppe.write2loop('finish', sys);  % finalize file

fprintf('TR = %.3f ms\n', toppe.getTRtime(1, 2, sys)*1e3);

% Create 'sequence stamp' file for TOPPE
% This file is listed in line 6 of the .entry file
toppe.preflightcheck(arg.entryFile, 'seqstamp.txt', sys);

% Write files to tar archive (for convenience).
system(sprintf('tar cf b0.tar %s seqstamp.txt scanloop.txt modules.txt *.mod', arg.entryFile));

toppe.utils.scanmsg(arg.entryFile);

% Play sequence in loop (movie) mode
%nModulesPerTR = 2;
%toppe.playseq(nModulesPerTR, sys, ...
%    'tpause', 0.05, ...
%    'nTRskip', 8);

return;


% Optional: Convert to Pulseq

% First add non-zero nChop to rf module to make Siemens happy
system('cp tipdown.mod tipdown_orig.mod');
[rf,gx,gy,gz,desc,paramsint16,paramsfloat,hdr] = toppe.readmod('tipdown.mod');
rf = [rf; zeros(48,1)];
toppe.writemod(sysGE, 'rf', rf, 'gz', gz, 'desc', desc, ...
	'nomflip', flip, ...  % needed to get right B1 scaling in .seq file
    'nChop', [48 48], 'ofname', 'tipdown.mod');

if true
	pulsegeq.ge2seq('b0.tar', sysGE, sysSiemens, 'seqFile', 'b0.seq');
else
	% test
	pulsegeq.ge2seq('b0.tar', sysGE, sysSiemens, 'seqFile', 'b0.seq', 'nt', 50);
	seq = mr.Sequence(sysSiemens);
	seq.read('b0.seq');
	b = seq.getBlock(1);
	%rep = seq.testReport;
	%fprintf([rep{:}]);
end
