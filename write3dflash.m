% Create B0 mapping sequence for GE (for TOPPE interpeter v5)
% Modified from github/HarmonizedMRI/Calibration/b0/GE/b04ge.m

%% PARAMETERS. EDIT THIS SECTION

addpath ~/github/toppeMRI/toppe/         % +toppe toolbox

doplot = false;   % plot sequence, and view in loop/movie mode

% GE system parameters.
% If desired (e.g., for waveform design or to reduce PNS), 
% maxSlew and maxGrad can be < system limits.
sys = toppe.systemspecs('maxSlew', 10, 'slewUnit', 'Gauss/cm/ms', ...
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
fatChemShift = 3.5;   % fat/water chemical shift (ppm)
fatFreq = -sys.gamma*1e4*sys.B0*fatChemShift*1e-6;   % Hz
deltaTE = [0 1000/abs(fatFreq)];   % (ms) acquire 2 or more images with different TE, for B0 mapping
deltaTE = [0];  % only acquire one image with 'minimum full' TE

% Options
arg.entryFile = 'toppeN.entry';
arg.scanFilePath = '/usr/g/research/pulseq/cal/b0/';
arg.tbw = 8;                     % RF pulse time-bandwidth product
arg.rfDur = 0.5;                   % RF pulse duration (ms)
arg.ftype = 'min';               % 'min': minimum-phase SLR pulse; 'ls': linear phase
arg.slabThick = 0.8*FOV(3);      % excited slab thickness
arg.rfSpoilSeed = 117;           % RF spoiling phase increment factor (degrees)
arg.exMod         = 'tipdown.mod';
arg.readoutMod    = 'readout.mod';
arg.nCyclesSpoil = 2;   % number of cycles of phase across voxel (along x and z)
arg.autoChop = true;             % If true, only acquire data on gradient plateau. See makegre.m

%% END OF PARAMETERS SECTION


% Since we are using the helper function 'makegre' below,
% the in-plane FOV and matrix size must be square.
if N(1) ~= N(2) | FOV(1) ~= FOV(2)
    error('In-plane FOV and matrix be square.');
end

nx = N(1);
ny = N(2);
nz = N(3);

% Write modules.txt
nModules = 2;
fid = fopen('modules.txt', 'wt');
fprintf(fid, 'Total number of unique cores\n');
fprintf(fid, '%d\n', nModules);
fprintf(fid, 'fname  duration(us)    hasRF?  hasDAQ?\n');
fprintf(fid, '%s\t0\t1\t0\t-1\n', arg.exMod);
fprintf(fid, '%s\t0\t0\t1\t-1\n', arg.readoutMod);
fclose(fid);

% Write TOPPE entry file (to be placed in /usr/g/research/pulseq/ on scanner).
% This can be edited by hand as needed after copying to scanner.
toppe.writeentryfile(arg.entryFile, ...
    'filePath', arg.scanFilePath, ...
    'b1ScalingFile', arg.exMod, ...
    'readoutFile', arg.readoutMod);

% Create excitation module (tipdown.mod)
[ex.rf, ex.g] = toppe.utils.rf.makeslr(flip, arg.slabThick, ...
    arg.tbw, arg.rfDur, nz*arg.nCyclesSpoil, sys, ...
    'ftype', arg.ftype, ...
    'spoilDerate', 0.8, ...
    'ofname', arg.exMod);

% Create data acquisition module (readout.mod).
% Here we use the helper function 'makegre' to do that, but that's not a requirement.
zres = FOV(3)/nz;
toppe.utils.makegre(FOV(1), nx, zres, sys, ... 
    'ofname', arg.readoutMod, ...
    'dwell', 16, ...      % ADC sample time (must be multiple of 2)
    'autoChop', arg.autoChop, ...   
    'ncycles', arg.nCyclesSpoil); 

% Write scanloop.txt
rfphs = 0;              % radians
rfSpoilSeed_cnt = 0;

toppe.write2loop('setup', sys, 'version', 5);  % initialize the file

for iz = 0:nz     % We use iz<1 for approach to steady-state
    fprintf('\b\b\b\b\b\b\b\b\b\b\b\b%d of %d', max(1,iz), nz);
    for iy = 1:ny
        for ite = 1:length(deltaTE)
            % y/z phase encode amplitudes. Turn off during approach to steady-state.
            % My convention is to start at (-kymax, -kzmax)
            a_gy = -((iy-1+0.5)-ny/2)/(ny/2) * (iz>0);  
            a_gz = -((iz-1+0.5)-nz/2)/(nz/2) * (iz>0);

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
system(sprintf('tar cf spgr.tar %s seqstamp.txt scanloop.txt modules.txt *.mod', arg.entryFile));

toppe.utils.scanmsg(arg.entryFile);

if doplot
    % Plot beginning of sequence
    toppe.plotseq(1, 2, sys);

    % Play sequence in loop (movie) mode
    nModulesPerTR = 2;
    toppe.playseq(nModulesPerTR, sys, 'tpause', 0.05, 'nTRskip', 2);
end
