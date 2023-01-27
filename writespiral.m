%function sospresto3d(sys, N, FOV, flip, nframes, varargin)
% Stack-of-spirals PRESTO fMRI sequence
% 
% WIP

% system hardware limits
sys = toppe.systemspecs('maxSlew', 13, 'slewUnit', 'Gauss/cm/ms', ...
    'maxGrad', 4, 'gradUnit', 'Gauss/cm', ...
    'myrfdel', 152, ... % psd_rf_wait
    'daqdel', 152, ...  % psd_grd_wait
    'gradient', 'hrmb', ... % xrm: MR750; hrmb: UHP; hrmw: Premier
    'B0', 3, ...        % Tesla
    'timessi', 200);    % us

res = 0.24;       % cm
N = [92 92 1];   % image matrix size
slThick = 0.5;    % cm
FOV = N*res;     % cm
FOV(3) = slThick;
fatsat = true;
fatFreqSign = -1;
nframes = 10;
nframes = 20;
flip = 15;  % degrees

toppeVersion = 5;

nleaf = 1;     % number of spirals (shots) per kz-encoding level
dsamp = 800;   % number of 4us samples in densely (fully) sampled center
Router = 3;    % undersampling factor outside densely sampled center
rewDerate = 0.8;  % derate slew rate during spiral rewinder by this factor (to control PNS)
inout = 'in';   % spiral 'in' or 'out'


% Write modules.txt
mods.ex = 'tipdown.mod';
mods.readout = 'readout.mod';
fid = fopen('modules.txt', 'wt');
fprintf(fid, 'Total number of unique cores\n');
fprintf(fid, '%d\n', length(fieldnames(mods)));
fprintf(fid, 'fname  duration(us)    hasRF?  hasDAQ?\n');
fprintf(fid, '%s\t0\t1\t0\n', mods.ex);
fprintf(fid, '%s\t0\t0\t1\n', mods.readout);
fclose(fid);


% Write entry file (to be placed in /usr/g/research/pulseq/ on scanner host)
% This can be edited by hand as needed after copying to scanner.
toppe.writeentryfile('toppeN.entry', ...
    'filePath', '/usr/g/research/pulseq/spiral/');


%% Create (balanced) stack-of-spirals readout. Isotropic resolution.
g = toppe.utils.spiral.makesosreadout(sys, N, FOV, nleaf, ...
    'dsamp', dsamp, ...
    'Router', Router, ...
    'ofname', 'tmp.mod', ...
    'inout', inout, ...
    'rewDerate', rewDerate);

system('rm tmp.mod');

% Create Router in-plane rotations, and store all rotations in the readout.mod file
% (as explicit/arbitrary waveforms).
% This makes the conversion to Pulseq easier (rotation information not retained in Pulseq file).
gc = g(:,1) + 1i*g(:,2);
nrots = Router;
n = size(g,1);
gx = zeros(n, nrots);
gy = zeros(n, nrots);
gz = zeros(n, nrots);

for irot = 1 : Router
    phi = 2*pi*(irot-1)/nrots;   % leaf rotation angle (radians)
    tmp = exp(1i*phi)*gc;
    gx(:,irot) = real(tmp);
    gy(:,irot) = imag(tmp);
    gz(:,irot) = g(:,3);
end

gx = toppe.makeGElength(gx);
gy = toppe.makeGElength(gy);
gz = toppe.makeGElength(gz);
toppe.writemod(sys, 'gx', gx, 'gy', gy, 'ofname', 'readout.mod'); % no gz blips

%readoutDur = seq.sys.raster*1e3*length(roInfo.sampWin);  % msec
%sampWin = roInfo.sampWin;


%% Create slice-selective excitation (tipdown.mod)

tbw = 6; dur = 2.5;
toppe.utils.rf.makeslr(flip, FOV(3), tbw, dur, 8, sys, ...
                            'spoilDerate', 0.5, ...
                            'type', 'st', ...
                            'ftype', 'ls', ...
                            'ofname', 'tipdown.mod');


%% Create scanloop.txt

rfphs = 0;              % radians
rfphsLast = rfphs;      % Phase of RF pulse from previous TR
daqphs = 0;
rf_spoil_seed_cnt = 0;

fprintf('Writing scanloop.txt for fMRI sequence\n');

toppe.write2loop('setup', sys, 'version', toppeVersion);

nrots = Router; 
ndisdaq = 2*nrots;

rf_spoil_seed = 117;
rf_spoil_seed_cnt = 0;

for iframe = (-ndisdaq+1):nframes
    if ~mod(iframe,10)
        fprintf([repmat('\b',1,20) sprintf('%d of %d', iframe+ndisdaq, nframes+ndisdaq)]);
    end

    nTrigPeriod = 18;  % play TTL out pulse every nTrigPeriod TRs
    if ~mod(iframe, nTrigPeriod)
        trigout = 1;
    else
        trigout = 0;
    end

    for irot = mod(iframe,nrots) + 1;  %1:nrots
        % rf excitation module (includes PRESTO gradients)
        toppe.write2loop('tipdown.mod', sys, 'RFphase', rfphs);

        % readout. Data is stored in 'slice', 'echo', and 'view' indeces.
        slice = 1;
        echo = 1;
        view = (max(iframe,1)-1)*nrots + irot;
        view = min(view, 100);
        phi = 2*pi*(irot-1)/nrots;                 % leaf rotation angle (radians)
        daqphs = rfphs;
        toppe.write2loop('readout.mod', sys, ...
            'DAQphase', daqphs, ...
            'slice', slice, 'echo', echo, 'view', view, ...
            'waveform', irot, ...
            'trigout', trigout, ...
            'Gamplitude', [1 1 0]');

        % update rf phase (RF spoiling)
        rfphsLast = rfphs;
        rfphs = rfphs + (rf_spoil_seed/180 * pi)*rf_spoil_seed_cnt ;  % radians
        rf_spoil_seed_cnt = rf_spoil_seed_cnt + 1;

    end
end
fprintf('\n');
toppe.write2loop('finish', sys);

return

%% Create 'sequence stamp' file for TOPPE.
% This file is listed in the 5th row in entryFile
% NB! The file entryFile must exist in the folder from where this script is called.
toppe.preflightcheck('toppeN.entry', 'seqstamp.txt', sys);

%% create tar file
system(sprintf("tar cf spiral2d.tar modules.txt scanloop.txt *.mod toppeN.entry seqstamp.txt %s", kspfile));

fprintf('Scan time for 3D spiral fmri sequence: %.2f min\n', toppe.getscantime(sys)/60);

%% display sequence
%toppe.playseq(2, sys, 'tpause', 1);

return

%% convert to Pulseq
if 0
% addpath ~/gitlab/toppe/pulseq/
cd tar
ge2seq('scan.tar');
seq = mr.Sequence();
seq.read('out.seq');
seq.plot_sg('TimeRange', [10 10.05]);
cd ..
end


