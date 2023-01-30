% Write 2D spiral sequence with TTL triggers,
% for testing Skope field camera using TOPPE interpreter v5.

% System hardware specifications.
% 'maxSlew' and 'maxGrad' are design parameters and can be less than
% physical limits.
sys = toppe.systemspecs('maxSlew', 13, 'slewUnit', 'Gauss/cm/ms', ...
    'maxGrad', 5, 'gradUnit', 'Gauss/cm', ...
    'psd_rf_wait', 90, ...       % microseconds
    'psd_grd_wait', 100, ...     % microseconds
    'gradient', 'xrm', ...       % xrm: MR750 (default); hrmw: Premier; hrmb: UHP
    'timessi', 100, ...
    'B0', 3);           % Tesla

% Write modules.txt
% Last column is TTL position within the module (us). 
trigpos = 1000; % us
mods.ex = 'tipdown.mod';
mods.readout = 'readout.mod';
fid = fopen('modules.txt', 'wt');
fprintf(fid, 'Total number of unique cores\n');
fprintf(fid, '%d\n', length(fieldnames(mods)));
fprintf(fid, 'fname  duration(us)  hasRF?  hasDAQ?  trigpos\n');
fprintf(fid, '%s\t0\t1\t0\t-1\n', mods.ex);
fprintf(fid, '%s\t0\t0\t1\t%d\n', mods.readout, trigpos);
fclose(fid);

% Create slice-selective excitation (tipdown.mod)
tbw = 6; dur = 2.5;
slThick = 0.5; % cm
flip = 15;  % degrees
toppe.utils.rf.makeslr(flip, slThick, tbw, dur, 8, sys, ...
                      'spoilDerate', 0.5, ...
                      'type', 'st', ...
                      'ftype', 'ls', ...
                      'ofname', mods.ex);

% Design 2D spiral (readout.mod)
res = 0.3;       % cm
N = [92 92 1];   % image matrix size
FOV = N*res;     % cm
nleaf = 1;     % number of spirals (shots) per kz-encoding level
g = toppe.utils.spiral.makesosreadout(sys, N, FOV, nleaf, ...
    'dsamp', 700, ...  % number of 4us samples in densely (fully) sampled center
    'Router', 3, ...   % undersampling factor outside densely sampled center
    'ofname', mods.readout, ...
    'inout', 'in', ...
    'rewDerate', 0.8);  % derate slew rate during spiral rewinder by this factor (to control PNS)

% Write entry file (to be placed in /usr/g/research/pulseq/ on scanner host)
% This can be edited by hand as needed after copying to scanner.
toppe.writeentryfile('loppe1.entry', ...
    'filePath', '/usr/g/research/pulseq/spiral/');

% Create scanloop.txt
rfphs = 0;              % radians
toppe.write2loop('setup', sys, 'version', 5);
nframes = 20;

for iframe = 1:nframes
    if ~mod(iframe, 10)
        fprintf([repmat('\b',1,20) sprintf('%d of %d', iframe, nframes)]);
    end

    nTrigPeriod = 15;  % play TTL out pulse every nTrigPeriod TRs
    if ~mod(iframe, nTrigPeriod)
        trigout = 1;
    else
        trigout = 0;
    end

    % rf excitation module
    toppe.write2loop(mods.ex, sys, 'RFphase', rfphs);

    % readout. Data is stored in 'slice', 'echo', and 'view' indeces.
    toppe.write2loop(mods.readout, sys, ...
        'DAQphase', rfphs, ...
        'slice', 1, 'echo', 1, 'view', iframe, ...
        'trigout', trigout, ...
        'textra', 10, ...    % add pause at end of TR (ms)
        'Gamplitude', [1 1 0]');
end

fprintf('\n');
toppe.write2loop('finish', sys);

% Create 'sequence stamp' file 
toppe.preflightcheck('toppe1.entry', 'seqstamp.txt', sys);

system(sprintf("tar cf spiral.tar modules.txt scanloop.txt *.mod toppe1.entry seqstamp.txt"));

% display sequence
toppe.plotseq(1, 2, sys);
%toppe.playseq(2, sys, 'tpause', 1);
