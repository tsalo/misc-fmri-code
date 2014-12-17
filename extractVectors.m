function extractVectors(spmFile, outDir)
% FORMAT extractVectors(spmFile, outDir)
% Creates vectors from an existing SPM.mat file.
%
%
% Created by Taylor Salo 141217

load(spmFile);

if ~exist(outDir, 'dir')
    mkdir(outDir);
end

for iSess = 1:length(SPM.Sess)
    onsets = {};
    durations = {};
    names = {};

    for jCond = 1:length(SPM.Sess(iSess).U)
        onsets = [onsets SPM.Sess(iSess).U(jCond).ons];
        durations = [durations SPM.Sess(iSess).U(jCond).dur];
        names = [names SPM.Sess(iSess).U(jCond).name{1}];
    end
    save([outDir '/vectors_' sprintf('%03d', iSess)], 'onsets', 'durations', 'names');
end

