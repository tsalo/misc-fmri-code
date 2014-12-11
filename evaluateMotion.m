function [outMotionData, stimulusCorrelatedMotion] = evaluateMotion(motionFile, vectorFile, tr)
% FORMAT [outMotionData, stimulusCorrelatedMotion] = evaluateMotion(motionFile, vectorFile, tr)
% Calculates framewise displacement and stimulus-correlated motion based on
% condition onsets/durations, motion information provided by SPM8's
% realignment tool, and TR.
% Author: Taylor Salo 141002-141209
%
%
% Inputs:
% motionFile:               String pointing to rp*.txt file outputted by
%                           SPM's realignment tool. Contains nVolumes X 6
%                           array.
% vectorFile:               String pointing to mat file where single
%                           block's vector information is stored.
% tr:                       Time to repetition, in seconds. Double.
%
% Outputs:
% outMotionData:            Motion parameter array from motionFile with
%                           framewise displacement appended as 7th column.
% stimulusCorrelatedMotion: Cell array of structures. One cell for each
%                           condition. Structure contains name field with
%                           condition name and corr field with correlation
%                           between expected BOLD response and motion.

motionData = importdata(motionFile);

% In Power et al. 2012, "rotational displacements were converted from
% degrees to millimeters by calculating displacement on the surface of a
% sphere of radius 50 mm, which is approximately the mean distance from the
% cerebral cortex to the center of the head." Rotation is given in radians
% by SPM8, and the formula for length of an arc from radians and radius is
% length = radians * radius.
motionDataDisp = motionData;
motionDataDisp(:, 4:6) = motionData(:, 4:6) .* 50;
motionDataDeriv = [0 0 0 0 0 0; diff(motionDataDisp)];
framewiseDisplacement = sum(abs(motionDataDeriv), 2);
motionDataDerivAndFd = [motionDataDeriv framewiseDisplacement];
[nScans, ~] = size(framewiseDisplacement);
outMotionData = [motionData framewiseDisplacement];
vectors = load(vectorFile);

% Quantify Stimulus-Correlated Motion.
% Convolution method from spm.martinpyka.de/?p=41
stimulusCorrelatedMotion = cell(size(vectors.onsets));
for iCond = 1:length(vectors.onsets)
    vec = zeros(nScans, 1);
    singleCond = (int16(vectors.onsets{iCond}) / tr) + 1;
    vec(singleCond) = 1;
    for jDur = 1:int16(vectors.durations{iCond} / tr)
        vec(singleCond + jDur) = 1;
    end
    
    xBF.dt = tr;
    xBF.name = 'hrf';
    bf = spm_get_bf(xBF);
    
    U.u = vec;
    U.name = {'reg'};
    convreg = spm_Volterra(U, bf.bf);
    
    stimulusCorrelatedMotion{iCond}.name = vectors.names{iCond};
    for jMot = 1:7
        stimulusCorrelatedMotion{iCond}.corr(jMot) = abs(corr(convreg, motionDataDerivAndFd(:, jMot)));
    end
end
end
