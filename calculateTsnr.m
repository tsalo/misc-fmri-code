function outputFiles = calculateTsnr(inputFiles)
% FORMAT outputFiles = calculateTsnr(inputFiles)
% Calculate temporal signal to noise ratio for volumes in inputFiles.
% 
%
% INPUT
% inputFiles:       Functional image filenames, either
%                   1. Character array (such as from SPM.xY.P) or cell
%                      array of strings.
%                   2. SPM structure.
%
% OUTPUT
% outputFiles:      Image(s) with voxel-wise mean / SD of input images.
%
% Does not include any implicit masking.
% PW 25/03/2012
% Aesthetic changes (to match Matlab conventions) by Taylor Salo 141209.

if isstruct(inputFiles)
    if isfield(inputFiles,'Sess')
        nSessions = length(inputFiles.Sess);
        for iSession = 1:nSessions
            rowIndex = inputFiles.Sess(iSession).row;
            sessionFiles = inputFiles.xY.P(rowIndex);
            [sessionPath, ~, ~] = spm_fileparts(inputFiles(1, :));
            sessionTsnrFilename = fullfile(sessionPath, sprintf('tsnr_sess_%03d', iSession));
			if iSession == 1
				outputFiles = docalc(sessionFiles, sessionTsnrFilename);
			else
				outputFiles = char(outputFiles, docalc(sessionFiles, sessionTsnrFilename));
			end
        end
    else
        error('Input variable is a structure, but does not have field Sess, so it probably is not a valid SPM structure')
    end
elseif iscellstr(inputFiles)
    inputFiles = char(inputFiles);
    [sessionPath, ~, ~] = spm_fileparts(inputFiles(1,:));
    tsnrFilename = fullfile(sessionPath, 'tsnr.nii');
    outputFiles = docalc(inputFiles, tsnrFilename);
elseif ischar(inputFiles)
    [sessionPath, ~, ~] = spm_fileparts(inputFiles(1,:));
    tsnrFilename = fullfile(sessionPath, 'tsnr.nii');
    outputFiles = docalc(inputFiles, tsnrFilename);
else
    error('Could not figure out input variable P')
end

end

function outputFiles = docalc(inputFiles, outputFilename)
V = spm_vol(inputFiles);
Vo = V(1);
Vo.fname = outputFilename;
Vo = spm_imcalc(V, Vo, 'mean(X) ./ std(X)', {1 0 0});
outputFiles = Vo.fname;
end


