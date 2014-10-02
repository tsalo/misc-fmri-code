function outFilename = erode_image(filename, THR, ERODE)
% FORMAT outFilename = erode_image(filename, THR, ERODE)
% Adapted from the Conn toolbox's mask erosion code to threshold and erode
% nifti images.
%
%
% filename: Nifti image to be binarized and eroded. String.
% THR:      Threshold to apply to file for binarization. Default is 0.5.
%           Double.
% ERODE:    Number of layers to erode. Default is 1. Double.
%
%
% Modified by Taylor Salo 140902 from code written by Alphonso Nieto-Castanon

if ~exist('THR', 'var')
    THR = 0.5;
end
if ~exist('ERODE', 'var')
    ERODE = 1;
end

V0 = spm_vol(filename);
[X0,~] = spm_read_vols(V0);

% Conn code
idx1 = find(X0(:) > THR);
if isempty(idx1)
    error('No suprathreshold voxels in ROI file.');
end
[idxx, idxy, idxz] = ind2sub(size(X0), idx1);

idxt = find((idxx > ERODE) & (idxx < (size(X0, 1) + 1 - ERODE)) & ...
            (idxy > ERODE) & (idxy < (size(X0, 2) + 1 - ERODE)) & ...
            (idxz > ERODE) & (idxz < (size(X0, 3) + 1 - ERODE)));
for n1 = 1:length(idxt)
    if (sum(sum(sum(X0(idxx(idxt(n1)) + (-ERODE:ERODE), idxy(idxt(n1)) + (-ERODE:ERODE), idxz(idxt(n1)) + (-ERODE:ERODE)) < THR, 3), 2), 1)) >= 1
        idxt(n1) = 0;
    end
end
idxt = idxt(idxt > 0);
idx1 = idx1(idxt);
X1 = zeros(size(X0));
X1(idx1) = 1;
[fpath, ffile, fext] = fileparts(V0.fname);
outFilename = fullfile(fpath, ['e' ffile fext]);
V0.fname = outFilename;
spm_write_vol(V0, X1);
end
