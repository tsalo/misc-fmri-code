% Copies data in nifti image from location of one ROI to location of another, identically shaped/sized ROI.
% Take the movedRoi images that are created and use them in the second-level model.
subjects = {'subj1' 'subj2'};
subjectMasks = {'subj1_mask.nii' 'subj2_mask.nii'};
subjConFile = 'con_0001.img';
generalMask = 'dlpfc.nii';
outPath = '/whereverTheDataShouldGo/';

subjPath = '/firstHalfOfPath/';
subjSubPath = '/secondHalfOfPath/';

outPath = '/home/tsalo/';

for iSubj = 1:length(subjects)
    % Load subject's data
    vSubjCon = spm_vol([subjPath subjects{iSubj} subjSubPath subjConFile]);
    [ySubjCon, ~] = spm_read_vols(vSubjCon);
    yOut = ySubjCon;

    % Mask where you want to see it in the brain.
    vMniRoi = spm_vol(generalMask);
    [yMniRoi, ~] = spm_read_vols(vMniRoi);
    mniRoiIndex = find(yMniRoi == 1);
    [mniX, mniY, mniZ] = ind2sub(size(yMniRoi), mniRoiIndex);

    % Mask for a given subject.
    vSubjRoi = spm_vol(subjectMasks{iSubj});
    [ySubjRoi, ~] = spm_read_vols(vSubjRoi);
    subjRoiIndex = find(ySubjRoi == 1);

    % Create a barrier to make boundaries clear.
    minX = min(mniX) - 1;
    minY = min(mniY) - 1;
    minZ = min(mniZ) - 1;

    maxX = max(mniX) + 1;
    maxY = max(mniY) + 1;
    maxZ = max(mniZ) + 1;

    yOut(minX:maxX, minY:maxY, minZ:maxZ) = ones(size(yOut(minX:maxX, minY:maxY, minZ:maxZ))) .* 100;

    % Copy the values.
    yOut(mniRoiIndex) = ySubjCon(subjRoiIndex);

    % Write out new file.
    [~, fileName, fileSuffix] = fileparts(vSubjCon.fname);
    vSubjCon.fname = [outPath '/' subjects{iSubj} '_movedRoi_' fileName fileSuffix];

    vOut = spm_create_vol(vSubjCon);
    spm_write_vol(vOut, yOut);

    clear vOut vSubjCon yOut ySubjCon subjRoiIndex ySubjRoi
end
