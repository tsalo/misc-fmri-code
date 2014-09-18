function meanVal = meanImageVal(valuesFile, maskFile)
% FORMAT meanVal = meanImageVal(valuesFile, maskFile)
% Gives you the mean across a ResMS file.

if exist('maskFile', 'var')
    maskedVals = spm_summarise(valuesFile, maskFile);
    meanVal = mean(maskedVals);
else
    disp('Warning. No mask given. Extracting mean from whole volume.')
    Vvals = spm_vol(valuesFile);
    imageVals = spm_read_vols(Vvals);
    imageVals(isnan(imageVals)) = 0;
    meanVal = mean(imageVals(:));
end
end
