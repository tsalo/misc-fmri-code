function quickViewGif(inputFile, saveFile, coordinates, view)
% FORMAT quickViewGif(inputFile, saveFile, coordinates, view)
% Saves xjview screenshots cropped to fit specific views.
%
% inputFile:    Nifti file (con, spmT, etc) for which image is required.
% saveFile:     output path and filename (without extension) where image
%               will be saved.
% coordinates:  Coordinates. 1x3 double.
% view:         View. String. Options: 'sag' (sagittal), 'ax' (axial),
%               'cor' (coronal).

p = 0.001;  % P value.
k = 100;   % Minimum cluster size.

% Save image
if exist(inputFile,'file')
   fprintf(['\tFile Found: ' inputFile ': \n']);
else
   fprintf(['\tFile: ' inputFile ' DNE\n']);
   return
end

xjviewModified(inputFile, p, k, coordinates);
fg = spm_figure('FindWin','Graphics');
% Where to save screenshots
set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [.5 2.5 8.5 11.0]);
print(fg, '-noui', '-djpeg', saveFile); 
gf = imread([saveFile '.jpg']);
switch view
    case 'sag'
        gf_crop = imcrop(gf,[637, 72, 306, 259]);
        imwrite(gf_crop,[saveFile '.jpg']);
    case 'ax'
        gf_crop = imcrop(gf,[637, 341, 306, 259]);
        imwrite(gf_crop,[saveFile '.jpg']);
    case 'cor'
        gf_crop = imcrop(gf,[955, 72, 255, 259]);
        imwrite(gf_crop,[saveFile '.jpg']);
end
fprintf(['Saved ' saveFile '.jpg\n']);
close gcf
end
