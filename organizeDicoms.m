function organizeDicoms(dicomFolder, outFolder)
% FORMAT organizeDicoms(dicomFolder, outFolder)
% Copy and rename dicom files based on ProtocolName and SeriesNumber.
%
%
% Inputs:
% dicomFolder:  Single folder where all of the dicoms are for a single
%               appointment.
% outFolder:    Folder where dicoms will be organized (in generated
%               subfolders).

dicomFiles = dir([dicomFolder '/*.dcm']);

disp('Organizing dicoms.');
for iFile = 1:length(dicomFiles)
    dicomHeader = dicominfo([dicomFolder '/' dicomFiles(iFile).name]);
    seriesNumber = num2str(dicomHeader.SeriesNumber);
    protocolName = strrep(dicomHeader.ProtocolName, ' ', '_');
    
    outSubFolder = [seriesNumber '_' protocolName];
    
    if ~exist([outFolder '/' outSubFolder '/'], 'dir')
        outFolders{counter} = [outFolder '/' outSubFolder '/'];
        mkdir(outFolders{counter});
        counter = counter + 1;
    end
    
    system(['cp ' dicomFolder '/' dicomFiles(iFile).name ' ' outFolder '/' outSubFolder]);
end

disp('Renaming dicoms.');
for iFolder = 1:length(outFolders)
    seriesDicoms = dir([outFolders{iFolder} '/*.dcm']);
    for jFile = 1:length(seriesDicoms)
        system(['mv ' outFolders{iFolder} '/' seriesDicoms(jFile).name ' ' outFolders{iFolder} '/' sprintf('%04d', jFile)]);
    end
end
end
