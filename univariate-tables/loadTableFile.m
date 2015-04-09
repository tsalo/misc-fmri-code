function data3 = loadTableFile(textFile)
% FORMAT data3 = loadTableFile(textFile)
% Loads the Pickatlas text files.
% Dependencies: readCsv.m, removeEmptyCells.m
TextStruct = readCsv(textFile);
data = TextStruct{1}.col;
data2 = cellfun(@strsplit, data, repmat({'\t'}, size(data)), 'UniformOutput', 0);
data3 = cell(length(data2), length(data2{1}));

for iRow = 1:length(data2)
    data2{iRow} = removeEmptyCells(data2{iRow});
    data3(iRow, :) = data2{iRow};
end
end
