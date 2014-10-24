function writeCsv(data, loc, varargin)
% FORMAT writeCsv(data, loc, varargin)
% Description: This function accepts a data structure and and save location
% and file name e.g. /home/name/output.csv. The data structure is a sort of
% pseudo-standard. The structure should have an area for a .header and .col
% field. This is requisite. Essentially the format should mirror that of an
% excel file. I think this should work for any data set in the given
% format.
%
% TLDR
% This function saves structures as .csv files.
% data   => structure with .header and .col fields
% loc    => save location and filename
% B.R. Geib Winter 2012

default='w+'; nohead=0;
if ~isempty(varargin),
	default = varargin{1}; 
	if exist(loc, 'file') == 0 
		nohead = 0;
		default = 'w+';
	else
		nohead = 1;
	end
end
fid = fopen(loc, default);

% Print out the headers first (if it's a new file)
if nohead == 0
	for i = 1:length(data)
		current_class = class(data{i}.header);
		switch current_class
			case 'cell'
				fprintf(fid, [char(data{i}.header) ',']);
			case 'char'
				fprintf(fid, [data{i}.header ',']);
			otherwise
				fprintf(fid, [num2str(data{i}.header) ',']);
		end
	end
	fprintf(fid, '\n');
end

% Check the field col
if ~isfield(data{1}, 'col')
	fprintf(['Error writing: ' loc '\n']);
	fprintf('\tData does not appear to exist with the variable "data"\n');
	return;
end

% Examine the data structure
for j = 1:length(data{1}.col)
    for i = 1:length(data)
        if ~isfield(data{i}, 'col')
            fprintf(['Error writing: ' loc '\n']);
            fprintf('\tData does not appear to exist with the variable "data"\n');
			return;
        end
        % Not all columns are demanded to be the same length. If a column
        % is empty at the end we want it to just print as empty.
        L = length(data{i}.col);
        if j <= L
            % Determine the class of the {} or (), we don't know which it
            % should be
            current_class1 = class(data{i}.col(j));
            % Here we try out both classes, an error in one causes a
            % default occurance of the second.
            try
                % If the object is empty, print as such.
                if ~isempty(data{i}.col(j))
                    switch current_class1
                        case 'cell'
                            fprintf(fid, char(data{i}.col(j)));
                        case 'char'
                            fprintf(fid, data{i}.col(j));
                        otherwise
                            fprintf(fid, num2str(data{i}.col(j)));
                    end
                end
            catch err
                current_class2 = class(data{i}.col{j});
                if ~isempty(data{i}.col{j})
                    switch current_class2
                        case 'cell'
                            fprintf(fid, char(data{i}.col{j}));
                        case 'char'
                            fprintf(fid, data{i}.col{j});
                        otherwise
                            fprintf(fid, num2str(data{i}.col{j}));
                    end
                end
            end
        end
        fprintf(fid, ',');
    end
    fprintf(fid, '\n');
end

fclose(fid);
fprintf(['Saved: ' loc '\n']);
end
