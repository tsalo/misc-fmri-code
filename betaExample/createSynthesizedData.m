% createSynthesizedData
codeFile = mfilename('fullpath');
[codeDir, ~] = fileparts(codeFile);

spm('defaults', 'fMRI');

if ~exist([codeDir '/simulated_data/'], 'dir')
    mkdir([codeDir '/simulated_data/']);
end
if ~exist([codeDir '/first_level/'], 'dir')
    mkdir([codeDir '/first_level/']);
end

% Our base onsets and durations (0).
sessLength = 100; % in volumes
onsets{1} = [6 10 20 26 40 80 140 142 144 146]; % in seconds
durations{1} = 0;
names{1} = 'Impulses';
convNames = {'Regular', '3x Multiplier', '6 second Duration', '6 second Delay'};

save([codeDir '/vectors.mat'], 'names', 'onsets', 'durations')

% General settings
tr = 2;
xBF.dt = tr;
xBF.name = 'hrf';
bf = spm_get_bf(xBF);
orig_onsets = onsets{1} ./ tr;

% Create our regular convolved vector
vec = zeros(sessLength, 1);
singleCond = (int16(onsets{1}) / tr) + 1;
vec(singleCond) = 1;
for jDur = 1:int16(durations{1} / tr)
    vec(singleCond + jDur) = 1;
end

U.u = vec;
U.name = {'reg'};
convRegs{1} = spm_Volterra(U, bf.bf);

% Our 3x multiplier vector is just our regular convolved vector times 3
convRegs{2} = convRegs{1} .* 3;

% Create our 6 second duration convolved vector
durations{1} = 6;

vec = zeros(sessLength, 1);
singleCond = (int16(onsets{1}) / tr) + 1;
vec(singleCond) = 1;
for jDur = 1:int16(durations{1} / tr)
    vec(singleCond + jDur) = 1;
end

bf = spm_get_bf(xBF);

U.u = vec;
U.name = {'reg'};
convRegs{3} = spm_Volterra(U, bf.bf);

% Create our 6 second delay convolved vector
onsets{1} = onsets{1} + 6;
durations{1} = 0;

vec = zeros(sessLength, 1);
singleCond = (int16(onsets{1}) / tr) + 1;
vec(singleCond) = 1;
for jDur = 1:int16(durations{1} / tr)
    vec(singleCond + jDur) = 1;
end

xBF.dt = tr;
xBF.name = 'hrf';
bf = spm_get_bf(xBF);
U.u = vec;
U.name = {'reg'};
convRegs{4} = spm_Volterra(U, bf.bf);

% Create and save plot of our four convolved vectors
x_ax = 1:sessLength;
figure
plot(x_ax, convRegs{1}, 'b', x_ax, convRegs{2}, 'c',...
     x_ax, convRegs{3}, 'g', x_ax, convRegs{4}, 'r');
hx = graph2d.constantlineseries(orig_onsets, 'LineStyle', ':', 'Color', [1 0 1]);
changedependvar(hx, 'x');
legend(convNames{1}, convNames{2}, convNames{3}, convNames{4});
title('Convolved Vectors','FontWeight','bold');
print(gcf, '-djpeg', [codeDir '/convolvedVectors.jpeg']);
close(gcf);

% Create niftis of our simulated data, plus 100, in a block of noise
Y = zeros(40, 40, 40, length(convRegs{1}));
Y(5:35, 5:35, 5:35, :) = 100;
X = rand([40, 40, 40, length(convRegs{1})]);
mnX = mean(X(:));
X = X - mnX;
Y = Y + X;
for iConv = 1:length(convNames)
    Y(20 + iConv, 20, 20, :) = 100 + convRegs{iConv};
end

[x, y, z, nVolumes] = size(Y);

Vin = spm_vol([codeDir '/templateFile.nii']);

for iVolume = 1:nVolumes
    V(iVolume) = Vin(1);
    V(iVolume).fname = [codeDir '/simulated_data/volume_' sprintf('%03d', iVolume) '.nii'];
    V(iVolume).dim = [x y z];
    V(iVolume).private.dat.dim = [x y z nVolumes];
    V(iVolume).private.dat.fname = V(iVolume).fname;
    
    V2 = spm_create_vol(V(iVolume));
    spm_write_vol(V2, Y(:, :, :, iVolume));
    
    scans{iVolume} = [V(iVolume).fname ',1'];
end
clear V Y

% Run very basic first level GLM with regular vectors
load([codeDir '/templateJobfile.mat']);
matlabbatch{1}.spm.stats.fmri_spec.sess.scans = scans;
matlabbatch{1}.spm.stats.fmri_spec.sess.multi{1} = [codeDir '/vectors.mat'];
matlabbatch{1}.spm.stats.fmri_spec.dir{1} = [codeDir '/first_level/'];
matlabbatch{1}.spm.stats.fmri_spec.timing.RT = tr;
matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = names{1};
save([codeDir '/jobfile.mat'], 'matlabbatch');

spm_jobman('initcfg');
spm_jobman('run',matlabbatch);
close(gcf);

V = spm_vol([codeDir '/first_level/beta_0001.img']);
[Y, ~] = spm_read_vols(V);

outStruct{1}.header{1} = 'Convolution Version';
outStruct{2}.header{1} = 'Beta Value';
outStruct{3}.header{1} = 'Coordinates';
for iConv = 1:length(convNames)
    coords{iConv} = [20 + iConv, 20, 20] * V.mat(1:3, 1:3) + V.mat(1:3, 4)';
    outStruct{1}.col{iConv} = convNames{iConv};
    outStruct{2}.col{iConv} = Y(20 + iConv, 20, 20);
    outStruct{3}.col{iConv} = num2str(coords{iConv});
    fprintf([convNames{iConv} ':\t\t' num2str(coords{iConv}) '\n']);
end

writeCsv(outStruct, [codeDir '/betas.csv']);

% Now, just use SPM to plot those four coordinates to get an idea of how
% well the GLM fit the convolved vectors:
% Load the contrast in SPM Results (with p <= 1, to be safe).
% Navigate to one of the four coordinates.
% Hit "Plot" --> "Fitted Responses" --> "Adjusted" --> "scan or time" --> SAVE. THAT. FIGURE.
