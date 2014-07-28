function clustsim(path, pthr, athr, extras, open_file)
% FORMAT clustsim(path, pthr, athr, extras, open_file)
% Runs AFNI's 3dClustSim (Monte Carlo simulations) to determine optimal 
% cluster size threshold for second-level whole brain results.
% Uses 3dcalc to create sqrt_ResMS image, then uses 3dFWHMx to obtain
% smoothness of the noise in x, y, and z dimensions. Finally, it runs
% 3dClustSim (10000 iterations) to get a table of minimum cluster sizes for
% different p and alpha thresholds. 
%
%
% path:             Path to folder containing second level SPM.mat.
%                   String.
% pthr:             Optional. P threshold(s) at which clustsim will be run. 
%                   Double vector.
% athr:             Optional. Alpha value(s) for which clustsim will 
%                   present results. Double vector.
% extras:           Optional. Additional inputs (e.g. -nodec -2sided
%                   -quiet). String.
% open_file:        Optional. Chooses whether or not to open the
%                   results file at the end. 1 to open file, 0 to not.
    
if ~exist(path,'dir')
    fprintf('Directory does not exist: %s\n',path);
    return
end

if ~exist('pthr','var')
    fprintf('pthr is empty, setting to default\n');
    pthr = [0.05 0.01 0.005 0.001];
    %pthr = [0.001];
elseif ischar(pthr)
    fprintf('pthr is string, setting to default\n');
    pthr = [0.05 0.01 0.005 0.001];
end

if ~exist('athr','var')
    fprintf('athr is empty, setting to default\n');
    athr = 0.05;
elseif ischar(athr)
    fprintf('athr is string, setting to default\n');
    athr = 0.05;
end

if ~exist('extras','var')
    extras = '';
end

if ~exist(open_file, 'var')
    open_file = 0;
end

system(['/nfs/pkg64/afni/3dcalc -a ' path '/ResMS.hdr -expr "sqrt(a)" -prefix ' path '/sqrt_ResMS.nii -overwrite']);
system(['/nfs/pkg64/afni/3dFWHMx -mask ' path '/mask.hdr -input ' path '/sqrt_ResMS.nii -output ' path '/fwhm.txt -overwrite']);
fid = fopen([path '/fwhm.txt']);
fwhm = textscan(fid,'%f');
fwhmx = fwhm{1}(1);
fwhmy = fwhm{1}(2);
fwhmz = fwhm{1}(3);

delete([path '/fwhm.txt']);

system(['/nfs/pkg64/afni/3dClustSim -mask ' path '/mask.hdr -iter 10000 -pthr ' num2str(pthr) ' -athr ' num2str(athr) ' -fwhmxyz ' ... 
    num2str(fwhmx) ' ' num2str(fwhmy) ' ' num2str(fwhmz) ' ' extras ' -prefix ' path '/montecarlo']);

if open_file == 1
    system(['emacs22-gtk ' path '/montecarlo.NN1.1D']);
else
    fprintf('open_file is not 1, not opening output file.\n');
end
end

