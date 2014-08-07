function images = generate_spm_singletrial_newLSS(subject, spmFolder, outputdir, lump_conditions, settings)

% This script takes an existing SPM.mat file and converts it to a
% single-trial SPM.mat file, where each trial is modeled as a separate
% regressor. All other parameters (explicit mask, covariates, etc) are
% pulled in identically from the original SPM.mat file.
%
% Setting the 'estimate' flag to a non-zero value will automatically
% estimate the newly-generated SPM.mat file.
%
% A beta_info.mat file is generated, specifying condition information for
% each of the resulting betas.
%
% This script has two modes: multi-regressor and multi-model. The multi-
% regressor approach estimates a single model with all trials represented
% by individual regressors. The multi-model approach estimates a model for
% each individual trial, setting the first regressor to the individual
% trial and the second to all other trials. Beta images are then moved and
% renamed in a single betas directory (with the option to get rid of all of
% the extra files). The multi-regressor approach is similar to that
% described in Rissman et al. 2004 NI, and the multi-model approach is
% similar to the LS-S approach described in Mumford et al. 2012 NI. Note
% that the multi-model approach takes a long time to run in this
% implementation, so there's a lot of room for improvement.
%
% Requires SPM8.
%
% To-do list: contrast vectors for original conditions (for comparison
% against conventional model), option to estimate contrasts for each trial
% so that t-images can be used instead of betas, option to use
% files with a different prefix than original model (e.g., unsmoothed data)
%
% Note: The multi-regressor approach, as implemented here, has been
% debugged with my own data, but not the multi-model approach. I include
% the multi-model code here for completeness, but if you use it, debug and
% use with caution.
%
%
% Author: Maureen Ritchey, 10-2012
% Modded by Taylor Salo (140806) according to adjustments suggested by
% Jeanette Mumford. My focus was on LSS, so Rissman kind of fell by the
% wayside and I don't know if it still works. LSS now matches Turner et
% al. 2012 NI, where the design matrix basically matches the original
% design matrix (for a given block), except for the single trial being
% evaluated, which gets its own regressor. Also now you can lump multiple
% conditions. I also added some overwrite stuff so that it won't re-run
% existing subjects if you don't want it to.

%% USER INFORMATION
% Directory for the new model. If it does not exist, it will be created.
% stdir = '/whereverYouOutputTo/';

% Directory for the existing model
% modeldir = '/thatFirstPart/';
% modelSubDir = '/thatThirdPart/';
% spmFolder = [modeldir, subject, modelSubDir];

% Subjects to run
% subject = 'somebody';

% Flag for estimating model (0=no estimation, 1=estimation)
% estimate = 1;

% If desired, specify a condition name to lump together (no indiv trial
% regressors). This is useful if the original SPM.mat model included a
% 'catch-all' condition of no interest. Set to {'NONE'} if not desired.
% lump_conditions = {'BoringCondition1' 'BoringCondition2' 'BoringCondition3'};

% Type of model to run (1=multi-regressor approach, 2=multi-model approach)
% modeltype = 2;

% Option to discard extra files in multi-model approach to save space
% (1=discard, 0=leave them). Always keeps regs and covs files.
% discard_mm_files = 0;


%% MAIN CODE
modeltype = settings.modeltype;
estimate = settings.estimate;
discard_mm_files = settings.discard_mm_files;
overwrite = settings.overwrite;

% load pre-existing SPM file containing model information
fprintf('\nLoading previous model for %s:\n%s\n', subject, [spmFolder, '/SPM.mat']);
if exist([spmFolder '/SPM.mat'],'file')
    load([spmFolder '/SPM.mat'])
else
    error('Cannot find SPM.mat file.');
end

% find or create the output directory for the single-trial model
% outputdir = [stdir subject '/'];

if ~exist(outputdir, 'dir')
    fprintf('\nCreating directory:\n%s\n', outputdir);
    mkdir(outputdir)
end

% get model information from SPM file
fprintf('\nGetting model information...\n');
files = SPM.xY.P;
fprintf('Modeling %i timepoints across %i sessions.\n', size(files, 1), length(SPM.Sess));

% MULTI-MODEL APPROACH
if modeltype == 2
    % make trial directory
    betadir = [outputdir 'betas/'];
    if ~exist(betadir, 'dir')
        mkdir(betadir)
    end

    % set up beta information
    trialinfo = {'beta_number' 'session' 'condition' 'condition_rep' 'number_onsets' 'beta_name' 'trial_dir' 'beta_name'};
    counter = 1;

    % loop across sessions
    for iSess = 1:length(SPM.Sess)
        rows = SPM.Sess(iSess).row;
        sess_files = files(rows', :);
        sess_files = cellstr(sess_files);
        covariates = SPM.Sess(iSess).C.C;

        for jCond = 1:length(SPM.Sess(iSess).U)
            % check for special condition names to ignore
            if cellstrfind(SPM.Sess(iSess).U(jCond).name{1}, lump_conditions)
                if strcmp(lump_conditions{1}, 'NONE')
                else
                    fprintf('\nIgnoring conditions:\n')
                    for kCond = 1:length(lump_conditions)
                        fprintf('\t\t%s\n', lump_conditions{kCond});
                    end
                end
            % otherwise set up a model for each individual trial
            else
                for kCond = 1:length(SPM.Sess(iSess).U)
                    allOtherConds{kCond} = SPM.Sess(iSess).U(kCond).name{1};
                    allConds{kCond} = SPM.Sess(iSess).U(kCond).name{1};
                end
                allOtherConds(jCond) = [];
                otherDiffCondNames = allOtherConds;
                
                for jjCond = 1:length(SPM.Sess(iSess).U)
                    if jCond ~= jjCond
                        for jjjCond = 1:length(allOtherConds)
                            if strcmp(SPM.Sess(iSess).U(jjCond).name{1}, allOtherConds{jjjCond})
                                otherDiffCondOnsets{jjjCond} = SPM.Sess(iSess).U(jjCond).ons;
                                otherDiffCondDurations{jjjCond} = SPM.Sess(iSess).U(jjCond).dur;
                            end
                        end
                    end
                end
                
                for kTrial = 1:length(SPM.Sess(iSess).U(jCond).ons)
                    % Set onsets and durations. setdiff will reorder alphabetically/numerically,
                    % but that should not matter.
                    onsets = {};
                    durations = {};
                    names = {};
                    
                    singleName = [SPM.Sess(iSess).U(jCond).name{1} '_' sprintf('%03d', kTrial)];
                    otherSameCondName = ['OTHER_' SPM.Sess(iSess).U(jCond).name{1}];

                    singleOnset = SPM.Sess(iSess).U(jCond).ons(kTrial);
                    singleDuration = SPM.Sess(iSess).U(jCond).dur(kTrial);
                    [otherSameCondOnsets, index] = setdiff(SPM.Sess(iSess).U(jCond).ons, SPM.Sess(iSess).U(jCond).ons(kTrial));
                    otherSameCondDurations = SPM.Sess(iSess).U(jCond).dur(index);
                    
                    onsets = [onsets singleOnset otherSameCondOnsets otherDiffCondOnsets];
                    durations = [durations singleDuration otherSameCondDurations otherDiffCondDurations];
                    names = [names singleName otherSameCondName otherDiffCondNames];
                    
                    % make trial directory
                    trialdir = [outputdir 'Sess' sprintf('%03d', iSess) '/' singleName '/'];
                    if ~exist(trialdir,'dir')
                        mkdir(trialdir)
                    end

                    % add trial information
                    curinfo = {counter iSess SPM.Sess(iSess).U(jCond).name{1} kTrial...
                        length(SPM.Sess(iSess).U(jCond).ons(kTrial)) singleName trialdir ['Sess' sprintf('%03d', iSess) '_' singleName '.img']};
                    trialinfo = [trialinfo; curinfo];
                    

                    % save regressor onset files
                    regfile = [trialdir 'st_regs.mat'];
                    save(regfile, 'names', 'onsets', 'durations');

                    covfile = [trialdir 'st_covs.txt'];
                    dlmwrite(covfile, covariates, '\t');

                    % create matlabbatch for creating new SPM.mat file
                    matlabbatch = create_spm_init(trialdir, SPM);
                    matlabbatch = create_spm_sess(matlabbatch, 1, sess_files, regfile, covfile, SPM);

                    % run matlabbatch to create new SPM.mat file using SPM batch tools
                    if counter == 1
                        spm_jobman('initcfg')
                        spm('defaults', 'FMRI');
                    end
                    if overwrite || ~exist([trialdir 'beta_0001.img'], 'file')
                        fprintf('\nCreating SPM.mat file:\n%s\n\n', [trialdir 'SPM.mat']);
                        spm_jobman('serial', matlabbatch);
                        run_batches = 1;
                    else
                        fprintf('Exists: %s\n', [trialdir 'beta_0001.img']);
                        run_batches = 0;
                    end
                    counter = counter + 1;
                    
                    if estimate && run_batches % optional: estimate SPM model
                        fprintf('\nEstimating model from SPM.mat file.\n');
                        spmfile = [trialdir 'SPM.mat'];
                        matlabbatch = estimate_spm(spmfile);
                        spm_jobman('serial', matlabbatch);
                        clear matlabbatch

                        % copy first beta image to beta directory
                        copyfile([trialdir 'beta_0001.img'],[betadir 'Sess' sprintf('%03d', iSess) '_' singleName '.img']);
                        copyfile([trialdir 'beta_0001.hdr'],[betadir 'Sess' sprintf('%03d', iSess) '_' singleName '.hdr']);
                        
                        % discard extra files, if desired
                        if discard_mm_files == 1
                            prevdir = pwd;
                            cd(trialdir);
                            delete SPM*; delete *.hdr; delete *.img;
                            cd(prevdir);
                        end
                    end
                end
            end
        end
        
        % Make 4D image for each condition of interest in block.
        wantedConds = setdiff(allConds, lump_conditions);
        for jCond = 1:length(wantedConds)
            condVols = dir([betadir 'Sess' sprintf('%03d', iSess) '_' wantedConds{jCond} '*.img']);

            cellVols = struct2cell(condVols);
            cellVols = cellVols(1, :);
            for kVol = 1:length(cellVols)
                cellVols{kVol} = [betadir cellVols{kVol} ',1'];
            end
            images{jCond}{iSess} = [betadir '4D_' wantedConds{jCond} '_Sess' sprintf('%03d', iSess) '.nii'];
            matlabbatch{1}.spm.util.cat.name = [betadir '4D_' wantedConds{jCond} '_Sess' sprintf('%03d', iSess) '.nii'];
            matlabbatch{1}.spm.util.cat.vols = cellVols;
            matlabbatch{1}.spm.util.cat.dtype = 0;
            
            if overwrite || ~exist([betadir '4D_' wantedConds{jCond} '_Sess' sprintf('%03d', iSess) '.nii'], 'file')
                save([betadir '3Dto4D_jobfile.mat'], 'matlabbatch');
                spm_jobman('run', matlabbatch);
            else
                fprintf('Exists: %s\n', [betadir '4D_' wantedConds{jCond} '_Sess' sprintf('%03d', iSess) '.nii']);
            end
        end
    end

    % save beta information
    infofile = [betadir subject '_beta_info.mat'];
    save(infofile,'trialinfo');

% MULTI-REGRESSOR APPROACH
elseif modeltype == 1
    % set up beta information
    trialinfo = {'beta_number' 'session' 'condition' 'condition_rep' 'number_onsets' 'first_onset' 'beta_name'};
    counter = 1;

    % loop across sessions
    for iSess = 1:length(SPM.Sess)
        rows = SPM.Sess(iSess).row;
        sess_files = files(rows', :);
        sess_files = cellstr(sess_files);
        covariates = SPM.Sess(iSess).C.C;

        onsets = {};
        durations = {};
        names = {};

        for jCond = 1:length(SPM.Sess(iSess).U)
            % check for special condition names to lump together
            if cellstrfind(SPM.Sess(iSess).U(jCond).name{1}, lump_conditions)
                onsets = [onsets SPM.Sess(iSess).U(jCond).ons'];
                durations = [durations SPM.Sess(iSess).U(jCond).dur'];
                singleName = [SPM.Sess(iSess).U(jCond).name{1}];
                names = [names singleName];
                curinfo = {counter iSess SPM.Sess(iSess).U(jCond).name{1} [1] length(SPM.Sess(iSess).U(jCond).ons) SPM.Sess(iSess).U(jCond).ons(1) singleName};
                trialinfo = [trialinfo; curinfo];
                counter = counter + 1;
            % otherwise set up a regressor for each individual trial
            else
                for kTrial = 1:length(SPM.Sess(iSess).U(jCond).ons)
                    onsets = [onsets SPM.Sess(iSess).U(jCond).ons(kTrial)];
                    durations = [durations SPM.Sess(iSess).U(jCond).dur(kTrial)];
                    singleName = [SPM.Sess(iSess).U(jCond).name{1} '_' num2str(kTrial)];
                    names = [names singleName];
                    curinfo = {counter iSess SPM.Sess(iSess).U(jCond).name{1} kTrial length(SPM.Sess(iSess).U(jCond).ons(kTrial)) SPM.Sess(iSess).U(jCond).ons(kTrial) singleName};
                    trialinfo = [trialinfo; curinfo];
                    counter = counter + 1;
                end
            end
        end

        % save regressor onset files
        fprintf('Saving regressor onset files for Session %i: %i trials included\n', iSess, length(names));
        regfile = [outputdir 'st_regs_run' num2str(iSess) '.mat'];
        save(regfile,'names', 'onsets', 'durations');

        % save covariates (e.g., motion parameters) that were specified
        % in the original model
        covfile = [outputdir 'st_covs_run' num2str(iSess) '.txt'];
        dlmwrite(covfile, covariates, '\t');
        if ~isempty(covariates)
            for icov = 1:size(covariates, 2)
                curinfo = {counter iSess 'covariate' icov 1 0 strcat('covariate',num2str(icov))};
                trialinfo = [trialinfo; curinfo];
                counter = counter + 1;
            end
        end

        % create matlabbatch for creating new SPM.mat file
        if iSess == 1
            matlabbatch = create_spm_init(outputdir,SPM);
        end
        matlabbatch = create_spm_sess(matlabbatch,iSess,sess_files,regfile,covfile,SPM);

    end

    % save beta information
    infofile = [outputdir 'beta_info.mat'];
    save(infofile, 'trialinfo');

    % run matlabbatch to create new SPM.mat file using SPM batch tools
    fprintf('\nCreating SPM.mat file:\n%s\n', [outputdir 'SPM.mat']);
    spm_jobman('initcfg')
    spm('defaults', 'FMRI');
    spm_jobman('serial', matlabbatch);

    if estimate > 0 % optional: estimate SPM model
        fprintf('\nEstimating model from SPM.mat file.\n');
        spmfile = [outputdir 'SPM.mat'];
        matlabbatch = estimate_spm(spmfile);
        spm_jobman('serial', matlabbatch);
    end

    clear SPM matlabbatch
else
    error('Specify model type as 1 or 2');
end

clear SPM
end

%% SUBFUNCTIONS
function [matlabbatch] = create_spm_init(outputdir, SPM)
    % subfunction for initializing the matlabbatch structure to create the SPM
    matlabbatch{1}.spm.stats.fmri_spec.dir = {outputdir};
    matlabbatch{1}.spm.stats.fmri_spec.timing.units = SPM.xBF.UNITS;
    matlabbatch{1}.spm.stats.fmri_spec.timing.RT = SPM.xY.RT;
    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = SPM.xBF.T;
    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = SPM.xBF.T0;
    matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
    matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
    matlabbatch{1}.spm.stats.fmri_spec.volt = SPM.xBF.Volterra;
    matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
    if isempty(SPM.xM.VM)
        matlabbatch{1}.spm.stats.fmri_spec.mask = {''};
    else
        matlabbatch{1}.spm.stats.fmri_spec.mask = {SPM.xM.VM.fname};
    end
    matlabbatch{1}.spm.stats.fmri_spec.cvi = SPM.xVi.form;
end

function [matlabbatch] = create_spm_sess(matlabbatch, iSess, sess_files, regfile, covfile, SPM)
    % subfunction for adding sessions to the matlabbatch structure
    matlabbatch{1}.spm.stats.fmri_spec.sess(iSess).scans = sess_files; %fix this
    matlabbatch{1}.spm.stats.fmri_spec.sess(iSess).cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {});
    matlabbatch{1}.spm.stats.fmri_spec.sess(iSess).multi = {regfile};
    matlabbatch{1}.spm.stats.fmri_spec.sess(iSess).regress = struct('name', {}, 'val', {});
    matlabbatch{1}.spm.stats.fmri_spec.sess(iSess).multi_reg = {covfile};
    matlabbatch{1}.spm.stats.fmri_spec.sess(iSess).hpf = SPM.xX.K(iSess).HParam;
end

function [matlabbatch] = estimate_spm(spmfile)
    % subfunction for creating a matlabbatch structure to estimate the SPM
    matlabbatch{1}.spm.stats.fmri_est.spmmat = {spmfile};
    matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;
end

