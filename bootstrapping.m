%% Bootstrapping

% paths
addpath(genpath('/Users/hlutwak/Documents/MATLAB/'))
addpath(genpath('/Applications/Psychtoolbox'))
addpath('/Users/hlutwak/Documents/GitHub/VRopticflow/Analysis')
dataFolder = '/Users/hlutwak/Documents/GitHub/VRopticflow/Data';
figFolder = '/Users/hlutwak/Documents/GitHub/VRopticflow/Figures';
analysisFolder = '/Users/hlutwak/Documents/GitHub/VRopticflow/Analysis';

% names of files
S = dir(fullfile(dataFolder,'*.mat'));

% which subjects data to analyze
subjects = ["MG"]; %,"DL", "PL","MG", "SM", "IK", "JO", "KZ"]; %"MP","DL", "PL","MG", "SM", "IK", "JO", "KZ","IG"
s = 1;
data_set = 1; %1.5 deg threshold

stims =  ["full-1", "full-2"]; %, "monocular-2"]; %add "copy" to have pa.good_trials, and/or dconst and dsurround based on vertical eye movements
ideal_eye = 1; % use measurements of data_const and data_surr based on ideal eye movements, otherwise use eyetracking vertical movements
depth_est = 0;
depth_range = 1.05; % multiplicative

nboot = 10;
boot_dev_const = zeros(nboot, 1);
boot_dev_surr = zeros(nboot, 1);


%load appropriate files
%% bootsrap over condition
tic
for b = 1:nboot
data = [];
data_const = [];
data_surr = [];
count = 0;

for f = 1:length(S)
   subj = contains(S(f).name,subjects(s));
   sti = contains(S(f).name,stims) && contains(S(f).name, 'eyetracking');

        if subj && sti
            load(fullfile(dataFolder,S(f).name));
            display(fullfile(dataFolder,S(f).name));

            if ideal_eye
                % get all trials of this condition, resample and
                % recalculate number of correct/incorrect
                n_conditions = length(pa.speed)*length(pa.direction);
                conditions = fullfact([numel(pa.speed), numel(pa.direction)]);
                data_session = [];
                for cond = 1:n_conditions
                    idx_speed = find(pa.fullFactorial(:,3) == pa.speed(conditions(cond,1)));
                    idx_direction = find(pa.fullFactorial(:,4) == pa.direction(conditions(cond,2)));
                    idx = intersect(idx_speed, idx_direction);
                    if isfield(pa, 'goodTrials')
                        idx = intersect(idx, pa.goodTrials{data_set});
                    end

                    % resample 
                    outcome = eq(pa.LR(idx), pa.LRresponse(idx));
                    resampled_ncorrect = sum(datasample(outcome, length(outcome)));

                    data_session(cond,:) = [nan(1) nan(1) pa.speed(conditions(cond,1)), rad2deg(pa.direction(conditions(cond,2))), resampled_ncorrect, length(idx)];

                end
%                 pa.speed = [0.3];
%                 pa.direction = deg2rad(90);
                [dconst, dsurr] = DistanceToConstraint(ds, pa, depth_range, depth_est);
                data_session(:,1) = dconst(:);
                data_session(:,2) = dsurr(:);
                data = [data; data_session];
                data_const = [data(:,1) data(:,end-1:end)]; % to surround data_const = [data(:,1) data(:,end-1:end)];
                data_surr = [data(:,2) data(:,end-1:end)];
                count = count+1;


            else
                
                pcorrect = eq(pa.LR(pa.goodTrials{data_set}), pa.LRresponse(pa.goodTrials{data_set}));
                pcorrect = +pcorrect;
                nTrials = ones(size(pcorrect));
                
                data_const = [data_const; pa.dconst{data_set}; pcorrect; nTrials]';
                data_surr = [data_surr; pa.dsurr{data_set}; pcorrect; nTrials]';
            end
        end

end



    % disp([pa.subjectName, ' goodtrials = ', num2str(sum(data_const(:,end))/(length(data_const)*20)*100), '%']) %pa.nRepeats=20 except MP in one environment


    options             = struct;   % initialize as an empty struct
    options.sigmoidName = 'weibull'; %'weibull';
    options.expType     = '2AFC';   % choose 2-AFC as the paradigm of the experiment
    % this sets the guessing rate to .5 and
    % fits the rest of the parameters
    options.fixedPars = NaN(5,1);
    
    if ideal_eye
%         4 speeds, 6 directions
        options.dataColor = [255,153,255; 255,102,255; 255,51,255; 204,0,204;
            255,153,153; 255,102,102; 255,51,51; 204,0,0;
            255,204,153; 255,178,102; 255, 153, 51; 204,102,0;
            204,255,153; 178,255,102; 153,255,51; 102,204,0;
            153,255,255; 102,255,255; 51,255,255; 0,204,204;
            153,153,255; 102,102,255; 51,51,255; 0,0,204]/255;
    % options.poolxTol = 0.0025;
    else
        options.poolxTol = 0.005;

    end

    result_const = psignifit(data_const,options);
    result_surr = psignifit(data_surr,options);

    boot_dev_const(b) = result_const.deviance;
    boot_dev_surr(b) = result_surr.deviance;


end
toc


%% resample over whole dataset

data = [];
count = 0;

% get data
for f = 1:length(S)
   subj = contains(S(f).name,subjects(s));
   sti = contains(S(f).name,stims) && contains(S(f).name, 'eyetracking');

        if subj && sti
            load(fullfile(dataFolder,S(f).name));
            display(fullfile(dataFolder,S(f).name));

               goodTrials = pa.goodTrials{data_set};
               % want data with speed, direction, correct/incorrect
               outcome = eq(pa.LR(goodTrials), pa.LRresponse(goodTrials));
               data_session = [pa.fullFactorial(goodTrials,3), ...
                                pa.fullFactorial(goodTrials,4), outcome'];
               data = [data; data_session];

        end

end
data_const = [];
data_surr = [];
[dconst, dsurr] = DistanceToConstraint(ds, pa, depth_range, depth_est);

% for n boot
for b = 1:nboot

% resample data
resample_idx = datasample(1:length(data), length(data));

n_conditions = length(pa.speed)*length(pa.direction);
conditions = fullfact([numel(pa.speed), numel(pa.direction)]);

for cond = 1:n_conditions
    idx_speed = find(data(:,1) == pa.speed(conditions(cond,1)));
    idx_direction = find(data(:,2) == pa.direction(conditions(cond,2)));
    idx = intersect(idx_speed, idx_direction);

    data_conditions(cond,:) = [nan(1) nan(1) pa.speed(conditions(cond,1)), rad2deg(pa.direction(conditions(cond,2))), sum(data(idx,end)), length(idx)];
end

data_conditions(:,1) = dconst(:);
data_conditions(:,2) = dsurr(:);

data_const = [data_conditions(:,1) data_conditions(:,end-1:end)]; % to surround data_const = [data(:,1) data(:,end-1:end)];
data_surr = [data_conditions(:,2) data_conditions(:,end-1:end)];

    options             = struct;   % initialize as an empty struct
    options.sigmoidName = 'weibull'; %'weibull';
    options.expType     = '2AFC';   % choose 2-AFC as the paradigm of the experiment
    % this sets the guessing rate to .5 and
    % fits the rest of the parameters
    options.fixedPars = NaN(5,1);
    
    if ideal_eye
%         4 speeds, 6 directions
        options.dataColor = [255,153,255; 255,102,255; 255,51,255; 204,0,204;
            255,153,153; 255,102,102; 255,51,51; 204,0,0;
            255,204,153; 255,178,102; 255, 153, 51; 204,102,0;
            204,255,153; 178,255,102; 153,255,51; 102,204,0;
            153,255,255; 102,255,255; 51,255,255; 0,204,204;
            153,153,255; 102,102,255; 51,51,255; 0,0,204]/255;
    % options.poolxTol = 0.0025;
    else
        options.poolxTol = 0.005;

    end

    result_const = psignifit(data_const,options);
    result_surr = psignifit(data_surr,options);

boot_dev_const(b) = result_const.deviance;
boot_dev_surr(b) = result_surr.deviance;

end

