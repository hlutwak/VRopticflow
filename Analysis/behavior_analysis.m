%% Analysis code
% get behavior files, conccatenate same trials, run psignifit

% add psignifit toolbox
% addpath('/Users/hopelutwak/Documents/MATLAB/psignifit')
addpath(genpath('C:\Users\hlutw\OneDrive\Documents\MATLAB\psignifit-master'))
% assign data folder
% dataFolder = '/Users/hopelutwak/Documents/GitHub/VRopticflow/Data';
dataFolder = 'C:\Users\hlutw\OneDrive\Documents\GitHub\VRopticflow\Data';

% names of files
S = dir(fullfile(dataFolder,'*.mat'));

% which subjects data to analyze
subjects = ["HL"]; %"HL" "IK"
stims = "eyesim"; %"pilot"
depth_range = .1;

% loop over all subjects


for s  = 1:length(subjects)

    data = [];
    %load appropriate files
    count = 0;
    for f = 1:length(S)
        subj = contains(S(f).name,subjects(s));
        sti = contains(S(f).name,stims);
        
        if subj && sti
            load(fullfile(dataFolder,S(f).name));
            n_conditions = length(pa.speed)*length(pa.direction);
            conditions = fullfact([numel(pa.speed), numel(pa.direction)]);
            data_session = [];
            for cond = 1:n_conditions
                idx_speed = find(pa.fullFactorial(:,3) == pa.speed(conditions(cond,1)));
                idx_direction = find(pa.fullFactorial(:,4) == pa.direction(conditions(cond,2)));
                idx = intersect(idx_speed, idx_direction);
                data_session(cond,:) = [nan(1) nan(1) pa.speed(conditions(cond,1)), rad2deg(pa.direction(conditions(cond,2))), sum(eq(pa.LR(idx), pa.LRresponse(idx))), pa.nRepeats];
            end
            [dconst, dsurr] = DistanceToConstraint(ds, pa, depth_range);
            data_session(:,1) = dconst(:);
            data_session(:,2) = dsurr(:);
            data = [data; data_session];
            data_const = [data(:,1) data(:,end-1:end)];
            data_surr = [data(:,2) data(:,end-1:end)];
            count = count+1;
            
        end

    end

    
    options             = struct;   % initialize as an empty struct
    options.sigmoidName = 'weibull';   
    options.expType     = '2AFC';   % choose 2-AFC as the paradigm of the experiment
                                    % this sets the guessing rate to .5 and
                                    % fits the rest of the parameters
    options.fixedPars = NaN(5,1);                                
    options.fixedPars(5) = 0;       % fix eta (dispersion) to zero

    result_const = psignifit(data_const,options);
    result_surr = psignifit(data_surr, options);
    figure, plotPsych(result_const, options);
    title(['distance to constraint, depth range = ', num2str(depth_range)])
    figure, plotPsych(result_surr, options);
    title('distance to surround')
    
    

    
end