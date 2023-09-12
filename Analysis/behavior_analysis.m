%% Analysis code
% get behavior files, conccatenate same trials, run psignifit

% add psignifit toolbox
addpath('/Users/hopelutwak/Documents/MATLAB/psignifit')

% assign data folder
dataFolder = '/Users/hopelutwak/Documents/GitHub/VRopticflow/Data';

% names of files
S = dir(fullfile(dataFolder,'*.mat'));

% which subjects data to analyze
subjects = ["IK"]; %"ABC", "HL","MR", "KB", "KZ"
stims = "pilot";

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
            for cond = 1:n_conditions
                idx_speed = find(pa.fullFactorial(:,3) == pa.speed(conditions(cond,1)));
                idx_direction = find(pa.fullFactorial(:,4) == pa.direction(conditions(cond,2)));
                idx = intersect(idx_speed, idx_direction);
                data_session(cond,:) = [pa.speed(conditions(cond,1)), rad2deg(pa.direction(conditions(cond,2))), sum(eq(pa.LR(idx), pa.LRresponse(idx))), pa.nRepeats];
            end
            data = [data; data_session];
            count = count+1;
        end
    end

    
end