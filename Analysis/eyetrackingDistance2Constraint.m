%% incorporating eye movement from eye tracking into distance to constraint calculation

%% get eye tracking interval
clear
% close all

% addpath('/Users/hopelutwak/Documents/MATLAB/psignifit')
% addpath(genpath('/Users/hopelutwak/Documents/GitHub/VRopticflow/Analysis'))
% dataFolder = '/Users/hopelutwak/Documents/GitHub/VRopticflow/Data';

addpath('/Users/hlutwak/Documents/MATLAB/psignifit')
addpath(genpath('/Users/hlutwak/Documents/GitHub/VRopticflow/Analysis'))
dataFolder = '/Users/hlutwak/Documents/GitHub/VRopticflow/Data';

% 
% % which subjects data to analyze
subjects = "KZ"; %"HL" "IK"
stims = ["monocular-2"]; %["full-1", "full-2"]; %"pilot"

D=dir('Data/');

for d = 1:length(D)
    subj = contains(D(d).name, subjects);
    stim = contains(D(d).name, stims);
    eyeExist = D(d).isdir;
    eyetracking = contains(D(d).name, "eyetracking");
    if subj && stim && eyeExist
        filename = D(d).name;
        disp(filename)
    elseif subj && stim && eyetracking
        load(fullfile(dataFolder,D(d).name));
        disp("loading: " + D(d).name);
    end
end

% filename = '2024-02-26_15-15-38_MP-full-1';

gaze = readtable(['Data/', filename, '/gaze.csv']);
blinks = readtable(['Data/', filename, '/blinks.csv']);
evts = readtable(['Data/', filename, '/events.csv']);
% t= table2array(worldtime(:,end));
timestamps = table2array(gaze(:,3));
blink_start = table2array(blinks(:,4));
blink_end = table2array(blinks(:,5));
evt = table2array(evts(:,2));

x = table2array(gaze(:,9));
y = table2array(gaze(:,10));

if sum(y) == 0

    x = table2array(gaze(:,4));
    y = table2array(gaze(:,5));
end

figure, scatter(x,y), axis equal

% in datetime
date_time = datetime(timestamps/1e9, 'ConvertFrom', 'posixtime', 'TimeZone','local', 'Format', 'd-MMM-y HH:mm:ss:ms');
blink_start_time = datetime(blink_start/1e9, 'ConvertFrom', 'posixtime', 'TimeZone','local', 'Format', 'd-MMM-y HH:mm:ss:ms');
blink_end_time= datetime(blink_end/1e9, 'ConvertFrom', 'posixtime', 'TimeZone','local', 'Format', 'd-MMM-y HH:mm:ss:ms');
evts_time = datetime(evt/1e9, 'ConvertFrom', 'posixtime', 'TimeZone','local', 'Format', 'd-MMM-y HH:mm:ss:ms');


figure, plot(date_time,x, 'linewidth', 2), hold on, plot(date_time,y,'linewidth', 2)
yl = ylim;
hold on, line([oc.UTCtrialEnd; oc.UTCtrialEnd],[yl(1); yl(2)].*ones(size(oc.UTCtrialEnd)), 'color','r') 
hold on, line([oc.UTCtrialStart; oc.UTCtrialStart],[yl(1); yl(2)].*ones(size(oc.UTCtrialStart)), 'color','g') 
hold on, line([evts_time'; evts_time'], [yl(1); yl(2)].*ones(size(evts_time')), 'color','k') 


% offset
idx = find(strcmpi(evts.name, 'trialstart')==1); %find trial start event

offset = milliseconds(evts_time(idx) - oc.UTCtrialStart(2));

% subtract offset from eyetracking data
synced = date_time-milliseconds(offset);
synced_evts = evts_time -milliseconds(offset);

figure, plot(synced,x, 'linewidth', 2), hold on, plot(synced,y,'linewidth', 2)
yl = ylim;
hold on, line([oc.UTCtrialEnd; oc.UTCtrialEnd],[yl(1); yl(2)].*ones(size(oc.UTCtrialEnd)), 'color','r') 

hold on, line([oc.UTCtrialStart; oc.UTCtrialStart],[yl(1); yl(2)].*ones(size(oc.UTCtrialStart)), 'color','g') 
hold on, line([synced_evts'; synced_evts'], [yl(1); yl(2)].*ones(size(evts_time')), 'color','k') 

title('offset')

% detrend_y = detrend(y);
% detrend_x = detrend(x);
% figure, plot(synced,detrend_x, 'linewidth', 2), hold on, plot(synced,detrend_y,'linewidth', 2)
% yl = ylim;
% hold on, line([oc.UTCtrialEnd; oc.UTCtrialEnd],[yl(1); yl(2)].*ones(size(oc.UTCtrialEnd)), 'color','r') 
% 
% hold on, line([oc.UTCtrialStart; oc.UTCtrialStart],[yl(1); yl(2)].*ones(size(oc.UTCtrialStart)), 'color','g') 
% hold on, line([synced_evts'; synced_evts'], [yl(1); yl(2)].*ones(size(evts_time')), 'color','k') 
% title('detrend')

%
% butterworth filter
fs = 200; %round(1/seconds(mean(diff(synced))))
flow = .1;
fhigh = 60;
% [b,a]=butter(4,1/60);
[b,a] = butter(4,[flow, fhigh]/(fs/2));
% output_datax= filter(b,a,detrend_x);
% output_datay = filter(b,a,detrend_y);
output_datax= filtfilt(b,a,x);
output_datay = filtfilt(b,a,y); %zero phase shift

figure, set(gcf,'renderer','Painters'), plot(synced,output_datax, 'linewidth', 2), hold on, plot(synced,output_datay,'linewidth', 2)
yl = ylim;
hold on, line([oc.UTCtrialEnd; oc.UTCtrialEnd],[yl(1); yl(2)].*ones(size(oc.UTCtrialEnd)), 'color','r') 

hold on, line([oc.UTCtrialStart; oc.UTCtrialStart],[yl(1); yl(2)].*ones(size(oc.UTCtrialStart)), 'color','g') 
hold on, line([synced_evts'; synced_evts'], [yl(1); yl(2)].*ones(size(evts_time')), 'color','k') 

title('offset and filtered, filtfilt')

% figure
% freqz(b,a,[],fs)

% 
x = output_datax;
y= output_datay;


% eye data within intervals
trial_times = [];
eyetracking = [];
gaze_speed = NaN(1, pa.trialNumber-1);
start_position = NaN(2, pa.trialNumber-1);
end_position = NaN(2, pa.trialNumber-1);


% get eyetracking for calibration time
calib = find(synced < oc.UTCCalibrationEnd);
figure, hold on,scatter(x(calib), y(calib));
axis equal

% get calibration points
idx = find(strcmpi(evts.name, 'calibrationpoint')==1);
cali_times = synced_evts(idx);

for t = 2:pa.nTrials %pa.trialNumber %full set, change to pa.nTrials
    tf = isbetween(synced, oc.UTCtrialStart(t), oc.UTCtrialEnd(t));
    trial_times = [trial_times; synced(tf)];
    eyetracking = [eyetracking; x(tf), y(tf)];
    idx = find(tf);
    if ~isempty(idx)
    %     gaze_speed(t) = vecnorm([x(idx(1)) y(idx(1))]-[x(idx(end)) y(idx(end))])/.5;
        gaze_speed_thistrial = vecnorm([diff(x(tf)) diff(y(tf))]')./seconds(diff(synced(tf)))';
%         gaze_speed(t) = max(gaze_speed_thistrial);
        gaze_speed(t) = median(gaze_speed_thistrial);
        
        [peak, peak_idx] = max(y(idx));
        [trough, trough_idx] = min(y(idx));
        start_position(:,t) = [x(idx(peak_idx)); y(idx(peak_idx))];
        end_position (:,t) = [x(idx(trough_idx)); y(idx(trough_idx))];
        hold on, scatter(x(tf), y(tf))
    else
        warning(['recording ended at ' num2str(t)])
        break
    end
end

% remove trials with high gaze speed
reasonable_speed = find(gaze_speed<100);


% % find startpos by getting half second before first trial starts
% t = 1;
% [val, idx] = min(abs(synced - oc.UTCtrialStart(t)));
% interval = .5;
% [val, interval_start] = min(abs(synced(idx) - synced-seconds(.5)));
% startpos = [mean(x(interval_start:idx)), mean(y(interval_start:idx))];
% hold on, scatter(startpos(1), startpos(2), 100, 'filled', 'b')

% average start and end position
startpos = nanmean(start_position(:,reasonable_speed)');
endpos = nanmean(end_position(:,reasonable_speed)');
hold on, scatter(startpos(1), startpos(2), 'filled', 'g')
hold on, scatter(endpos(1), endpos(2), 'filled', 'r')


% find good trials, no large horizontal variations

good_trials = [];
bad_trials = [];
target_comparison = [];
response_comparison = [];
% startpos = [5.4, 15.4]; %guess MP [5,-1.45]; DL [6,14] [8,12], MG [5.4, 15.4]
% % endpos = [5,11]; % DL [5,11], MG [5.2, 11.2]
hdist_tolerance_deg = 0.75; % original 1.5, this one 0.75

if hdist_tolerance_deg <1
    data_set = 2;
else
    data_set = 1;
end

% find 
 figure, set(gcf,'renderer','Painters')
for t = 2:pa.nTrials %pa.trialNumber %full set, change to pa.nTrials
    tf = isbetween(synced, oc.UTCtrialStart(t), oc.UTCtrialEnd(t));
    idx = find(tf);
    xvals = [startpos(1), endpos(1)]; % get horizontal vals for start and end pos
    distX  = abs(x(tf) - xvals); % calculate distance between either
    if sum(tf)>0 && ~sum(find(distX>hdist_tolerance_deg))
        good_trials = [good_trials, t];
%         hold on, scatter(x(tf), y(tf))
        hold on, plot(x(tf), y(tf), '.-', 'LineWidth', 3, 'color',[0, 0, 0])
    elseif sum(tf)>0 && sum(find(distX>hdist_tolerance_deg))
        bad_trials = [bad_trials, t];
%         hold on, scatter(x(tf), y(tf), 50,[0.5, 0.5, 0.5])
%         hold on, plot(x(tf), y(tf), '.-', 'LineWidth', 3)

        hold on, plot(x(tf), y(tf),'color',[0.5, 0.5, 0.5], 'LineWidth', 3)

        % does the target go in the direction of the eye movement? does
        % this align with response?
        % what direction did the eye movement go 
        [val,val_idx] = max(abs(diff(x(tf)))); % peak horizontal velocity
        eye_direction = sign(x(idx(val_idx)) - x(idx(1)));
        target_comparison = [target_comparison pa.LR(t) == eye_direction];
        response_comparison = [response_comparison pa.LRresponse(t) == eye_direction];
    end
    
end

hold on, scatter(startpos(1), startpos(2), 50, 'filled', 'g')
hold on, scatter(endpos(1), endpos(2), 50, 'filled', 'r')
axis equal
set(gca, 'FontSize', 16)

disp(['num good trials = ', num2str(length(good_trials)), ' out of ', num2str(pa.nTrials), ', ' num2str(length(good_trials)/pa.nTrials*100), '%'])


pa.good_trials = good_trials;


% baseDir = pwd;
% pa.dataFile = fullfile(baseDir, 'Data', [pa.subjectName '-' num2str(pa.block) '-' pa.date '-eyetracking' '.mat']);
% save(pa.dataFile, 'pa', 'ds', 'kb','oc');
% 
%% try for multiple trials
% const_dev = zeros(1,10);
% surr_dev = zeros(1,10);

% for redo = 1:10
tic
load('theta.mat'); % theortical theta
dconst_overTrials = [];
dsurr_overTrials = [];
depth_range = 0.05;
depth_est = 0;

for t = good_trials %pa.trialNumber %full set, change to pa.nTrials
    tf = isbetween(synced, oc.UTCtrialStart(t), oc.UTCtrialEnd(t)+milliseconds(50)); %in case being inbetween trialstart and end cuts off too much eyetracking data
    trial_times = [];
    eyetracking = [];
    theta_thisTrial = [];
    trial_times = [trial_times; synced(tf)];
    eyetracking = [eyetracking; x(tf), y(tf)];
    idx = find(tf);
%     gaze_speed(t) = vecnorm([x(idx(1)) y(idx(1))]-[x(idx(end)) y(idx(end))])/.5;
%     start_position(:,t) = [x(idx(1)); y(idx(1))];
%     end_position (:,t) = [x(idx(end)); y(idx(end))];
%     hold on, scatter(x(tf), y(tf))
    
    interp_times = trial_times(1):milliseconds(1/ds.frameRate*1000):trial_times(end);
%         if length(theta)>length(interp_times)
%             
%             interp_times = trial_times(1):milliseconds(1/ds.frameRate*1000):trial_times(end)+100*milliseconds(1/ds.frameRate*1000);
%         end
    interp_x = interp1(trial_times, eyetracking(:,1), interp_times);
    interp_y = interp1(trial_times, eyetracking(:,2), interp_times);
    

    theta_thisTrial = theta(1)- deg2rad(interp_y(1:length(theta)));
    trial = t;
    
    [dconst, dsurr] = DistanceToConstraint(ds, pa, depth_range, depth_est, theta, trial);
    dconst_overTrials = [dconst_overTrials dconst];
    dsurr_overTrials = [dsurr_overTrials dsurr];
    if mod(t,50) == 0
        display(t)
    end
end
toc

%%
% save trials
old_goodTrials = pa.goodTrials;
pa.goodTrials = {};
pa.goodTrials{1} = old_goodTrials;
pa.goodTrials{data_set} = good_trials;

old_dconst = pa.data_const(:,1)';
pa.dconst = {};
pa.dconst{1} = old_dconst;
pa.dconst{data_set} = dconst_overTrials;

old_dsurr = pa.data_surr(:,1)';
pa.dsurr = {};
pa.dsurr{1} = old_dsurr;
pa.dsurr{data_set} = dsurr_overTrials;


% psignifit with new dconst
% matrix of results
pcorrect = eq(pa.LR(pa.goodTrials{data_set}), pa.LRresponse(pa.goodTrials{data_set}));
pcorrect = +pcorrect;
nTrials = ones(size(pcorrect));

% do this for both blocks
data_const = [pa.dconst{data_set}; pcorrect; nTrials]';
data_surr = [pa.dsurr{data_set}; pcorrect; nTrials]';

all = 1:pa.nTrials;
idx = ismember(all, good_trials);
removed_trials = all(~idx);

old_removed_trials = pa.removed_trials;
pa.removed_trials = {};
pa.removed_trials{1} = old_removed_trials;
pa.removed_trials{data_set} = removed_trials;

pa.data_const = data_const;
pa.data_surr = data_surr;
pa.removed_trials = removed_trials;
 
% save(pa.dataFile, '-struct', 'MPdata');
pa.baseDir = pwd;
pa.dataFile = fullfile(pa.baseDir, 'Data', [pa.subjectName '-' num2str(pa.block) '-' pa.date '-eyetracking' '.mat']);
save(pa.dataFile, 'pa', 'ds', 'kb','oc');

%% after running both blocks, clear and load the first block
% do this for first block
data_set = 2;
pcorrect = eq(pa.LR(pa.goodTrials{data_set}), pa.LRresponse(pa.goodTrials{data_set}));
pcorrect = +pcorrect;
nTrials = ones(size(pcorrect));

data_const = [pa.dconst{data_set}; pcorrect; nTrials]';
data_surr = [pa.dsurr{data_set}; pcorrect; nTrials]';

% load second block (don't delete anything) and do this
pcorrect = eq(pa.LR(pa.goodTrials{data_set}), pa.LRresponse(pa.goodTrials{data_set}));
pcorrect = +pcorrect;
nTrials = ones(size(pcorrect));

data_const = [data_const; [pa.dconst{data_set}; pcorrect; nTrials]'];
data_surr = [data_surr; [pa.dsurr{data_set}; pcorrect; nTrials]'];

%% then fit over both blocks
%
depth_range = 1.05;
depth_est = 0;

options             = struct;   % initialize as an empty struct
options.sigmoidName = 'weibull';   
options.expType     = '2AFC';   % choose 2-AFC as the paradigm of the experiment
                                % this sets the guessing rate to .5 and
                                % fits the rest of the parameters
options.fixedPars = NaN(5,1);  
options.poolxTol = 0.005;
% options.fixedPars(5) = 0;       % fix eta (dispersion) to zero
result_const = psignifit(data_const,options);
result_surr = psignifit(data_surr,options);


figure, plotPsych(result_const, options);
title(['depth range = ', num2str(depth_range)])
figure, plotPsych(result_surr, options);
title(['depth range = ', num2str(depth_range)])


display(['const dev = '  num2str(result_const.deviance)])
display(['surr dev = '  num2str(result_surr.deviance)])

% const_dev(redo) = result_const.deviance;
% surr_dev(redo) = result_surr.deviance;

% end

%% archive 
% result_const.depth_range = depth_range;
% result_surr.depth_range = depth_range;
% 
% JO_processed.full.result_surr = result_surr;
% JO_processed.full.result_const = result_const;
% 
% % save('JO_processed', 'JO_processed')
% 
% %% save figs
% s = 1;
% figFolder = '/Users/hopelutwak/Documents/GitHub/VRopticflow/Figures';
% stim = extractBefore(stims(1),'-');
% figname = [subjects(s)+'_const_eyetracking_'+stim+'.eps'];
% 
% 
% figname = [subjects(s)+'_surr_eyetracking_'+stim+'.eps'];
% 
% saveas(gcf, fullfile(figFolder, figname), 'epsc')
% 
% 
% %% archive
% %% look at trajectory for a trial (SKIP)
% for t= 4 %pa.trialNumber %full set, change to pa.nTrials
%     tf = isbetween(synced, oc.UTCtrialStart(t), oc.UTCtrialEnd(t));
%     trial_times = [trial_times; synced(tf)];
%     eyetracking = [eyetracking; x(tf), y(tf)];
%     idx = find(tf);
% 
% figure, for ii = 1:length(idx),hold on, scatter(x(idx(ii)), y(idx(ii)))
% pause(.05), end
% end
% 
% %% interpolate eye data one trial (SKIP)
% interp_times = trial_times(1):milliseconds(1/ds.frameRate*1000):trial_times(end);
% interp_x = interp1(trial_times, eyetracking(:,1), interp_times);
% interp_y = interp1(trial_times, eyetracking(:,2), interp_times);
% 
% figure, plot(trial_times,eyetracking(:,1),'o',interp_times,interp_x,':.');
% hold on, plot(trial_times,eyetracking(:,2),'o',interp_times,interp_y,':.');
% 
% % startpos = [5.4, 15.4]; %guess MP [5,-1.45]; DL [6,14], MG [5.4, 15.4]
% % endpos = [5.2, 11.2]; % DL [5,11], MG [5.2, 11.2]
% 
% % subtract off what start and end positions are
% interp_x = interp_x - startpos(1);
% interp_y = interp_y - startpos(2);
% 
% % for now try just changing theta for y
% load('theta.mat')
% theta = theta(1)- deg2rad(interp_y(1:length(theta)));
% 
% trial = 3;
% 
% % calculate for this trial
% 
% [dconst, dsurr] = DistanceToConstraint(ds, pa, .05, theta, trial);
% 
% 
% %%
% % JO_processed.full2.const = data_const;
% % JO_processed.full2.surr = data_surr;
% % JO_processed.full2.good_trials = good_trials;
% % JO_processed.full2.removed_trials = removed_trials;
% % 
% % save('JO_processed', 'JO_processed')
% 
% %% do this once you've gone through both blocks
% JO_processed.full.const = [JO_processed.full1.const; JO_processed.full2.const]
% JO_processed.full.surr = [JO_processed.full1.surr; JO_processed.full2.surr]
% 
% 
% JO_processed.dots.const = [JO_processed.dots1.const; JO_processed.dots2.const]
% JO_processed.dots.surr = [JO_processed.dots1.surr; JO_processed.dots2.surr]
% 
% JO_processed.monocular.const = [JO_processed.monocular1.const; JO_processed.monocular2.const]
% JO_processed.monocular.surr = [JO_processed.monocular1.surr; JO_processed.monocular2.surr]
% 
% 
% data_const = JO_processed.full.const;
% data_surr = JO_processed.full.surr;
% 
% data_const = JO_processed.dots.const;
% data_surr = JO_processed.dots.surr;
% 
% data_const = JO_processed.monocular.const;
% data_surr = JO_processed.monocular.surr;