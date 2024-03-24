%% incorporating eye movement from eye tracking into distance to constraint calculation

%% get eye tracking interval
addpath('/Users/hopelutwak/Documents/MATLAB/psignifit')
addpath(genpath('/Users/hopelutwak/Documents/GitHub/VRopticflow/Analysis'))
dataFolder = '/Users/hopelutwak/Documents/GitHub/VRopticflow/Data';


% names of files
S = dir(fullfile(dataFolder));
% 
% % which subjects data to analyze
% subjects = ["DL"]; %"HL" "IK"
% stims = ["full-1"]; %["full-1", "full-2"]; %"pilot"

D=dir('Data/');
filename = '2024-03-06_11-52-42_MG-full-1';

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


%% offset
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


%% butterworth filter
[b,a]=butter(4,1/50);
output_datax=filter(b,a,x);
output_datay = filter(b,a,y);

figure, plot(synced,output_datax, 'linewidth', 2), hold on, plot(synced,output_datay,'linewidth', 2)
yl = ylim;
hold on, line([oc.UTCtrialEnd; oc.UTCtrialEnd],[yl(1); yl(2)].*ones(size(oc.UTCtrialEnd)), 'color','r') 

hold on, line([oc.UTCtrialStart; oc.UTCtrialStart],[yl(1); yl(2)].*ones(size(oc.UTCtrialStart)), 'color','g') 
hold on, line([synced_evts'; synced_evts'], [yl(1); yl(2)].*ones(size(evts_time')), 'color','k') 

x = output_datax;
y= output_datay;


%% eye data within intervals
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
    gaze_speed(t) = vecnorm([x(idx(1)) y(idx(1))]-[x(idx(end)) y(idx(end))])/.5;
    start_position(:,t) = [x(idx(1)); y(idx(1))];
    end_position (:,t) = [x(idx(end)); y(idx(end))];
    hold on, scatter(x(tf), y(tf))
    
end

%% find good trials, no large horizontal variations

good_trials = [];
startpos = [5.4, 15.4]; %guess MP [5,-1.45]; DL [6,14] [8,12], MG [5.4, 15.4]
% endpos = [5,11]; % DL [5,11], MG [5.2, 11.2]

% figure
for t = 2:pa.nTrials %pa.trialNumber %full set, change to pa.nTrials
    tf = isbetween(synced, oc.UTCtrialStart(t), oc.UTCtrialEnd(t));
    
    distX  = abs(x(tf) - startpos(1));
    if ~sum(find(distX>2))
        good_trials = [good_trials, t];
        hold on, scatter(x(tf), y(tf))
    end
    
end

set(gca, 'FontSize', 16)

%% interpolate eye data
interp_times = trial_times(1):milliseconds(1/ds.frameRate*1000):trial_times(end);
interp_x = interp1(trial_times, eyetracking(:,1), interp_times);
interp_y = interp1(trial_times, eyetracking(:,2), interp_times);

figure, plot(trial_times,eyetracking(:,1),'o',interp_times,interp_x,':.');
hold on, plot(trial_times,eyetracking(:,2),'o',interp_times,interp_y,':.');

startpos = [5.4, 15.4]; %guess MP [5,-1.45]; DL [6,14], MG [5.4, 15.4]
endpos = [5.2, 11.2]; % DL [5,11], MG [5.2, 11.2]

% subtract off what start and end positions are
interp_x = interp_x - startpos(1);
interp_y = interp_y - startpos(2);

% for now try just changing theta for y
load('theta.mat')
theta = theta(1)- deg2rad(interp_y(1:length(theta)));

trial = 4;

% calculate for this trial

[dconst, dsurr] = DistanceToConstraint(ds, pa, .05, theta, trial);


%% try for multiple trials
tic
load('theta.mat');
dconst_overTrials = [];
dsurr_overTrials = [];
subsetTrials = good_trials;
for t = subsetTrials %pa.trialNumber %full set, change to pa.nTrials
    tf = isbetween(synced, oc.UTCtrialStart(t), oc.UTCtrialEnd(t)+milliseconds(50)); %in case being inbetween trialstart and end cuts off too much eyetracking data
    trial_times = [];
    eyetracking = [];
    theta_thisTrial = [];
    trial_times = [trial_times; synced(tf)];
    eyetracking = [eyetracking; x(tf), y(tf)];
    idx = find(tf);
    gaze_speed(t) = vecnorm([x(idx(1)) y(idx(1))]-[x(idx(end)) y(idx(end))])/.5;
    start_position(:,t) = [x(idx(1)); y(idx(1))];
    end_position (:,t) = [x(idx(end)); y(idx(end))];
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
    
    [dconst, dsurr] = DistanceToConstraint(ds, pa, .05, theta, trial);
    dconst_overTrials = [dconst_overTrials dconst];
    dsurr_overTrials = [dsurr_overTrials dsurr];
    if mod(t,50) == 0
        display(t)
    end
end
toc


%% figure out how to save extra field to pa with dconst/surr values
% MPdata = load(pa.dataFile);
pa.goodTrials = good_trials;
pa.dconst = dconst_overTrials;
pa.dsurr = dsurr_overTrials;


% save(pa.dataFile, '-struct', 'MPdata');
pa.baseDir = pwd;
pa.dataFile = fullfile(pa.baseDir, 'Data', [pa.subjectName '-' num2str(pa.block) '-' pa.date '-copy' '.mat']);
save(pa.dataFile, 'pa', 'ds', 'kb','oc');

%% psignifit with new dconst

% matrix of results
pcorrect = eq(pa.LR(pa.goodTrials), pa.LRresponse(pa.goodTrials));
pcorrect = +pcorrect;
nTrials = ones(size(pcorrect));

options             = struct;   % initialize as an empty struct
options.sigmoidName = 'weibull';   
options.expType     = '2AFC';   % choose 2-AFC as the paradigm of the experiment
                                % this sets the guessing rate to .5 and
                                % fits the rest of the parameters
options.fixedPars = NaN(5,1);                                
% options.fixedPars(5) = 0;       % fix eta (dispersion) to zero

% do this for both blocks
data_const = [pa.dconst; pcorrect; nTrials]';
data_surr = [pa.dsurr; pcorrect; nTrials]';

%
% data_const2 = [pa.dconst; pcorrect; nTrials]';
% data_surr2 = [pa.dsurr; pcorrect; nTrials]';
% 
% data_const = [data_const; data_const2];
% data_surr = [data_surr; data_surr2];

result_const = psignifit(data_const,options);
result_surr = psignifit(data_surr,options);

figure, plotPsych(result_const, options);
figure, plotPsych(result_surr, options);