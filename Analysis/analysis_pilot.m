%% analysis
addpath(genpath('Data'))

%% plotting HMD data
% load oc.modelViewDataLeft
% load('test-fixed-20230322T142931-0')
figure
plot3(oc.HMD(3:4:end,4), oc.HMD(1:4:end,4),...
    oc.HMD(2:4:end,4))
axis equal
xlabel('z')
ylabel('x')
zlabel('y')
title('HMD tracking')

figure
plot3(oc.modelViewDataLeft(3:4:end,4), oc.modelViewDataLeft(1:4:end,4),...
    oc.modelViewDataLeft(2:4:end,4))
axis equal
xlabel('z')
ylabel('x')
zlabel('y')
title('camera view tracking')


% figure, plot3(pa.zPosition, pa.xPosition,ones(size(pa.xPosition))*pa.floorHeight+pa.paddleHalfHeight)
% axis equal
% xlabel('z')
% ylabel('x')
% zlabel('y')

% how much did you travel
vecnorm(oc.modelViewDataLeft(end-3:end-1,4)-oc.modelViewDataLeft(1:3,4))

%% 3D velocities tested


% c = [255,153,153; 255,51,51; 204,0,0; 153,0,0;
%     153,255,153; 0,255,0; 0,204,0; 0,102,0;
%     153,204,255; 51,153,255; 0,128,255; 0,76,153;
%     204,153,255; 178,102,255; 153,51,255; 102,0,204]/255;
% c = [255,153,153; 255,51,51; 204,0,0; 
%                     153,255,153; 0,255,0; 0,204,0;
%                     153,204,255; 51,153,255; 0,128,255; 
%                     204,153,255; 178,102,255; 153,51,255;
%                     255, 204, 153; 255,153,51; 204,102,0]/255;
c = options.dataColor;
%            
factorial = fullfact([length(pa.speed), length(pa.direction)]);
fullfacts = [(pa.speed(factorial(:,1)).*cos(pa.direction(factorial(:,2)))); -(pa.speed(factorial(:,1)).*sin(pa.direction(factorial(:,2))))];

figure, hold on,scatter(fullfacts(1,:), -fullfacts(2,:),[],c,'filled'), axis equal
xlim([-max(pa.speed)*1.2, max(pa.speed)*1.2])
ylim([-max(pa.speed)*1.2, max(pa.speed)*1.2])


%% screencap
s = size(screenCap);

figure
set(gcf,'position',[0, 0, 1280, 1600])

for ii = 1:s(3)
    imagesc(screenCap(:,:,ii)), colormap(gray)
    axis equal
    pause(1/10)
end

%% eye tracking data
D=dir('Data/');
% filename = '2023-10-11_14-57-35-1e54f077';
% filename = '2023-10-12_12-25-27-004bb3c5'; % fix/obj dist 2,1.5
% filename = '2023-11-02_14-35-26-21023a71ET';
%filename = '2024-02-13_17-47-47-d08cbe84HLcali2';
% filename = '2024-02-21_12-36-02_hl_fulleyetracking-fixed-15trans-20240221t123559-0-012e0633';
% filename = '2024-02-21_13-15-23_hl_fulleyetrackingdots-fixed-20240221t131520-0-44d050a3';
% filename = '2024-02-21_14-44-58_hl_fulleyetracking_eyesim-fixed-20240221t144454-0-f2d8da04';
filename = '2024-02-26_15-15-38_MPfull1';

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

% % in UTC time
% figure, plot(timestamps,x, 'linewidth', 2), hold on, plot(timestamps,y,'linewidth', 2)
% yl = ylim;
% exp_timestamps = posixtime(oc.UTCtrialStart)*1e9;
% hold on, line([exp_timestamps; exp_timestamps],[yl(1); yl(2)].*ones(size(oc.UTCtrialStart)), 'color','k') 

% in datetime
date_time = datetime(timestamps/1e9, 'ConvertFrom', 'posixtime', 'TimeZone','local', 'Format', 'd-MMM-y HH:mm:ss:ms');
blink_start_time = datetime(blink_start/1e9, 'ConvertFrom', 'posixtime', 'TimeZone','local', 'Format', 'd-MMM-y HH:mm:ss:ms');
blink_end_time= datetime(blink_end/1e9, 'ConvertFrom', 'posixtime', 'TimeZone','local', 'Format', 'd-MMM-y HH:mm:ss:ms');
evts_time = datetime(evt/1e9, 'ConvertFrom', 'posixtime', 'TimeZone','local', 'Format', 'd-MMM-y HH:mm:ss:ms');


figure, plot(date_time,x, 'linewidth', 2), hold on, plot(date_time,y,'linewidth', 2)
yl = ylim;
% hold on, line([oc.UTCtrialStart; oc.UTCtrialStart],[yl(1); yl(2)].*ones(size(oc.UTCtrialStart)), 'color','g') 
% hold on, line([oc.UTCtrialStart; oc.UTCtrialStart]+seconds(.5),[yl(1); yl(2)].*ones(size(oc.UTCtrialStart)), 'color','r') 

% hold on, line([oc.UTCtrialEnd; oc.UTCtrialEnd]-seconds(pa.targetMotionDuration),[yl(1); yl(2)].*ones(size(oc.UTCtrialEnd)), 'color','g') 
hold on, line([oc.UTCtrialEnd; oc.UTCtrialEnd],[yl(1); yl(2)].*ones(size(oc.UTCtrialEnd)), 'color','r') 

% hold on, line([blink_start_time'; blink_start_time'], [yl(1); yl(2)].*ones(size(blink_start_time')),'color','#77AC30');
% hold on, line([blink_end_time'; blink_end_time'], [yl(1); yl(2)].*ones(size(blink_end_time')),'color','#A2142F');

% time in between start and end
% milliseconds(oc.UTCtrialEnd - oc.UTCtrialStart)
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
[b,a]=butter(4,1/25);
output_datax=filter(b,a,x);
output_datay = filter(b,a,y);

figure, plot(synced,output_datax, 'linewidth', 2), hold on, plot(synced,output_datay,'linewidth', 2)
yl = ylim;
hold on, line([oc.UTCtrialEnd; oc.UTCtrialEnd],[yl(1); yl(2)].*ones(size(oc.UTCtrialEnd)), 'color','r') 

hold on, line([oc.UTCtrialStart; oc.UTCtrialStart],[yl(1); yl(2)].*ones(size(oc.UTCtrialStart)), 'color','g') 
hold on, line([evts_time'; evts_time'], [yl(1); yl(2)].*ones(size(evts_time')), 'color','k') 

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
    [v,ii] = max(y(idx));
    start_position(:,t) = [x(idx(ii)); y(idx(ii))];
    end_position (:,t) = [x(idx(end)); y(idx(end))];
    hold on, scatter(x(tf), y(tf))
end

% figure, hold on, scatter(eyetracking(:,1), eyetracking(:,2))
% axis equal

% figure
% xlim([700, 1000])
% ylim([600, 800])
% axis equal
% 
% for t = 1:length(eyetracking)-100
%     plot(eyetracking(t:t+100,1), eyetracking(t:t+100,2))
%        pause(.01)
% end

% remove trials with high gaze speed
reasonable_speed = find(gaze_speed<10);

% average start and end position
startpos = nanmean(start_position(:,reasonable_speed)');
endpos = nanmean(end_position(:,reasonable_speed)');
hold on, scatter(startpos(1), startpos(2), 'filled', 'g')
hold on, scatter(endpos(1), endpos(2), 'filled', 'r')

%% getting trials that are outside degree radius of fixation
% load('track_theta')
% expected_diff = max(rad2deg(track_theta)) - min(rad2deg(track_theta));
% expected_diff - vecnorm(startpos - endpos)
% looks reasonable compared to simulation?

% go through each trial again, sample line that goes from average start
% point and endpoint with as many timepoints as that trial, for each
% timestamp if eye tracking goes outside of degree radius of the line, flag
% the trial

trial_times = [];
eyetracking = [];
dstart = [];
dend = [];
% start_position = NaN(2, pa.trialNumber-1);
% end_position = NaN(2, pa.trialNumber-1);

bad_trials = [];

% avg_trials = find(gaze_speed>5 & gaze_speed<6);
for trial = 2:pa.nTrials %600 % 2:pa.nTrials %full set, change to pa.nTrials
    tf = isbetween(synced, oc.UTCtrialStart(trial), oc.UTCtrialEnd(trial));
   
    trial_times = synced(tf);
    eyetracking = [x(tf), y(tf)];
    idx = find(tf);
     
%     averageEyePath=[linspace(startpos(1),endpos(1),length(eyetracking)).' linspace(startpos(2),endpos(2),length(eyetracking)).'];
    [peak,ii] = max(eyetracking(:,2));
    dstart = [dstart vecnorm(eyetracking(ii,:) - startpos)];
    dend = [dend vecnorm(eyetracking(end,:)-endpos)];
    if dstart(end) >1 || dend(end) >1
        bad_trials = [bad_trials, trial];
    end
    
%     for t = 1:length(eyetracking)
%         d = vecnorm(eyetracking(t,:) - averageEyePath(t,:));
%         if d>1
%             bad_trials = [bad_trials, trial];
%             break
%         end
%     end
end



figure, hold on, plot(eyetracking(:,1), eyetracking(:,2), 'LineWidth', 5), axis equal
% hold on, scatter(averageEyePath(:,1), averageEyePath(:,2))

% hold on, scatter(eyetracking(t,1), eyetracking(t,2))
% hold on, scatter(averageEyePath(t,1), averageEyePath(t,2))
% legend('eyetracking', 'average eye path', 'first break eyetracking', 'first break eye path')

eye_pos_trial = [x(tf) y(tf)];
figure, histogram(dstart, 'BinWidth', .5);
hold on, histogram(dend, 'BinWidth', .5);

%% get fixation over trial intervals

start_times = timetable(oc.UTCtrialStart', ones(numel(oc.UTCtrialStart),1));
end_times = timetable([oc.UTCtrialStart+seconds(0.5)]', ones(numel(oc.UTCtrialStart),1));

eyetracking = timetable(date_time, x, y);
TT = synchronize(start_times, end_times, eyetracking);
synched = table2array(TT);
idx_start = find(synched(:,1) == 1);
idx_end = find(synched(:,2) ==1);

stim_interval_idx = [];
n_practice = 1;
for t = n_practice+1:pa.nTrials %pa.nTrials %cut out practice trials
    stim_interval_idx = [stim_interval_idx idx_start(t)+1:idx_end(t)-1];
end
figure, scatter(synched(stim_interval_idx, 3), synched(stim_interval_idx, 4))
figure, plot(TT.Time(stim_interval_idx), TT.x(stim_interval_idx))
hold on, plot(TT.Time(stim_interval_idx), TT.y(stim_interval_idx))
hold on, line([oc.UTCtrialStart; oc.UTCtrialStart],[yl(1); yl(2)].*ones(size(oc.UTCtrialStart)), 'color','g') 
hold on, line([oc.UTCtrialStart; oc.UTCtrialStart]+seconds(.5),[yl(1); yl(2)].*ones(size(oc.UTCtrialStart)), 'color','r') 

%check variablility of fixation over time
% figure, scatter(x_trials, y_trials)
% interval_duration = (timestamps(trial_interval(end))-timestamps(trial_interval(1)))/1e9;
cov(TT.x(stim_interval_idx), TT.y(stim_interval_idx))

%% load behavioral file
load('Data\HL_pilot-fixed-20230726T160928-0.mat')

% have to take data from stimulus with eye simulation
change_angle = rad2deg(max(track_theta) - min(track_theta));

% fixations
fixations = readtable('Data/2023-07-26_16-09-57-7a8b312d/fixations.csv');

% event check
events = readtable('Data/2023-08-30_16-12-45-6bae8328/events.csv');
ev_timestamps = table2array(events(:,2));

t_ne = uint64(1582650648869329937);
NS = 1e9;
right_over = mod(t_ne, NS);
left_over = t_ne - right_over;
d = datetime( double(left_over)/NS, 'convertfrom', 'posixtime', 'Format', 'dd-MMM-uuuu HH:mm:ss.SSSSSSSSS') + seconds(double(right_over)/NS)


%% response data
% 90 is forward, 270 is backward
pcorrect = NaN(length(pa.speed),length(pa.direction));
data = NaN(length(pa.speed), 3, length(pa.direction));
for speed = 1:length(pa.speed)
    idx_speed = find(pa.fullFactorial(:,3) == pa.speed(speed));
    for direction = 1:length(pa.direction)
        idx_dir = find(pa.fullFactorial(:,4) == pa.direction(direction));
        idx = intersect(idx_speed, idx_dir);
        pcorrect(speed,direction) = sum(eq(pa.LR(idx), pa.LRresponse(idx)))/pa.nRepeats;
        data(speed, :, direction) = [pa.speed(speed), sum(eq(pa.LR(idx), pa.LRresponse(idx))), pa.nRepeats];
    end
end

%% psignifit 
% addpath(genpath('C:\Users\hlutw\OneDrive\Documents\MATLAB\psignifit-master'));
 addpath('/Users/hopelutwak/Documents/MATLAB/psignifit')
addpath(genpath('/Users/hopelutwak/Documents/GitHub/VRopticflow/Analysis'));
% addpath(genpath('C:\Users\hlutw\OneDrive\Documents\GitHub\VRopticflow\Analysis'));
% set options for psychometric functions
options             = struct;   % initialize as an empty struct
options.sigmoidName = 'weibull';   
options.expType     = '2AFC';   % choose 2-AFC as the paradigm of the experiment
                                % this sets the guessing rate to .5 and
                                % fits the rest of the parameters
options.fixedPars = NaN(5,1);                                
options.fixedPars(5) = 0;       % fix eta (dispersion) to zero
options.dataColor = [255,153,255; 255,102,255; 255,51,255; 204,0,204;
                    255,153,153; 255,102,102; 255,51,51; 204,0,0;
                    255,204,153; 255,178,102; 255, 153, 51; 204,102,0;
                    204,255,153; 178,255,102; 153,255,51; 102,204,0;
                    153,255,255; 102,255,255; 51,255,255; 0,204,204;
                    153,153,255; 102,102,255; 51,51,255; 0,0,204]/255;
                
figure, hold on

%loop through stim conditions and get threshold and plot curves
n_conditions = length(pa.direction);
thresholds = zeros(1,n_conditions);
CI = zeros(2, n_conditions);
for ii = 1:n_conditions
    result = psignifit(data(:,:,ii),options);
    thresholds(ii)= exp(result.Fit(1));
    CI(:,ii) =exp(result.conf_Intervals(1,:,1))'; %get CI for threshold (first row) at 95% (first layer)

    subplot(1,n_conditions,ii)
    hold on
    plotPsych(result);
    ylim([0 1])
    title(num2str(pa.direction(ii)))
%     sgtitle(stim) 
end

% reshape matrix so it's in terms of distance to constraint
C=permute(data,[1 3 2]);
C = reshape(C,[],size(data,2),1);
% relace speed column with distance to constraint variable
% C(:,1) = [0.005; 0.0026; 0.0103; 0.0051; 0.0092; 00.0046];

% obj 3m, depth range .2
% a = [ 0.0135     0.0710    0.0167   0.0840; 0.0027     0.0131    0.0027   0.0132; 2e-5       0.0052    0.0012   0.0036;    2e-5  0.0024    0.0004   0.0014];

% obj 3m, depth range .1m
% a = [0.0145   0.0495   0.0112   0.0559; 0.0040   0.0120   0.0022   0.0129; 0.0013   0.0053   0.0010   0.0055; 0.0001   0.0021   0.0004   0.0020];

% obj 3m, depth range .2?
% a = [0.0044 0.0047; 0.0021 0.0022; 0.0010 0.0011];

% obj 2m, depth range .2
% a = [0.0300    0.0792    0.0234    0.0928;  0.0029    0.0153    0.0038    0.0124;  2e-5     0.0069    0.0018    0.0038;  2e-5     0.0033    0.0008    0.0018];

% obj 2m, depth range .2, direction  250  270  290  350, speeds 0.1500    0.0750    0.0375    0.0187
% a = [0.0064    0.0007    0.0102    0.0229; 0.0030   6.2774e-05  0.0031     0.0090;  0.0014     6.3316e-05    0.0015     0.0041; 0.0006     6.3573e-05    0.0007     0.0021];
 % - distance to constraint point
%  a = [ 0.0080    0.0160    0.0285  0.0382;  0.0038    0.0080    0.0135  0.0191;  0.0019    0.0041    0.0068  0.0097; 0.0011    0.0022    0.0035  0.0050];
% speeds Speeds 0.5  0.25  .125    0.0625  m/s, .2m depth estimate
% a = [0.0250    0.0292    0.0374    0.0768    0.1020; 0.0102    0.0116    0.0066    0.0265    0.0457; 0.0008    0.0052    0.0002    0.0073    0.0178; 6.4579e-05      0.0025    6.2959e-05     0.0025    0.0072];
%.1 distance
% a = [    0.0281    0.0315    0.0576    0.0901    0.1100;
%     0.0141    0.0116    0.0196    0.0401    0.0551;
%     0.0049    0.0052    0.0048    0.0144    0.0244;
%     0.0002    0.0025    0.0002    0.0041    0.0097];
% corrected height
%to constraint
% a =    [ 0.0315    0.0405    0.0567    0.0883    0.1011;
%     0.0101    0.0065    0.0090    0.0187    0.0233;
%     0.0051    0.0032    0.0044    0.0090    0.0112;
%     0.0024    0.0016    0.0023    0.0046    0.0057];

% a = [    0.0052    0.0026    0.0048    0.0097    0.0144;
%     0.0025    0.0012    0.0002    0.0021    0.0041;
%     0.0012    0.0006    0.0001    0.0007    0.0013];
% distance to constraint
% a = [    0.0065    0.0090    0.0137    0.0187    0.0233;
%     0.0032    0.0044    0.0067    0.0090    0.0112;
%     0.0016    0.0023    0.0034    0.0046    0.0057];
% distance .2m
% a = [0.0052    0.0026    0.0048    0.0097    0.0144;
%     0.0025    0.0012    0.0002    0.0021    0.0041;
%     0.0012    0.0006    0.0001    0.0007    0.0013];


% 
% a = [    0.0281    0.0315    0.0576    0.0901    0.1100;
%     0.0141    0.0116    0.0196    0.0401    0.0551;
%     0.0049    0.0052    0.0048    0.0144    0.0244;
%     0.0002    0.0025    0.0002    0.0041    0.0097;
%     0.0052    0.0026    0.0048    0.0097    0.0144;
%     0.0025    0.0012    0.0002    0.0021    0.0041;
%     0.0012    0.0006    0.0001    0.0007    0.0013];
% options.dataColor = [255,153,255; 255,102,255; 255,51,255; 204,0,204;
%                     255,153,153; 255,102,102; 255,51,51; 204,0,0;
%                     255,204,153; 255,178,102; 255, 153, 51; 204,102,0;
%                     204,255,153; 178,255,102; 153,255,51; 102,204,0;
%                     153,255,255; 102,255,255; 51,255,255; 0,204,204;
%                     153,153,255; 102,102,255; 51,51,255; 0,0,204]/255;
[dconst, dsurr] = DistanceToConstraint(ds, pa, .1);
a = dsurr;
C(:,1) = a(:);

% run psignifit
result = psignifit(C,options);
figure, plotPsych(result, options);
thresh = exp(result.Fit(1));
dev = result.deviance;
% options.dataColor = repmat([0,0,1], length(C),1);
% ** will only work of edit psignifit's plotPsych!!

% options.dataColor = [255,153,153; 255,51,51; 204,0,0; 153,0,0; 
%                      153,255,153; 0,255,0; 0,204,0; 0,102,0;
%                      153,204,255; 51,153,255; 0,128,255; 0,76,153;
%                      204,153,255; 178,102,255; 153,51,255; 102,0,204]/255;
% colors = []
% for ii = 1:numel(pa.direction)
%     c = 
% 

% 4 speeds, 6 directions
options.dataColor = [255,153,255; 255,102,255; 255,51,255; 204,0,204;
                    255,153,153; 255,102,102; 255,51,51; 204,0,0;
                    255,204,153; 255,178,102; 255, 153, 51; 204,102,0;
                    204,255,153; 178,255,102; 153,255,51; 102,204,0;
                    153,255,255; 102,255,255; 51,255,255; 0,204,204;
                    153,153,255; 102,102,255; 51,51,255; 0,0,204]/255;
 % 5 speeds, 6 directions 
%                 options.dataColor = [255,204,255; 255,153,255; 255,102,255; 255,51,255; 204,0,204;
%                     255,204,204; 255,153,153; 255,102,102; 255,51,51; 204,0,0;
%                     255,204,153; 255,204,153; 255,178,102; 255, 153, 51; 204,102,0;
%                     229,255,204; 204,255,153; 178,255,102; 153,255,51; 102,204,0;
%                     205,255,255; 153,255,255; 102,255,255; 51,255,255; 0,204,204;
%                     204,229,255; 153,153,255; 102,102,255; 51,51,255; 0,0,204]/255;
% options.dataColor = [224,224,224; 160,160,160; 128,128,128; 64,64,64;
%                       255,153,153; 255,102,102; 255,51,51; 204,0,0;
%                       255,255,153; 255,255,51; 204,204,0; 153,153,0;
%                       153,204,255; 51,153,255; 0,128,255; 0,102,204;
%                       204,153,255; 178,102,255; 153,51,255; 102,0,204;
%                       255,51,51; 204,0,0; 153,0,0;
%                       255,153,51; 255,128,0; 204,102,0;
%                       204,204,0; 153,153,0; 102,102,0;
%                       153,255,153; 0,255,0; 0,204,0;
%                       0,128,255; 0,102,204; 0,76,153]/255;
% options.dataColor = jet(size(result.data,1));

% 90  250     270     290  350
%     250 260 270 280 290
%  0.5  0.25  .125    0.0625
%            0.1250   0.0625    0.0313

%% iterate over different values of distance to const
distances  = linspace(.025, .1, 5);
dev = zeros(1,length(distances));
for d = 1:length(distances)
    [dconst, dsurr] = DistanceToConstraint(ds, pa, distances(d));
    a = dconst;
    C(:,1) = a(:);

    % run psignifit
    result = psignifit(C,options);
    figure, plotPsych(result, options);
%     thresh = exp(result.Fit(1));
    dev(d) = result.deviance;
end


%% using psychtoolbox

[alpha,beta,thresh92] = FitWeibTAFC(C(:,1),C(:,2),C(:,3)-C(:,2))

x  = exp(linspace(log(min(result.data(:,1))),log(max(result.data(:,1))),1000));
y = 1 - exp(-(x./beta)).^alpha;

figure, plot(x, wblcdf(x, alpha, beta))
hold on, scatter(C(:,1), C(:,2)./C(:,3))

% qq plot
load lightbulb
lightbulb = [lightbulb(lightbulb(:,3)==0,1),...
    lightbulb(lightbulb(:,3)==0,2)];
fluo = [lightbulb(lightbulb(:,2)==0,1)];
insc = [lightbulb(lightbulb(:,2)==1,1)];
pd = makedist('Weibull');
figure
qqplot(fluo,pd)

%example weibull
x = linspace(0,30);
plot(x,wblcdf(x,10,1),'DisplayName','A=10, B=1')
hold on
plot(x,wblcdf(x,10,2),'DisplayName','A=10, B=2')
plot(x,wblcdf(x,10,4),'DisplayName','A=10, B=4')
plot(x,wblcdf(x,10,0.5),'DisplayName','A=10, B=0.5')
plot(x,wblcdf(x,10,0.25),'DisplayName','A=10, B=0.25')
plot(x,wblcdf(x,5,1),'DisplayName','A=5, B=1')
hold off
legend('show','Location','southeast')
xlabel('x')
ylabel('cdf')


% %function y = Weibull(p,x)
% %y = Weibull(p,x)
% %
% %Parameters:  p.b slope
% %             p.t threshold yeilding ~80% correct
% %             x   intensity values.
% 
% g = 0.5;  %chance performance
% e = (.5)^(1/3);  %threshold performance ( ~80%)
% 
% %here it is.
% k = (-log( (1-e)/(1-g)))^(1/p.b);
% y = 1- (1-g)*exp(- (k*x/p.t).^p.b);
%% random
% exp2 = data;
% exp1 = data;

% C2 = C;
% C1 = C;
% C = [C1;C2];
