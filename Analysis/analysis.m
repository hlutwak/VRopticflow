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


c = [255,153,153; 255,51,51; 204,0,0; 153,0,0;
    153,255,153; 0,255,0; 0,204,0; 0,102,0;
    153,204,255; 51,153,255; 0,128,255; 0,76,153;
    204,153,255; 178,102,255; 153,51,255; 102,0,204]/255;
           
factorial = fullfact([length(pa.speed), length(pa.direction)]);
fullfacts = [(pa.speed(factorial(:,1)).*cos(pa.direction(factorial(:,2)))); -(pa.speed(factorial(:,1)).*sin(pa.direction(factorial(:,2))))];

figure, hold on,scatter(fullfacts(1,:), -fullfacts(2,:),'filled'), axis equal
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
D=dir('Data');
% gaze = readtable('Data\2023-07-26_16-09-57-7a8b312d\gaze.csv');
% worldtime = readtable('Data\2023-07-26_16-09-57-7a8b312d\world_timestamps.csv');
gaze = readtable('Data/2023-08-30_16-12-45-6bae8328/gaze.csv');
% t= table2array(worldtime(:,end));
timestamps = table2array(gaze(:,3));
x = table2array(gaze(:,4));
y = table2array(gaze(:,5));

figure, scatter(x,y), axis equal

figure, plot(timestamps,x, 'linewidth', 2), hold on, plot(timestamps,y,'linewidth', 2)

trial_interval = [2636:15260];
x_trials = x(trial_interval);
y_trials = y(trial_interval);

%check variablility of fixation over time
figure, scatter(x_trials, y_trials)
interval_duration = (timestamps(trial_interval(end))-timestamps(trial_interval(1)))/1e9
cov(x_trials, y_trials)

% load behavioral file
load('Data\HL_pilot-fixed-20230726T160928-0.mat')

% have to take data from stimulus with eye simulation
change_angle = rad2deg(max(track_theta) - min(track_theta));

% fixations
fixations = readtable('Data/2023-07-26_16-09-57-7a8b312d/fixations.csv');

% event check
events = readtable('Data/2023-08-30_16-12-45-6bae8328/events.csv');
ev_timestamps = table2array(events(:,2));

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
addpath(genpath('C:\Users\hlutw\OneDrive\Documents\MATLAB\psignifit-master'));
%  addpath('/Users/hopelutwak/Documents/MATLAB/psignifit')
% set options for psychometric functions
options             = struct;   % initialize as an empty struct
options.sigmoidName = 'weibull';   
options.expType     = '2AFC';   % choose 2-AFC as the paradigm of the experiment
                                % this sets the guessing rate to .5 and
                                % fits the rest of the parameters
options.fixedPars = NaN(5,1);                                
options.fixedPars(5) = 0;       % fix eta (dispersion) to zero

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

% a = [0.0300    0.0792    0.0234    0.0928;  0.0029    0.0153    0.0038    0.0124;  2e-5     0.0069    0.0018    0.0038;  2e-5     0.0033    0.0008    0.0018];
a = [0.0145   0.0495   0.0112   0.0559; 0.0040   0.0120   0.0022   0.0129; 0.0013   0.0053   0.0010   0.0055; 0.0001   0.0021   0.0004   0.0020];

% a = [0.0300    0.0792    0.0234    0.0928;  0.0029    0.0153    0.0038    0.0124;  2e-5     0.0069    0.0018    0.0038;  2e-5     0.0033    0.0008    0.0018];

% a = [0.0044 0.0047; 0.0021 0.0022; 0.0010 0.0011];
C(:,1) = a(:);
% 
% exp2 = data;
% exp1 = data;

% C2 = C;
% C1 = C;
% C = [C1;C2];

% run psignifit
result = psignifit(C,options);
% options.dataColor = repmat([0,0,1], length(C),1);
% ** will only work of edit psignifit's plotPsych!!

% options.dataColor = [255,153,153; 255,51,51; 204,0,0; 153,0,0; 
%                      153,255,153; 0,255,0; 0,204,0; 0,102,0;
%                      153,204,255; 51,153,255; 0,128,255; 0,76,153;
%                      204,153,255; 178,102,255; 153,51,255; 102,0,204]/255;
%  
figure, plotPsych(result, options);

