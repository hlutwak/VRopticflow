function [ds,pa] = SetupParameters(ds)

% Establishes the parameters for the Oculus 3D Pong experiment using the
% 2nd generation Oculus (DK2)

% At the start of the session, should only update lines 11 and 12


%% Basic experimental specs 

pa.subjectName = 'HL_pilot';
pa.feedbackFlag = 1;  % 0 --> no feedback;  1 --> only sound-based feedback;  2 --> visual + sound-based feedback

pa.criterion = 8*(pi/180);  % this is the criterion for the distance between the mid-point of the paddle and the midpoint of the ball; distances smaller than this are hits - any overlap counts

% Sound parameters
[yMiss, pa.freqMiss] = psychwavread('Sounds/Swish.wav');
pa.waveDataMiss = yMiss';
pa.nrchannelsMiss = size(pa.waveDataMiss,1); % Number of rows == number of channels.

[yHit, pa.freqHit] = psychwavread('Sounds/cowbell.wav');
pa.waveDataHit = yHit';
pa.nrchannelsHit = size(pa.waveDataHit,1); % Number of rows == number of channels.

% Setup sounds
InitializePsychSound(1);
pa.handleHit = PsychPortAudio('Open', [], [], 0, pa.freqHit, pa.nrchannelsHit);
PsychPortAudio('FillBuffer', pa.handleHit, pa.waveDataHit);
pa.handleMiss = PsychPortAudio('Open', [], [], 0, pa.freqMiss, pa.nrchannelsMiss);
PsychPortAudio('FillBuffer', pa.handleMiss, pa.waveDataMiss);



pa.date = datestr(now,30);
s = RandStream('mt19937ar','Seed','shuffle'); % this sets up a new state
RandStream.setGlobalStream(s); % sets the stream according to the state
defaultStream = RandStream.getGlobalStream; % for 64-bit; getCurrentStream for old 32-bit
pa.savedState = defaultStream.State; %% these two lines allow us to recover the subject's experience
pa.baseDir = pwd;


%% parameters for the texture surround and virtual 'room'
pa.translation = 1; %m/s forward


pa.fixationDiskRadius = 2.5; % deg
pa.fixationHalfSquareSize = 1;  % deg - just the 'box' in which the fixation lines will be drawn
pa.fixationHalfSquare = 0.0157; % m 
pa.fixationLineLength = 0.0079; % m 
    
pa.apertureRadius = 7.5;  % deg
   
pa.floorHeight = -1; % m
pa.floorWidth = 6;
pa.ceilingHeight = 1.5; % m 

pa.gazeangle = atan(-pa.floorHeight/(pa.floorWidth/2)); %angle camera is looking towards the ground

    
%% parameters for the adjustable paddle

pa.paddleHalfWidth = 0.075;% m
pa.paddleHalfHeight = 0.075;% m
pa.paddleHalfDepth = 0.075;% m
pa.aboveground = 0.15; %0.15;
% pa.paddleHeightFactor = 1;% 0.0057 m 
pa.paddleAngle = 0; % deg - start at the rightward position
pa.shiftPaddle = 0.25;
pa.rotationSpeed = pa.shiftPaddle; % this will update according to the observer's response
pa.acceleratePaddle = 0.008;
pa.paddleOrbitShift = 0.1185;% m
 

%% parameters for the target - will establish the speed distributions below

pa.targetMotionDuration = .5; % 1 s
pa.targetContrast = [1 0.15 0.075]; % fully-visible target, 15% and 7.5% contrast
pa.targetRadius = .25; % deg  0.25;
pa.targetSize = 0.025;% m
pa.fixationSize = pa.targetSize/2;
pa.nball = 50; %number of randomly placed objects
pa.ndots = pa.nball; % number of dots in dot task

%% experimental structure/design

pa.trialNumber = 0; % gotta start somewhere

% Set up a full factorial design 
pa.nRepeats = 20; % each target contrast condition gets pa.nRepeats trials - 75*3 = 225 per block 
% pa.speed = [0.5, 0.4, 0.3, 0.2,0.1]; %speeds m/s
% pa.speed = [0.5, 0.15, 0.075,  0.0375]; %speeds m/s
pa.speed = [0.5,.15,0.075, 0.0375]; %speeds m/s


pa.direction = deg2rad([90, 210, 270, 330, 352]) ; %(0 is radially out horizontally, 90 is forward/up, 270 is backwards/down x-z/x-y plane)

% pa.direction = deg2rad([90 230, 260, 290]) ; %(0 is to the right, 90 is forward, 270 is backwards)
factorial = fullfact([length(pa.speed), length(pa.direction)]); 
factorial = repmat(factorial,pa.nRepeats,1); % repeat the full factorial design nRepeats times
pa.fullFactorial = NaN(size(factorial));


% *Gaussian* vx and vz sampling each from a distribution with mean=0 and
% stdev=0.02m - the 0.061 m/s and -0.061 m/s would then be 3 stdevs from the mean,
% allowing us to cut off only a small percentage of speeds that come up
% for rtr=1:pa.nTrials 
%     % sample the vx and vz independently
%         inRange = 0;
%         while inRange==0
%             pa.xSpeedValues(rtr) = 0.005; %.*randn(1,1); % mean=0, stdev=0.02 
%             pa.zSpeedValues(rtr) = 0.005; %.*randn(1,1); % mean=0, stdev=0.02
%             if vecnorm([pa.xSpeedValues(rtr) pa.zSpeedValues(rtr)])<.0075
%                 inRange = 1;
%             end
%         end
% end

% pa.fullFactorial(:,1) = pa.targetContrast(pa.fullFactorial(:,1)); % instead of having this be indices of the contrast values it instead will have the actual contrast values - much nicer for data analysis
pa.fullFactorial(:,1) = (pa.speed(factorial(:,1)).*cos(pa.direction(factorial(:,2))));
pa.fullFactorial(:,2) = -(pa.speed(factorial(:,1)).*sin(pa.direction(factorial(:,2))));
pa.fullFactorial(:,3) = pa.speed(factorial(:,1)); % keep track in terms of speed and direction, not just vector values
pa.fullFactorial(:,4) = pa.direction(factorial(:,2));
pa.nTrials = size(pa.fullFactorial,1);
pa.fullFactorial = pa.fullFactorial(randperm(pa.nTrials),:); % Randomly permute order of trial presentation
%pa.fullFactorial(end+1,:) = pa.fullFactorial(1,:); % repeat the first trial because it is effectively a "junk" trial

pa.LR = randi([0,1],1,pa.nTrials)*2-1;  %-1*ones(1, pa.nTrial)
pa.LRresponse = NaN(1,pa.nTrials);

pa.response = []; % Eventual response matrix
pa.quitFlag = 0; % don't give up
end




