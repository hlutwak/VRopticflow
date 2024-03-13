function [ds,pa] = SetupParameters(ds)

% Establishes the parameters for the Oculus 3D Pong experiment using the
% 2nd generation Oculus (DK2)

% At the start of the session, should only update lines 11 and 12


%% Basic experimental specs 

pa.subjectName = 'JO-monocular';
pa.block = 1;
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

    
%% parameters for the adjustable paddle

pa.paddleHalfWidth = 0.075;% m
pa.paddleHalfHeight = 0.075;% m
pa.paddleHalfDepth = 0.075;% m
pa.smallFloorWidth = 1;
pa.aboveground = .15; %0.15;
pa.objectdist = 2;
pa.fixationdist = 3;
% pa.paddleHeightactor = 1;% 0.0057 m 
pa.paddleAngle = 0; % deg - start at the rightward position
pa.shiftPaddle = 0.25;
pa.rotationSpeed = pa.shiftPaddle; % this will update according to the observer's response
pa.acceleratePaddle = 0.008;
pa.paddleOrbitShift = 0.1185;% m
 
pa.gazeangle = atan(-pa.floorHeight/pa.fixationdist); %angle camera is looking towards the ground

%% parameters for the target - will establish the speed distributions below

pa.targetMotionDuration = .5; % 1 s
pa.targetContrast = [1 0.15 0.075]; % fully-visible target, 15% and 7.5% contrast
pa.targetRadius = .25; % deg  0.25;
pa.targetSize = 0.025;% m .025
pa.fixationSize = pa.targetSize;
pa.nball = 50; %number of randomly placed objects
pa.ndots = pa.nball; % number of dots in dot task

%% experimental structure/design

pa.trialNumber = 0; % gotta start somewhere

% Set up a full factorial design 
pa.nRepeats = 20; % each target contrast condition gets pa.nRepeats trials - 75*3 = 225 per block 

if pa.block == 1
    pa.nPractice = 11; %first trial is "junk" trial
else
    pa.nPractice = 6;
end
% pa.speed = [0.5, 0.4, 0.3, 0.2,0.1]; %speeds m/s
% pa.speed = [0.5, 0.35, 0.2, .1, .05]; %speeds m/s
pa.speed = [.5,.3,.1,.05];
pa.practice_speed = 1;
% pa.speed = [.5 .5*(2^-1) .5*(2^-3) .5*(2^-4)]; %speeds m/s

pa.direction = deg2rad([75 90 105 255 270 285]); %(0 is radially out horizontally, 90 is forward/up, 270 is backwards/down x-z/x-y plane)

% pa.direction = deg2rad([90 230, 260, 290]) ; %(0 is to the right, 90 is forward, 270 is backwards)
factorial = fullfact([length(pa.speed), length(pa.direction)]); 
practice_directions = randsample(factorial(:,2), pa.nPractice);
factorial = repmat(factorial,pa.nRepeats,1); % repeat the full factorial design nRepeats times
pa.fullFactorial = NaN(size(factorial));

% create factorial matrix, xspeed, yspeed, speed and direction
pa.fullFactorial(:,1) = (pa.speed(factorial(:,1)).*cos(pa.direction(factorial(:,2))));
pa.fullFactorial(:,2) = -(pa.speed(factorial(:,1)).*sin(pa.direction(factorial(:,2))));
pa.fullFactorial(:,3) = pa.speed(factorial(:,1)); % keep track in terms of speed and direction, not just vector values
pa.fullFactorial(:,4) = pa.direction(factorial(:,2));
nTrials = length(pa.fullFactorial);
pa.fullFactorial = pa.fullFactorial(randperm(nTrials),:); % Randomly permute order of trial presentation

% add 10 practice trials at high speed and randomly chosen of the set
% directions
practice = NaN(pa.nPractice, 4);
practice(:,1) = (pa.practice_speed*ones(pa.nPractice,1).*cos(pa.direction(practice_directions)'));
practice(:,2) = (pa.practice_speed*ones(pa.nPractice,1).*sin(pa.direction(practice_directions)'));
practice(:,3) = pa.practice_speed*ones(pa.nPractice,1);
practice(:,4) = pa.direction(practice_directions);

pa.fullFactorial = [practice; pa.fullFactorial];

%pa.fullFactorial(end+1,:) = pa.fullFactorial(1,:); % repeat the first trial because it is effectively a "junk" trial
pa.nTrials = size(pa.fullFactorial,1);

pa.LR = randi([0,1],1,pa.nTrials)*2-1;  %-1*ones(1, pa.nTrial)
pa.LRresponse = NaN(1,pa.nTrials);

pa.response = []; % Eventual response matrix
pa.quitFlag = 0; % don't give up
end




