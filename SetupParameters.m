function [ds,pa] = SetupParameters(ds)

% Establishes the parameters for the Oculus 3D Pong experiment using the
% 2nd generation Oculus (DK2)

% At the start of the session, should only update lines 11 and 12


%% Basic experimental specs 

pa.subjectName = 'test';
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
pa.gazeangle = deg2rad(15); %angle camera is looking at the ground
pa.translation = .5; %m/s forward

pa.fixationDiskRadius = 2.5; % deg
pa.fixationHalfSquareSize = 1;  % deg - just the 'box' in which the fixation lines will be drawn
pa.fixationHalfSquare = 0.0157; % m 
pa.fixationLineLength = 0.0079; % m 
    
pa.apertureRadius = 7.5;  % deg
   
pa.floorHeight = -1; % m
pa.ceilingHeight = 1.5; % m 
    
%% parameters for the adjustable paddle

pa.paddleHalfWidth = 0.05;% m
pa.paddleHalfHeight = 0.25;% m
pa.paddleHalfDepth = 0.05;% m
pa.paddleHeightFactor = 0.0057;% m 
pa.paddleAngle = 0; % deg - start at the rightward position
pa.shiftPaddle = 0.25;
pa.rotationSpeed = pa.shiftPaddle; % this will update according to the observer's response
pa.acceleratePaddle = 0.008;
pa.paddleOrbitShift = 0.1185;% m
 

%% parameters for the target - will establish the speed distributions below

pa.targetMotionDuration = 2; % 1 s
pa.targetContrast = [1 0.15 0.075]; % fully-visible target, 15% and 7.5% contrast
pa.targetRadius = .25; % deg  0.25;
pa.targetSize = 0.025;% m
pa.fixationSize = pa.targetSize/2;
pa.nball = 50; %number of randomly placed objects
  

%% experimental structure/design

pa.trialNumber = 0; % gotta start somewhere

% Set up a full factorial design 
pa.nRepeats = 75; % each target contrast condition gets pa.nRepeats trials - 75*3 = 225 per block 
pa.fullFactorial = fullfact(length(pa.targetContrast)); 
pa.fullFactorial = repmat(pa.fullFactorial,pa.nRepeats,1); % repeat the full factorial design nRepeats times
pa.nTrials = size(pa.fullFactorial,1);
pa.LR = randi([0,1],1,pa.nTrials)*2-1;
pa.LRresponse = NaN(1,pa.nTrials);

% *Gaussian* vx and vz sampling each from a distribution with mean=0 and
% stdev=0.02m - the 0.061 m/s and -0.061 m/s would then be 3 stdevs from the mean,
% allowing us to cut off only a small percentage of speeds that come up
for rtr=1:pa.nTrials 
    % sample the vx and vz independently
        inRange = 0;
        while inRange==0
            pa.xSpeedValues(rtr) = 0.05; %.*randn(1,1); % mean=0, stdev=0.02 
            pa.zSpeedValues(rtr) = 0.05; %.*randn(1,1); % mean=0, stdev=0.02
            if vecnorm([pa.xSpeedValues(rtr) pa.zSpeedValues(rtr)])>.01
                inRange = 1;
            end
        end
end

pa.fullFactorial(:,1) = pa.targetContrast(pa.fullFactorial(:,1)); % instead of having this be indices of the contrast values it instead will have the actual contrast values - much nicer for data analysis
pa.fullFactorial(:,2) = pa.xSpeedValues';
pa.fullFactorial(:,3) = pa.zSpeedValues';
pa.fullFactorial = pa.fullFactorial(randperm(pa.nTrials),:); % Randomly permute order of trial presentation
pa.fullFactorial(end+1,:) = pa.fullFactorial(1,:); % repeat the first trial because it is effectively a "junk" trial
pa.nTrials = size(pa.fullFactorial,1);

pa.response = []; % Eventual response matrix
pa.quitFlag = 0; % don't give up
end



