       %                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        14-Dec-2015  jf Written. Derived  from OculusSDK2PongDemo_Fixed.m
%  6-Jan-2016  jf Edited to improve lag  ged condition performance   
% 14-Jan-2016 jf  Switched over to Windows pla tform and optimized the code
 % for timing and stimulus presentation - including now measured gamma
 % correction   
 % 19-Aug-2016 jf Added a few modifications: randomized padd le start angle
% on each trial, feedback options, random/variable lag for lagged condition
% Dec-Jan-2019 - JF updated code to  work with CV1; minor  changes to call
 % projection matrices, added in  CV1-specific FOV and ot her parameters
% Jan 2023 - HL changed to detecting object motion experi ment
%% Important note about coding of angles in Oculus space: 

% In the Oculus, angles are coded ccw - so, straight righ 1t = 0 deg,
% directly in front of fix ation = 90 deg, straight left = 180 deg, directly
% behind fixation = 270 deg        


%% Basic per subject inputs:

% Use SetupDisplay.m to establish experimental condition as active, fixed, or lagged (line 51, ds.experimentType)
% Use SetupParameters.m to establish participant number for the data file names (line 11, pa.subjectName and line 12, pa.feedbackFlag)

%% IMPORTANT!  The Oculus must be plugged in and turned on *before* starting MATLAB
% The individual must be wearing the device prior to starting the
% experimental program to achieve proper frame rate
% Also the Oculus VR runtime version '0.5.0.1' *must* be installed for PTB
% to properly interact with/recognize the Oculus DK2, otherwise, use the
% latest version of the runtime for the CV1 (note: cannot have 2 versions
% on one machine)
% addpath(genpath('C:\Users\hlutw\Documents\MATLAB\Psychtoolbox-3-master\Psychtoolbox\PsychOpenGL\MOGL\wrap'));

clear all;
close all;

global DEBUG_FLAG KEYBOARD_FLAG
DEBUG_FLAG = 1 ; %1
KEYBOARD_FLAG = 0; %1 % there is a call to 'keyboard' in SetupDisplay that was breaking the code when using the hmd so I created a flag to turn it on/off - not certain what it is for (JF)

if DEBUG_FLAG
    Screen('Preference', 'SkipSyncTests', 1); % For debugging
end

% Set up Psychtoolbox for OpenGL 3D rendering support and initialize the
% mo gl OpenGL for Matlab/Octave wrapper:
global GL; % GL data structure needed for all OpenGL programs
InitializeMatlabOpenGL(1);
PsychDefaultSetup(2); % the input of 2 means: execute the AssertOpenGL command, execute KbName('UnifyKeyNames') routine, AND unifies the color mode and switches from 0-255 to 0.0-1.0 - color part only impacts the current function or script, not ones that are called
 
addpath(genpath([pwd filesep() 'Tools'])); % contains 'isodd.m' and 'oneoverf.m' for texture rendering

% Initialize screen, experimental parameters, etc.
[ds,oc] = SetupDisplay(); %(oc); % set up the  display, based on the DK2
[ds,pa] = SetupParameters(ds); % set up the experimental parameters for this session
[ds,pa] = CreateTextures(ds, pa); % create the surround & paddle face textures as well as the ceiling, floor, and walls of the virtual room - just needs to be done once
kb = SetupKeyboard(); % get the keyboard info for the participant's responses
% ListenChar(2);
if ~DEBUG_FLAG
    HideCursor(ds.screenId);
    ListenChar(2); % Stop making keypresses show up in the matlab scripts and
    %command window - it's really annoying and can cause all sorts of problems! % CSB: debug
end

%% Run eye tracking
pause(.5) 
finishedCalibration = 0;
readyToBegin = 0;   

if ds.eyetracking
    device = SetupEyetracker();
    disp(device.phone_name) 
    disp(device.module_serial)
    
    if device.module_serial == string(missing) || ~exist('device')
        error('Neon not connected')
    end
    
    recording_id = device.recording_start();
    disp(['Started recording with id', string(recording_id)]);
else
    device = [];
end

while ~finishedCalibration && ~readyToBegin
       
    % Camera position when using head tracking + HMD: (according to SuperShapeDemo.m)
    globalPos = [0, 0, 0]; % x,y,z  % in meters - just put something in here for now, will likely be much larger later for viewing the tv/'real' world - the demos use large values too
    heading = 0; % yaw
    ds.globalHeadPose = PsychGetPositionYawMatrix(globalPos, heading); % initialize observer's start position to the default camera position specified above
    
    if isempty(ds.hmd) % Oculus not connected
        % This is already done in SetupDisplay
        % load DefaultHMDParameters.mat;
        % oc.defaultState = defaultState;
        % Some stuff needs to be done here to get a proper initialState
        % There is something messed up with using both initial and default
        % state
        % can we simplify using only one or the other?
        % oc.initialState = defaultState.initialState;    
        oc.initialState = oc.defaultState;
        oc.initialState.modelView{1} = oc.defaultState.modelViewDataLeft;
        oc.initialState.modelView{2} = oc.defaultState.modelViewDataRight;
        
    else % Oculus connected
        oc.initialState = PsychVRHMD('PrepareRender', ds.hmd, ds.globalHeadPose);  % get the state of the hmd now
    end
    
    for renderPass = 0:1 % loop over eyes
        ds.renderPass = renderPass;
        Screen('SelectStereoDrawBuffer',ds.w,ds.renderPass);
        Screen('BeginOpenGL',ds.w);
                            
        % Setup camera position and orientation for this eyes view:
        glMatrixMode(GL.PROJECTION)
        glLoadMatrixd(ds.projMatrix{renderPass + 1});
        
%         modelView = oc.initialState.modelView{ds.renderPass + 1}; % Use per-eye modelView matrices
        modelView = [1 0 0 0; 0 1 0 0; 0 0 1 -ds.viewingDistance; 0 0 0 1];
        glLoadMatrixd(modelView); 
        

        glClearColor(.5,.5,.5,1); % gray background
        
        
        glClear(); % clear the buffers - must be done for every frame
        glColor3f(1,1,1);
        
            % fixation target
            glPushMatrix;
            glTranslatef(-.5,0, pa.objectdist);
            glCallList(ds.incorrect);
            glPopMatrix;
            
            % testing eye calibration
%         w = .5;
%         dist = -2;
%         coordpos = [0,0,dist; -w,w,dist; w,w,dist; w,-w,dist; -w,-w,dist];
%         ncoordpts = size(coordpos,1);
%                 for b = 1:ncoordpts
%                     glPushMatrix;
%                     glTranslatef(coordpos(b,1),coordpos(b,2),coordpos(b,3)); 
%                     glCallList(ds.highcontrastTarget);
%                     glPopMatrix;
%                 end
            
            
              
                
        Screen('EndOpenGL', ds.w);
%         Screen('DrawText',ds.w,'Eye Calibration. Press SPACE to end.',(ds.textCoords(1)-200*ds.renderPass)-100,ds.textCoords(2),[1 1 1]);
         %ds.pixelsPerDegree*ds.hFOV/4
%           calibx = ds.xc + round([0,ds.windowRect(4)/6, -ds.windowRect(4)/6, -ds.windowRect(4)/6,ds.windowRect(4)/6])- [200 * ds.renderPass]; %
%           caliby = ds.yc + round([0,ds.windowRect(3)/6, ds.windowRect(3)/6, -ds.windowRect(3)/6,-ds.windowRect(3)/6]);%
%           calibrationGridDeg = round([0,ds.hFOV/4, ds.hFOV/4, -ds.hFOV/4,-ds.hFOV/4; 0,ds.vFOV/4, ds.vFOV/4, -ds.vFOV/4, -ds.vFOV/4]);

          calibx = ds.xc + round([0,100, -100, -100, 100])- [200 * ds.renderPass]; %
          caliby = ds.yc + round([0,100, 100, -100, -100]);%
          
        calibrationGridDeg = [calibx-ds.xc+[200 * ds.renderPass];caliby-ds.yc]/12; %ds.pixelsPerDegree
        
        Screen('DrawDots', ds.w,[calibx(1); caliby(1)], 10,[1,1,1]);
        Screen('DrawDots', ds.w,[calibx(2); caliby(2)], 10,[1,0,0]);
        Screen('DrawDots', ds.w,[calibx(3); caliby(3)], 10,[0,1,0]);
        Screen('DrawDots', ds.w,[calibx(4); caliby(4)], 10,[0,0,1]);
        Screen('DrawDots', ds.w,[calibx(5); caliby(5)], 10,[0,0,0]);




          

%         Screen('DrawText',ds.w,'o',ds.xc,ds.yc,[1 1 1]);

    end
    
    Screen('DrawingFinished', ds.w);
    ds.vbl = Screen('Flip', ds.w);
     
    [kb.keyIsDown, kb.secs, kb.keyCode] = KbCheck(-1); % query the keyboard
    
    if ds.eyetracking 
        if kb.keyIsDown && kb.keyCode(kb.reorientKey) %use up arrow key
            device.send_event("calibrationpoint") 
            oc.UTCCalibrationPoint = datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss:ms');
        end
    end
    
    if kb.keyIsDown && kb.keyCode(kb.spacebarKey)
        finishedCalibration=1; 
        if ds.eyetracking && pa.trialNumber == 2 
            device.send_event("calibrationdone")   
        end
        oc.UTCCalibrationEnd = datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss:ms');
    
    end 
end
%% Start the experiment - opening screen, getting the participant set
pause(.5) 
readyToBegin = 0; 
while ~readyToBegin && finishedCalibration% confirm everything's ready to go
    
    % Camera position when using head tracking + HMD: (according to SuperShapeDemo.m)
    globalPos = [0, 0, 0]; % x,y,z  % in meters - just put something in here for now, will likely be much larger later for viewing the tv/'real' world - the demos use large values too
    heading = 0; % yaw
    ds.globalHeadPose = PsychGetPositionYawMatrix(globalPos, heading); % initialize observer's start position to the default camera position specified above
    
    if isempty(ds.hmd) % Oculus not connected
        % This is already done in SetupDisplay
        % load DefaultHMDParameters.mat;
        % oc.defaultState = defaultState;
        % Some stuff needs to be done here to get a proper initialState
        % There is something messed up with using both initial and default
        % state
        % can we simplify using only one or the other?
        % oc.initialState = defaultState.initialState;    
        oc.initialState = oc.defaultState;
        oc.initialState.modelView{1} = oc.defaultState.modelViewDataLeft;
        oc.initialState.modelView{2} = oc.defaultState.modelViewDataRight;
        
    else % Oculus connected
        oc.initialState = PsychVRHMD('PrepareRender', ds.hmd, ds.globalHeadPose);  % get the state of the hmd now
    end
    
    for renderPass = 0:1 % loop over eyes
        ds.renderPass = renderPass;
        Screen('SelectStereoDrawBuffer',ds.w,ds.renderPass);
        Screen('BeginOpenGL',ds.w);
                            
        % Setup camera position and orientation for this eyes view:
        glMatrixMode(GL.PROJECTION)
        glLoadMatrixd(ds.projMatrix{renderPass + 1});
        
        modelView = oc.initialState.modelView{ds.renderPass + 1}; % Use per-eye modelView matrices
        glLoadMatrixd(modelView); 
          
        Screen('EndOpenGL', ds.w);
        Screen('DrawText',ds.w,'Ready to start the experiment? Press SPACE to confirm.',(ds.textCoords(1)-200*ds.renderPass)-200,ds.yc,[1 1 1]);

    end
    
    Screen('DrawingFinished', ds.w);
    ds.vbl = Screen('Flip', ds.w);
    
    [kb.keyIsDown, kb.secs, kb.keyCode] = KbCheck(-1); % query the keyboard
    if kb.keyIsDown && kb.keyCode(kb.spacebarKey)
        readyToBegin=1;
    end
end


%% Participant is ready, so let's go
pause(1);
ds.tElapsed = 0;
ds.fCount = 0;

[ds, pa, kb, oc] = SetupNewTrial(ds, pa, kb, oc);
ds.vbl = pa.trialOnset;
tStart = ds.vbl;
pa.experimentOnset = ds.vbl;
oc.UTCtrialStart=[];
oc.UTCtrialEnd =[];

breakTime = 0;  % participants are running in the task
kb.nextTrialKey = 0;
track_trial = 0;
track = [];
track_theta = [];
track_dtheta = [];

while (pa.trialNumber <= pa.nTrials) && ~kb.keyCode(kb.escapeKey) % wait until all of the trials have been completed or the escape key is pressed to quit out
    
    % Get HMD state
    if isempty(ds.hmd) % Oculus is not connected - will display a poor imitation of the Oculus rift on your main computer screen
        state = oc.defaultState; % just set to a default, non-updating viewpoint
    else   % Oculus is connected - uses PTB's code + openGL code to display in the HMD
        % Track and predict head position and or ientation, retrieve modelview
        % camera matrices for rendering of each  eye. Apply some global transformation
        % to returned camera matrices. In this case a translation + rotation, as defined
        % by the PsychGetPositionYawMatrix() h elper function:
        state = PsychVRHMD('PrepareRender', ds.hmd, ds.globalHeadPose);  % Mark the start of the rendering cycle for a new 3D rendered stereoframe. Return a struct 'state' which contains various useful bits of information for 3D stereoscopic rendering of a scene, based on head tracking data
        oc.HMD = [oc.HMD; state.modelView{1}];
    end

    if pa.trialNumber>track_trial %dont' update head position during a trial
%        
       if pa.trialNumber>1
           eye.eyeIndex = 0;
%            eye.modelView = oc.modelViewDataRight(1:4,:);
%            eye.modelView = [1 0 0 0; 0 1 0 0; 0 0 1 -ds.viewingDistance; 0 0 0 1];
           eye.modelView = originaleye.modelView;

       else
           if isempty(ds.hmd)
               eye.modelView = oc.defaultState.modelViewDataRight;
                R = [1 0 0; 0 cos(pa.gazeangle) -sin(pa.gazeangle); 0 sin(pa.gazeangle) cos(pa.gazeangle)];
                   eye.modelView = [1 0 0 0; 0 1 0 0; 0 0 1 -ds.viewingDistance; 0 0 0 1];
                   eye.modelView(1:3,1:3) = eye.modelView(1:3,1:3)*R;
                   originaleye = eye;
                   theta = pa.gazeangle;
                   track_theta = [track_theta,theta];
           else %if hmd is connected
               eye = PsychVRHMD('GetEyePose', ds.hmd, ds.renderPass, ds.globalHeadPose);
               if ~ds.trackingFlag
                   R = [1 0 0; 0 cos(pa.gazeangle) -sin(pa.gazeangle); 0 sin(pa.gazeangle) cos(pa.gazeangle)];
                   eye.modelView = [1 0 0 0; 0 1 0 0; 0 0 1 -ds.viewingDistance; 0 0 0 1];
                   eye.modelView(1:3,1:3) = eye.modelView(1:3,1:3)*R;
                   originaleye = eye;
                   theta = pa.gazeangle;
                   track_theta = [track_theta,theta];
               else
               end
           end
        end
        track_trial = track_trial+1;
        
    end  
    
    % Render the scene separately for each eye:
    for renderPass = 0:1 %0 left, 1 right eye
        ds.renderPass = renderPass;
        
        if isempty(ds.hmd) % hmd not connected
            % Get head position from keyboard input
            
            eye.eyeIndex = ds.renderPass; % We are switching eye index, but not the eye.modelview here

%             % CSB: this inits camera matrix if it is empty
%             % BR: This should happen outside of the loop init'd on line 131
%             if ~isfield(eye,'modelView') % checks if camera matrix exists. it shouldn't on first trial and it will be init'd.
%                 if ds.renderPass==0 % drawing left eye
%                     eye.modelView = oc.defaultState.modelViewDataLeft;
%                 elseif ds.renderPass==1 % drawing right eye
%                     eye.modelView =  oc.defaultState.modelViewDataRight;
% %                     eye.modelView(1,4) =  eye.modelView(1,4)+100; % CSB: debug
%                 end
%             end
            
            [pa, kb, eye] = GetKeyboardHeadmotion(pa, ds, kb, eye);  % query the keyboard to allow the observer to rotate the paddle and eventually lock in his/her response to initiate a new trial
            oc.modelViewDataRight = [oc.modelViewDataRight; eye.modelView];
            
        else % hmd connected
            % Query which eye to render in this ds.renderPass, and query its
            % eyePose vector for the predicted eye position to use for the virtual
            % camera rendering that eyes view. The returned pose vector actually
            % describes tracked head pose, ie. HMD position and orientation in space.
            
            

            if ds.trackingFlag
                eye = PsychVRHMD('GetEyePose', ds.hmd, ds.renderPass, ds.globalHeadPose);
            else
                [pa, kb, eye] = GetKeyboardHeadmotion(pa,ds,kb,eye);

            end
            
            % this is for saving purposes to recreate participants' head motion

            if ds.renderPass % drawing right eye
                oc.modelViewDataRight = [oc.modelViewDataRight; eye.modelView];
            else % drawing left eye
                oc.modelViewDataLeft = [oc.modelViewDataLeft; eye.modelView];
            end
            
            % the 'active' condition just takes whatever the current state of the oculus is
                if ~ds.trackingFlag % 'fixed' condition loop - don't update the scene with tracked head motion, just use the default state
                    if ds.renderPass
%                         eye.modelView = oc.modelViewDataRight(1:4,:); % comes back from the initial call...will not update the scene based on head tracking
                    else 
%                         eye.modelView = oc.modelViewDataLeft(1:4,:);
                    end
                    state.tracked = 2;
                    eye.eyeIndex = ds.renderPass;
                elseif ds.trackingLag>0 && ds.lagNow==1 % 'lagged' condition loop - go back ds.trackingLag time steps if available
                    if (length(oc.modelViewDataRight)/4)>((ds.trackingLag+1)*4) && (length(oc.modelViewDataLeft)/4)>((ds.trackingLag+1)*4) % make sure enough time steps have passed to go back in time
                        stepbackind = length(oc.modelViewDataLeft) - ((ds.trackingLag+1)*4) + 1; % if there are enough time steps, figure out which to go back to
                        if ds.renderPass==0 % drawing left eye
                            eye.modelView = oc.modelViewDataLeft(stepbackind:stepbackind+3,:);
                        elseif ds.renderPass==1 % drawing right eye
%                             eye.modelView = oc.modelViewDataRight(stepbackind:stepbackind+3,:);
                        end
                    end
                end
        end
        
        Screen('SelectStereoDrawbuffer', ds.w, eye.eyeIndex); % Select 'eyeIndex' to render (left- or right-eye):
        modelView = eye.modelView; % Extract modelView matrix for this eye:

        Screen('BeginOpenGL', ds.w); % Manually reenable 3D mode in preparation of eye draw cycle
                    
%         % Setup camera position and orientation for this eyes view:
        glMatrixMode(GL.PROJECTION)
        glLoadMatrixd(ds.projMatrix{renderPass + 1});

        glMatrixMode(GL.MODELVIEW);
        glLoadMatrixd(modelView);  
        

        if ds.dotfield
            glClearColor(0,0,0,1); % gray background
        else
            glClearColor(.5,.5,.5,1); % gray background
        end
        
        glClear(); % clear the buffers - must be done for every frame
        glColor3f(1,1,1);
        
        glPushMatrix;
        
        if ds.binocular || (~ds.binocular && ds.renderPass)
            % fixation target
            glPushMatrix;
            glTranslatef(0,pa.floorHeight+pa.fixationSize,-pa.fixationdist);
            glCallList(ds.highcontrastTarget);
            glPopMatrix;
            
%             glPushMatrix;
%             glTranslatef(0,pa.floorHeight+pa.fixationSize,-pa.floorWidth/2+1);
%             glCallList(ds.dot);
%             glPopMatrix;


% testing eye calibration
%         w = .5;
%         dist = -2;
%         coordpos = [0,0,dist; -w,w,dist; w,w,dist; w,-w,dist; -w,-w,dist];
%         ncoordpts = size(coordpos,1);
%                 for b = 1:ncoordpts
%                     glPushMatrix;
%                     glTranslatef(coordpos(b,1),coordpos(b,2),coordpos(b,3)); 
%                     glCallList(ds.highcontrastTarget);
%                     glPopMatrix;
%                 end
            
            
            
        end
        
       
          
            
         %% Experiment Logic  
       
        if ds.vbl <  pa.trialOnset + pa.targetMotionDuration % if current time < present until time: draw target, 1 s target motion
            
            oc.trial_startflag = [1 oc.trial_startflag];
            if oc.trial_startflag(1)-oc.trial_startflag(2) == 1
                oc.UTCtrialStart = [oc.UTCtrialStart datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss:ms')];
                if ds.eyetracking && pa.trialNumber == 2
                    device.send_event("trialstart")
                end
            end
            
            % Target position 
            % make the onset delayed
            delay = 0;
            if ds.vbl <  pa.trialOnset + delay
                t = 0;
            else
                t = ds.vbl-pa.trialOnset-delay;
            end
            xPosition =pa.xSpeed.*t;  %pa.xSpeed.*t;
            zPosition = pa.zSpeed.*t;  %pa.zSpeed.*t
            
            if ds.binocular || (~ds.binocular && ds.renderPass)
                glPushMatrix;
                %x-z plane
                glTranslatef((xPosition+.5)*(pa.LR(pa.trialNumber)),pa.floorHeight+pa.paddleHalfHeight+pa.aboveground,zPosition-pa.objectdist); % shift the target to its position along its trajectory for this frame

                if ds.dotfield
                    glCallList(ds.fixation);
                else
                    glCallList(ds.paddleList);
                end

                    glPopMatrix;

                % stationary object
                glPushMatrix;
                glTranslatef(.5*(-pa.LR(pa.trialNumber)),pa.floorHeight+pa.paddleHalfHeight+pa.aboveground,-pa.objectdist); 
                if ds.dotfield
                    glCallList(ds.fixation);
                else
                    glCallList(ds.paddleList);
                end
                glPopMatrix;


                % place random stationary objects
                if ~ds.control
                    for b = 1:pa.nball
                        glPushMatrix;
                        glTranslatef(pa.positions(1,b),pa.floorHeight+pa.paddleHalfHeight+pa.positions(2,b),pa.positions(3,b)); 
                        glRotatef(pa.rotations(b,1), pa.rotations(b,2), pa.rotations(b,3), pa.rotations(b,4));
                    if ds.dotfield
                        glCallList(ds.fixation);
                    else
                        glCallList(ds.paddleList);
                    end
                    glPopMatrix;
                    end
                end
                
                if ds.dotfield
                    
                    if ~ds.control
                        for b = 1:pa.nball*5 %fill in floor with dots
                            glPushMatrix;
                            glTranslatef(pa.dotpositions(1,b),pa.dotpositions(2,b),pa.dotpositions(3,b));
                            glCallList(ds.fixation);
                            glPopMatrix;
                        end
                    
                    else
                        
                        for b = 1:pa.nball %fill in floor with dots
                            glPushMatrix;
                            glTranslatef(pa.surroundpositions(1,b),pa.surroundpositions(2,b),pa.surroundpositions(3,b));
                            glCallList(ds.fixation);
                            glPopMatrix;
                        end
                    end
                    
                end
            end
            
              

            if ds.eyesimulation % simulate eye rotation along with translation
%                 eye.modelView(3,4) =  pa.translation.*(ds.vbl-pa.trialOnset); %eye.modelView(3,4) + %eye.modelView(2,4) - %(ds.vbl-pa.trialOnset)
                track = [track (ds.vbl-pa.trialOnset)];
                    theta = atan(-pa.floorHeight./(pa.fixationdist-pa.translation*(ds.vbl-pa.trialOnset))); % update theta for observer fixating at a point at the ground in front of them, fixation m away
                    dtheta = theta - track_theta(end);
                    track_dtheta = [track_dtheta, dtheta];
                    R = [1 0 0; 0 cos(dtheta) -sin(dtheta); 0 sin(dtheta) cos(dtheta)]; % viewing angle is changed upon each update, so only rotate the difference in what the angle should be
                    eye.modelView(1:3,1:3) = eye.modelView(1:3,1:3)*R;
                    track_theta = [track_theta, theta];
                    eye.modelView(2,4) = -pa.translation.*(ds.vbl-pa.trialOnset)*sin(theta); %adjust so translation is strictly foward, absolute world coordinate
                    eye.modelView(3,4) =  pa.translation.*(ds.vbl-pa.trialOnset); %eye.modelView(3,4) + %eye.modelView(2,4) - %(ds.vbl-pa.trialOnset)

            else % just do translation
                eye.modelView(3,4) =  pa.translation.*(ds.vbl-pa.trialOnset); %eye.modelView(3,4) + %eye.modelView(2,4) - %(ds.vbl-pa.trialOnset)
                track = [track (ds.vbl-pa.trialOnset)];
                eye.modelView(2,4) = -pa.translation.*(ds.vbl-pa.trialOnset)*sin(pa.gazeangle); %adjust so translation is strictly foward
            
            end
            

            pa.responseOnset = ds.vbl; % start the timer on the response time
            

        elseif ds.vbl >  pa.trialOnset + pa.targetMotionDuration && ~oc.trialendflag
            oc.UTCtrialEnd = [oc.UTCtrialEnd datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss:ms')];

            oc.trialendflag = 1;
            
        elseif ~kb.responseGiven && oc.trialendflag % show paddle and allow observers to adjust its position - no time constraints - they press the space bar to lock in their response and start a new trial
            %2, % debugging flag
            

            pa.modelView = eye.modelView; % CSB: init camera position based on default or previous camera position. June 7th 2018
            [pa, kb] = GetResponse(pa, ds, kb);  % query the keyboard to allow the observer to rotate the paddle and eventually lock in his/her response to initiate a new trial
            eye.modelView = pa.modelView; % CSB: update camera position based on keys pressed. June 7th 2018
%             
%             theta = pi*pa.xS peed.*(ds.vbl-pa.trialOnset); 
%             track = [track theta];
%             rotationMatrixYawRight = [cos(theta) 0 sin(theta); 0 1 0; -sin(theta) 0 cos(theta);];
%             pa.modelView = eye.modelView;
%             pa.modelView(1,4) = pa.modelView(1,4) + pa.xSpeed.*(ds.vbl-pa.trialOnset);
%             track = [track pa.modelView(1,4)];
%             eye.modelView = pa.modelView;
            
            pa.feedbackOnset = ds.vbl;
            

        elseif (kb.responseGiven && pa.feedbackFlag==1 && pa.feedbackGiven==0) || (kb.responseGiven && pa.feedbackFlag==2 && pa.feedbackGiven==0) % sound feedback only
            %4, % debugging flag
            
            reportedAngle = angle(exp(1i*deg2rad(pa.paddleAngle(pa.thisTrial))));  % converts the paddle angle from 0-360 to -pi - + pi
            
            % display appropriate object
            if ~pa.feedbackGiven

                pa.feedbackGiven = 1;
                pa.waitTime = ds.vbl;
                
                % start any lag now while waiting for subject to initiate new trial
                if ds.trackingFlag==1 && ds.trackingLagStart>0 % screen will lag in updating with head movements during target motion
                    ds.trackingLag = randi(39)-1; % lag in frames: random lag between 0 and 500 ms
                    ds.lagNow = 1;
                else
                    ds.trackingLag = 0;
                    ds.lagNow = 0;
                end
                
            end
            
        elseif (kb.responseGiven && pa.feedbackFlag==1 && pa.feedbackGiven==1 && ds.vbl <= pa.waitTime+0.5) || (kb.responseGiven && pa.feedbackFlag==2 && pa.feedbackGiven==1 && ds.vbl <= pa.waitTime+.5)


            if ds.binocular || (~ds.binocular && ds.renderPass)
                if pa.LRresponse(pa.trialNumber) == pa.LR(pa.trialNumber)  %does the response match the target side
                        glPushMatrix;
                        glTranslatef(0,pa.floorHeight+pa.fixationSize,-(pa.fixationdist)); % display green sphere
                        glCallList(ds.correct);
                        glPopMatrix;
                else
                        glPushMatrix;
                        glTranslatef(0,pa.floorHeight+pa.fixationSize,-(pa.fixationdist)); % display red sphere
                        glCallList(ds.incorrect);
                        glPopMatrix;
                end
            end
                
        elseif (kb.responseGiven && pa.feedbackFlag==0) || (kb.responseGiven && pa.feedbackFlag==1 && pa.feedbackGiven==1) || (kb.responseGiven && pa.feedbackFlag==2 && pa.feedbackGiven==1) % done, set up for the next trial (i.e., determine the new random trajectory)
            %                 5, % debugging flag
            
            [ds, pa, kb, oc] = SetupNewTrial(ds, pa, kb, oc);
            track_trial = pa.trialNumber-1;
            
        end
        %% End of experiment logic
        
        % Update view with keyboard responses
%         glMatrixMode(GL.MODELVIEW); 
%         glLoadMatrixd(modelView);
%         
        glPopMatrix; % reset back to the origin
        moglDrawDots3D(ds.w, [inf inf inf], 10, [1 1 1 1], [], 2); % 'hack' to fix transparency issue - same one we've used in all 3D pong experiments - definitely works!
        
        % Distance debugging code
%         glCallList(ds.debugcylinder); % call the pre-compiled list that binds the textures to the paddle faces
        
%         glBindTexture(GL.TEXTURE_2D,ds.roomwall_texid); % was suggested to bind textures before/outside of call lists rather than in - doesn't buy us anything from what I can tell though
%         glCallList(ds.wallTexture);
%         glBindTexture(GL.TEXTURE_2D,ds.ceiling_texid);
%         gl CallList(ds.ceilingTexture);
         if ds.binocular || (~ds.binocular && ds.renderPass)
             
             if ~ds.dotfield
                 glBindTexture(GL.TEXTURE_2D,ds.floor_texid);
                 glCallList(ds.floorTexture);
             end
            
         end
        
%         glBindTexture(GL.TEXTURE_2D,ds.wall_texid); 
%         glCallList(ds.surroundTexture); % 1/f noise texture surround -  comes from CreateTexturesforSDK2.m
%           
            glBindTexture(GL.TEXTURE_2D, 0);
        
        
        % Manually disable 3D mode before switching to other eye or to flip:
        Screen('EndOpenGL', ds.w);
        
        % Compute simulation time for this draw cycle:
        ds.tElapsed = (ds.vbl - tStart) * 1;
        
        % Repeat for renderPass of other eye
    end 
    
    % Head position tracked in the HMD?
    if ~isempty(ds.hmd)
        if ~bitand(state.tracked, 2) && ds.trackingFlag==1
            % Nope, user out of cameras view frustum. Tell it like it is:
            DrawFormattedText(ds.w, 'Vision based tracking lost\nGet back into the cameras field of view!', 'center', 'center', [1 0 0]);
        end
    end
    
    % Stimulus ready. Show it on the HMD. We don't clear the color buffer here,
    % as this is done in the next iteration via glClear() call anyway:
    Screen('DrawingFinished', ds.w);
    
    Screen('Flip', ds.w,[],[],1);%, [], [], 2);%, ds.vbl + (1-0.5) * ds.ifi);
    ds.vbl = GetSecs();
    ds.fCount = ds.fCount + 1;
end
% Calculate average framerate:
ds.fps = ds.fCount / (ds.vbl - tStart); % uncomment to print out at end of run

pa.dataFile = fullfile(  pa.baseDir, 'Data', [pa.subjectName '-' num2str(pa.block) '-' pa.date '.mat']);
save(pa.dataFile, 'pa', 'ds', 'kb','oc');

if ds.eyetracking
    device.recording_stop_and_save()
end

% Calculate average framerate:

% Done (or quit out). Save data (in pa.response) and other relevant parameters/variables, close screen and exit
PsychVRHMD('SetAutoClose', ds.hmd, 2);
PsychPortAudio('Close', pa.handleHit);
PsychPortAudio('Close', pa.handleMiss);
Priority(0);
ShowCursor(ds.screenId);
ListenChar(1); % Let's start listening to the keyboard again
sca;
