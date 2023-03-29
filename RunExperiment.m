   % function RunExperiment
  
% 14-Dec-2015  jf Written. Derived from OculusSDK2PongDemo_Fixed.m   
% 6-Jan-2016  jf Edited to improve lagg      ed condition performance
% 14-Jan-2016 jf  Switched over to Windows platform and optimized the code
% for timing and stimulus presentation - including now measured gamma
% correction
% 19-Aug-2016 jf Added a few modifications: randomized paddle start angle
% on each trial, feedback options, random/variable lag for lagged condition
% Dec-Jan-2019 - JF updated code to  work with CV1; minor changes to call of
% projection matrices, added in  CV1-specific FOV and other parameters

%% Important note about coding of angles in Oculus space:

% In the Oculus, angles are coded ccw - so, straight right = 0 deg,
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
clear all;
close all;

global DEBUG_FLAG KEYBOARD_FLAG
DEBUG_FLAG = 1; %1
KEYBOARD_FLAG = 0; %1 % there is a call to 'keyboard' in SetupDisplay that was breaking the code when using the hmd so I created a flag to turn it on/off - not certain what it is for (JF)
 
if DEBUG_FLAG
    Screen('Preference', 'SkipSyncTests', 1); % For debugging
end

% Setup Psychtoolbox for OpenGL 3D rendering support and initialize the
% mogl OpenGL for Matlab/Octave wrapper:
global GL; % GL data structure needed for all OpenGL programs
InitializeMatlabOpenGL(1);
PsychDefaultSetup(2); % the input of 2 means: execute the AssertOpenGL command, execute KbName('UnifyKeyNames') routine, AND unifies the color mode and switches from 0-255 to 0.0-1.0 - color part only impacts the current function or script, not ones that are called

addpath(genpath([pwd filesep() 'Tools'])); % contains 'isodd.m' and 'oneoverf.m' for texture rendering

% Initialize screen, experimental parameters, etc.
[ds,oc] = SetupDisplay(); %(oc); % set up the  display, based on the DK2 
[ds,pa] = SetupParameters(ds); % set up the experimental parameters for this session
[ds,pa] = CreateTextures(ds, pa); % create the surround & paddle face textures as well as the ceiling, floor, and walls of the virtual room - just needs to be done once
kb = SetupKeyboard(); % get the keyboard info for the participant's responses
ListenChar(2);
if ~DEBUG_FLAG
    HideCursor(ds.screenId); 
    ListenChar(2); % Stop making keypresses show up in the matlab scripts and
    %command window - it's really annoying and can cause all sorts of problems! % CSB: debug
end

%% Start the experiment - opening screen, getting the participant set
 
readyToBegin = 0; 
while ~readyToBegin % confirm everything's ready to go
    
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
        Screen('DrawText',ds.w,'Ready to start the experiment? Press SPACE to confirm.',(ds.textCoords(1)-200*ds.renderPass)-100,ds.textCoords(2),[1 1 1]);

    end
    
    Screen('DrawingFinished', ds.w);
    ds.vbl = Screen('Flip', ds.w);
    
    [kb.keyIsDown, kb.secs, kb.keyCode] = KbCheck(-1); % query the keyboard
    if kb.keyIsDown && kb.keyCode(kb.spacebarKey)
        readyToBegin=1;
    end
end


%% Participant is ready, so let's go
ds.tElapsed = 0;
ds.fCount = 0;

[ds, pa, kb] = SetupNewTrial(ds, pa, kb);
ds.vbl = pa.trialOnset;
tStart = ds.vbl;
pa.experimentOnset = ds.vbl;

pa.block = 0;
breakTime = 0;  % participants are running in the task
kb.nextTrialKey = 0;
track_trial = 0;
track = [];

while (pa.trialNumber <= pa.nTrials) && ~kb.keyCode(kb.escapeKey) % wait until all of the trials have been completed or the escape key is pressed to quit out
    
    % Get HMD state
    if isempty(ds.hmd) % Oculus is not connected - will display a poor imitation of the Oculus rift on your main computer screen
        state = oc.defaultState; % just set to a default, non-updating viewpoint
    else   % Oculus is connected - uses PTB's code + openGL code to display in the HMD
        % Track and predict head position and orientation, retrieve modelview
        % camera matrices for rendering of each eye. Apply some global transformation
        % to returned camera matrices. In this case a translation + rotation, as defined
        % by the PsychGetPositionYawMatrix() helper function:
        state = PsychVRHMD('PrepareRender', ds.hmd, ds.globalHeadPose);  % Mark the start of the rendering cycle for a new 3D rendered stereoframe. Return a struct 'state' which contains various useful bits of information for 3D stereoscopic rendering of a scene, based on head tracking data
        oc.HMD = [oc.HMD; state.modelView{1}];
    end

    if pa.trialNumber>track_trial %dont' update head position during a trial
       positions = -ds.floorWidth/2+2*ds.floorWidth/2*rand(2,pa.nball); %uniform random positions across floor
       positions(2,:) = positions(2,:)-ds.floorWidth/2;
       exclude = [0;-ds.floorWidth/2]; %x,z coordinate where fixation target is
       while sum(vecnorm(positions-exclude)<pa.fixationSize+pa.paddleHalfWidth)>0 %at least one position overlaps with fixation
           idx = find(vecnorm(positions-exclude)<pa.fixationSize+pa.paddleHalfWidth);
           positions(:,idx) = -ds.floorWidth/2+2*ds.floorWidth/2*rand(2,length(idx));
       end
       
       if pa.trialNumber>1
           eye = originaleye;
           
       else
           eye = PsychVRHMD('GetEyePose', ds.hmd, ds.renderPass, ds.globalHeadPose);
           R = [1 0 0; 0 cos(pa.gazeangle) -sin(pa.gazeangle); 0 sin(pa.gazeangle) cos(pa.gazeangle)];
           eye.modelView = [1 0 0 0; 0 1 0 0; 0 0 1 -ds.viewingDistance; 0 0 0 1];
           eye.modelView(1:3,1:3) = eye.modelView(1:3,1:3)*R; 
           originaleye = eye;
       end

       track_trial = track_trial+1;
    end
    
    % Render the scene separately for each eye:
    for renderPass = 0:1 %0 left, 1 right eye
        ds.renderPass = renderPass;
        
        if isempty(ds.hmd) % hmd not connected
            % Get head position from keyboard input
            
            eye.eyeIndex = ds.renderPass; % We are switching eye index, but not the eye.modelview here

            % CSB: this inits camera matrix if it is empty
            % BR: This should happen outside of the loop init'd on line 131
            if ~isfield(eye,'modelView') % checks if camera matrix exists. it shouldn't on first trial and it will be init'd.
                if ds.renderPass==0 % drawing left eye
                    eye.modelView = oc.defaultState.modelViewDataLeft;
                elseif ds.renderPass==1 % drawing right eye
                    eye.modelView =  oc.defaultState.modelViewDataRight;
%                     eye.modelView(1,4) =  eye.modelView(1,4)+100; % CSB: debug
                end
            end
            
            [pa, kb, eye] = GetKeyboardHeadmotion(pa, ds, kb, eye);  % query the keyboard to allow the observer to rotate the paddle and eventually lock in his/her response to initiate a new trial

        else % hmd connected
            % Query which eye to render in this ds.renderPass, and query its
            % eyePose vector for the predicted eye position to use for the virtual
            % camera rendering that eyes view. The returned pose vector actually
            % describes tracked head pose, ie. HMD position and orientation in space.
            

            [pa, kb, eye] = GetKeyboardHeadmotion(pa,ds,kb,eye);
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
        

        glClearColor(.5,.5,.5,1); % gray background
        glClear(); % clear the buffers - must be done for every frame
        glColor3f(1,1,1);
        
        glPushMatrix;
        
        % fixation target
        glPushMatrix;
        glTranslatef(0,pa.floorHeight+pa.fixationSize,-ds.floorWidth/2);
        glCallList(ds.highcontrastTarget);
        glPopMatrix;
        
       
          
            
         %% Experiment Logic  
        if ds.vbl <  pa.trialOnset + pa.targetMotionDuration % if current time < present until time: draw target, 1 s target motion

            % Target position 
            % make the onset delayed
            delay = 0.5;
            if ds.vbl <  pa.trialOnset + delay
                t = 0;
            else
                t = ds.vbl-pa.trialOnset-delay;
            end
            xPosition =pa.xSpeed.*t;  %pa.xSpeed.*t;
            zPosition = pa.zSpeed.*t;  %pa.zSpeed.*t
            
            glPushMatrix;
            glTranslatef(xPosition+.5*(pa.LR(pa.trialNumber)),pa.floorHeight+pa.paddleHalfHeight,zPosition-ds.floorWidth/2); % shift the target to its position along its trajectory for this frame
%             pa.xPosition = [pa.xPosition, xPosition+.5*(pa.LR(pa.trialNumber))];
%             pa.zPosition = [pa.zPosition, zPosition-ds.floorWidth/2];
%             if pa.targetContrast==1
                glCallList(ds.paddleList);
                glPopMatrix;
%             elseif pa.targetContrast==0.15
%                 glCallList(ds.midcontrastTarget);
%             elseif pa.targetContrast==0.075
%                 glCallList(ds.lowcontrastTarget); 
%             end
%              
            
            % stationary object
            glPushMatrix;
            glTranslatef(.5*(-pa.LR(pa.trialNumber)),pa.floorHeight+pa.paddleHalfHeight,-ds.floorWidth/2); 
            glCallList(ds.paddleList);
            glPopMatrix;
            
            
            % place random stationary objects
            for b = 1:pa.nball
            	glPushMatrix;
                glTranslatef(positions(1,b),pa.floorHeight+pa.paddleHalfHeight,positions(2,b)); 
                glCallList(ds.paddleList);
                glPopMatrix;
            end
            
                
%             glBindTexture(GL.TEXTURE_2D,ds.floor_texid);
%             glCallList(ds.floorTexture);

            eye.modelView(3,4) =  pa.translation.*(ds.vbl-pa.trialOnset); %eye.modelView(3,4) + %eye.modelView(2,4) - %(ds.vbl-pa.trialOnset)
            track = [track (ds.vbl-pa.trialOnset)];
            eye.modelView(2,4) = - pa.translation.*(ds.vbl-pa.trialOnset)*sin(pa.gazeangle); %adjust so translation is strictly foward

            pa.responseOnset = ds.vbl; % start the timer on the response time
            

        elseif ~kb.responseGiven % show paddle and allow observers to adjust its position - no time constraints - they press the space bar to lock in their response and start a new trial
            %2, % debugging flag
            
%             glPushMatrix;
%             % here's where we draw the paddle
%             glRotatef(pa.paddleAngle(pa.thisTrial),0,-1,0); % position the paddle along the orbit according to the angle specified by the input (observer adjustment)
%             glTranslatef(pa.paddleOrbitShift-pa.paddleHalfWidth*2,pa.floorHeight+pa.paddleHalfHeight,0); % shift the paddle physically out to the orbit
%             glCallList(ds.paddleList); % call the pre-compiled list that binds the textures to the paddle faces
%             glPopMatrix;
            
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
            
        elseif kb.responseGiven && pa.feedbackFlag==2 && pa.feedbackGiven==0 && ds.vbl < (pa.feedbackOnset + pa.timeToPaddle - pa.targetMotionDuration) % visual feedback
            
            %3, % debugging flag
            
            % show the feedback (rest of the trajectory)
            % Draw Paddle
            glPushMatrix;
            glRotatef(pa.paddleAngle(pa.thisTrial),0,-1,0); % position the paddle along the orbit according to the angle specified by the input (observer adjustment)
            glTranslatef(pa.paddleOrbitShift-pa.paddleHalfWidth*2-.5,pa.floorHeight+pa.paddleHalfHeight,0); % shift the paddle physically out to the orbit
            glCallList(ds.paddleList); % call the pre-compiled list that binds the textures to the paddle faces
            glPopMatrix;
            
            % In this section, we are display both an opaque paddle + a
            % translucent target (sphere) so we need to do a bit of
            % openGL work to make the translucent sphere display
            % properly
            glDepthMask(GL.FALSE); % disable changes to the depth buffer
            glBlendFunc(GL.SRC_ALPHA, GL.ONE); % set the alpha to that of the target contrast without influence of the paddle or other scene elements
            
            % Target position - we only speed up the feedback if it's
            % going to take longer than 10 s and we speed it up
            % incrementally given time
            xPosition = (pa.speedUpFlag*pa.xSpeed).*(ds.vbl-pa.feedbackOnset+pa.targetMotionDuration);
            zPosition = (pa.speedUpFlag*pa.zSpeed).*(ds.vbl-pa.feedbackOnset+pa.targetMotionDuration);
            
            glTranslatef(xPosition-.5,pa.floorHeight+pa.targetSize,zPosition); % shift the target to its position along its trajectory for this frame
            
            if pa.targetContrast==1
                glCallList(ds.highcontrastTarget);
            elseif pa.targetContrast==0.15
                glCallList(ds.midcontrastTarget);
            elseif pa.targetContrast==0.075
                glCallList(ds.lowcontrastTarget);
            end
            
            glDepthMask(GL.TRUE); % resume the ability to make changes to the depth buffer for proper rendering of remaining components
            glBlendFunc(GL.SRC_ALPHA, GL.ONE_MINUS_SRC_ALPHA); % restore the proper blending function
  
            
             
        elseif (kb.responseGiven && pa.feedbackFlag==1 && pa.feedbackGiven==0) || (kb.responseGiven && pa.feedbackFlag==2 && pa.feedbackGiven==0) % sound feedback only
            %4, % debugging flag
            
            reportedAngle = angle(exp(1i*deg2rad(pa.paddleAngle(pa.thisTrial))));  % converts the paddle angle from 0-360 to -pi - + pi
            
            % display appropriate object
            if ~pa.feedbackGiven
%                 if pa.LRresponse(pa.trialNumber) == pa.LR(pa.trialNumber)  %does the response match the target side
%                         glPushMatrix;  
%                         glTranslatef(0,pa.floorHeight+pa.targetSize,-ds.floorWidth/2-0.5); % display green sphere
%                         glCallList(ds.correct);
%                         glPopMatrix;
%                 else
%                         glPushMatrix;
%                         glTranslatef(0,pa.floorHeight+pa.targetSize,-ds.floorWidth/2-0.5); % display red sphere
%                         glCallList(ds.incorrect);
%                         glPopMatrix;
%                 end
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
            
            
        elseif (kb.responseGiven && pa.feedbackFlag==1 && pa.feedbackGiven==1 && ds.vbl <= pa.waitTime+1) || (kb.responseGiven && pa.feedbackFlag==2 && pa.feedbackGiven==1 && ds.vbl <= pa.waitTime+.75)
            
%             % Draw Paddle
%             glPushMatrix; 
%             glRotatef(pa.paddleAngle(pa.thisTrial),0,-1,0); % position the paddle along the orbit according to the angle specified by the input (observer adjustment)
%             glTranslatef(pa.paddleOrbitShift-pa.paddleHalfWidth*2,0,0); % shift the paddle physically out to the orbit
%             %                 if renderPass == pa.randeye(pa.trialNumber)
%             glCallList(ds.paddleList); % call the pre-compiled list that binds the textures to the paddle faces
%             %                 end
%             glPopMatrix;
%             
%             % In this section, we are display both an opaque paddle + a
%             % translucent target (sphere) so we need to do a bit of
%             % openGL work to make the translucent sphere display
%             % properly
%             glDepthMask(GL.FALSE); % disable changes to the depth buffer
%             glBlendFunc(GL.SRC_ALPHA, GL.ONE); % set the alpha to that of the target contrast without influence of the paddle or other scene elements
%             glTranslatef(xPosition,0,zPosition); % shift the target to its position along its trajectory for this frame
%             %                 if renderPass == pa.randeye(pa.trialNumber)
%             if pa.feedbackFlag==2
%                 if pa.targetContrast==1
%                     glCallList(ds.highcontrastTarget);
%                 elseif pa.targetContrast==0.15
%                     glCallList(ds.midcontrastTarget);
%                 elseif pa.targetContrast==0.075
%                     glCallList(ds.lowcontrastTarget);
%                 end
%             end
%             glDepthMask(GL.TRUE); % resume the ability to make changes to the depth buffer for proper rendering of remaining components
%             glBlendFunc(GL.SRC_ALPHA, GL.ONE_MINUS_SRC_ALPHA); % re-store the proper blending function
%             

                if pa.LRresponse(pa.trialNumber) == pa.LR(pa.trialNumber)  %does the response match the target side
                        glPushMatrix;
                        glTranslatef(0,pa.floorHeight+pa.fixationSize,-(ds.floorWidth/2)); % display green sphere
                        glCallList(ds.correct);
                        glPopMatrix;
                else
                        glPushMatrix;
                        glTranslatef(0,pa.floorHeight+pa.fixationSize,-(ds.floorWidth/2)); % display red sphere
                        glCallList(ds.incorrect);
                        glPopMatrix;
                end
                
        elseif (kb.responseGiven && pa.feedbackFlag==0) || (kb.responseGiven && pa.feedbackFlag==1 && pa.feedbackGiven==1) || (kb.responseGiven && pa.feedbackFlag==2 && pa.feedbackGiven==1) % done, set up for the next trial (i.e., determine the new random trajectory)
            %                 5, % debugging flag
            
            [ds, pa, kb, oc] = SetupNewTrial(ds, pa, kb, oc);
            track_trial = pa.trialNumber-1;
            
        end
        %% End of experiment logic
        
        % Update view with keyboard responses
        glMatrixMode(GL.MODELVIEW); 
        glLoadMatrixd(modelView);
        
        glPopMatrix; % reset back to the origin
        moglDrawDots3D(ds.w, [inf inf inf], 10, [1 1 1 1], [], 2); % 'hack' to fix transparency issue - same one we've used in all 3D pong experiments - definitely works!
        
        % Distance debugging code
%         glCallList(ds.debugcylinder); % call the pre-compiled list that binds the textures to the paddle faces
        
%         glBindTexture(GL.TEXTURE_2D,ds.roomwall_texid); % was suggested to bind textures before/outside of call lists rather than in - doesn't buy us anything from what I can tell though
%         glCallList(ds.wallTexture);
%         glBindTexture(GL.TEXTURE_2D,ds.ceiling_texid);
%         glCallList(ds.ceilingTexture);
        
            glBindTexture(GL.TEXTURE_2D,ds.floor_texid);
            glCallList(ds.floorTexture);
        
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

pa.dataFile = fullfile(pa.baseDir, 'Data', [pa.subjectName '-' ds.experimentType '-' pa.date '-' num2str(pa.block) '.mat']);
save(pa.dataFile, 'pa', 'ds', 'kb','oc');

% Calculate average framerate:
fps = ds.fCount / (ds.vbl - tStart), % uncomment to print out at end of run
% Done (or quit out). Save data (in pa.response) and other relevant parameters/variables, close screen and exit
PsychVRHMD('SetAutoClose', ds.hmd, 2);
PsychPortAudio('Close', pa.handleHit);
PsychPortAudio('Close', pa.handleMiss);
Priority(0);
ShowCursor(ds.screenId);
ListenChar(1); % Let's start listening to the keyboard again
sca;
