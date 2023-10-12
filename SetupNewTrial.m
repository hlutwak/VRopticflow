function [ds, pa, kb, oc] = SetupNewTrial(ds, pa, kb, oc)

% add device to argument if sending events

% do this before every trial to get the next trial's starting info (e.g.,
% target motion parameters)


if pa.trialNumber>0 % if it's past the first trial, wait for the Up Arrow key to be pressed after the feedback to initiate a new trial
    [kb.keyIsDown, kb.secs, kb.keyCode] = KbCheck(-1);
    
    if kb.keyIsDown
        
        if kb.keyCode(kb.reorientKey)
  
            pa.trialNumber = pa.trialNumber + 1;
            
            pa.thisTrial = pa.trialNumber;

            
            if pa.trialNumber <= pa.nTrials
                
                pa.positions = -pa.floorWidth/2+2*pa.floorWidth/2*rand(3,pa.nball); %uniform random positions across floor, random height within range of y position of moving target
                pa.positions(3,:) = pa.positions(3,:)-pa.floorWidth/2;
                pa.positions(2,:) = pa.aboveground*rand(1,pa.nball);
                fixposition = [0; pa.floorHeight;-pa.fixationdist]; %x,z coordinate where fixation target is
                exclude = [0 -.5 0.5; pa.floorHeight pa.floorHeight pa.floorHeight; -pa.fixationdist -pa.objectdist -pa.objectdist];
% 
                    for eachpoint = 1:3
                        withinexclude = sum(vecnorm(pa.positions-exclude(:,eachpoint))<pa.fixationSize+pa.paddleHalfWidth*2);
    % 
                        while withinexclude>0 %at least one position overlaps with fixation
                            idx = find(vecnorm(pa.positions-exclude(:,eachpoint))<pa.fixationSize+pa.paddleHalfWidth);
                            pa.positions(:,idx) = -pa.floorWidth/2+2*pa.floorWidth/2*rand(2,length(idx));
                        end
    %                     
                    end

                
                rotangle = [0, 90, 180, 270];
                axs = [1 0 0; 0 1 0; 0 0 1];
                pa.rotations = [rotangle(randi(length(rotangle),pa.nball,1))' axs(randi(3,pa.nball,1),:)];
                
                
                if ds.dotfield
                    pa.dotpositions = -pa.floorWidth/2+2*pa.floorWidth/2*rand(3,pa.nball*5); %uniform random positions across floor, random height within range of y position of moving target
                    pa.dotpositions(3,:) = pa.dotpositions(3,:)-pa.floorWidth/2;
                    pa.dotpositions(2,:) = pa.floorHeight;
                    
                    if ds.control
                        surround = [-1.5, pa.floorHeight, -(pa.objectdist+2.25); 1.5, pa.floorHeight, -(pa.objectdist+2.25)]';
                        surround_idx = find(vecnorm(pa.dotpositions - surround(:,1))<2 | vecnorm(pa.dotpositions - surround(:,2))<2);
                        pa.surroundpositions = pa.dotpositions(:,surround_idx);
                        pa.nball = length(surround_idx);
                    end
                end
                
                pa.xSpeed = pa.fullFactorial(pa.thisTrial,1);
                pa.zSpeed = pa.fullFactorial(pa.thisTrial,2);
                
                pa.timeToPaddle = (pa.paddleOrbitShift-pa.paddleHalfWidth*3.8) ./ sqrt(pa.xSpeed.^2 + pa.zSpeed.^2);
                
                % speed up the really really slow feedback - this happens
                % on maybe 10% of trials?
                if pa.timeToPaddle > 2
                    pa.speedUpFlag = pa.timeToPaddle / 2; % just make the feedback no longer than 2 s long
                else
                    pa.speedUpFlag = 1;
                end
                pa.timeToPaddle = (pa.paddleOrbitShift-pa.paddleHalfWidth*3.8) ./ sqrt((pa.speedUpFlag*pa.xSpeed).^2 + (pa.speedUpFlag*pa.zSpeed).^2);
                
                if pa.trialNumber==1
                    pa.paddleAngle(pa.thisTrial) = pa.paddleAngle; % deg
                    pa.paddleAngleInitial = pa.paddleAngle(pa.thisTrial);
                else
                    pa.paddleAngle(pa.thisTrial) = randi(361)-1; % deg  - assigns a random integer angle between 0 and 360
                    pa.paddleAngleInitial = pa.paddleAngle(pa.thisTrial);
                end
                
                pa.targetContrast = pa.fullFactorial(pa.thisTrial,1);
                
                kb.responseGiven = 0;
                pa.feedbackGiven = 0;
    
                pa.trialOnset = ds.vbl; % the trial starts....NOW
                oc.trialStart = [oc.trialStart ds.vbl];
                oc.trial_startflag = 0;
                oc.trialendflag = 0;
%                 oc.UTCtrialStart = [oc.UTCtrialStart datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss:ms')];
%                 if ds.eyetracking
%                     device.send_event("trialStart")
%                 end

                
            end
        end
    end
    
else % if it is just the first trial, start right up after the subject presses 'space' to confirm that they are ready to go
    
    pa.trialNumber = pa.trialNumber + 1;
    
    pa.thisTrial = pa.trialNumber;

    
    if pa.trialNumber <= pa.nTrials
        
        pa.xSpeed = pa.fullFactorial(pa.thisTrial,1);
        pa.zSpeed = pa.fullFactorial(pa.thisTrial,2);
       
                pa.positions = -pa.floorWidth/2+2*pa.floorWidth/2*rand(3,pa.nball); %uniform random positions across floor, random height within range of y position of moving target
                pa.positions(3,:) = pa.positions(3,:)-pa.floorWidth/2;
                pa.positions(2,:) = pa.aboveground*rand(1,pa.nball);
                exclude = [0; pa.floorHeight;-pa.floorWidth/2]; %x,z coordinate where fixation target is
                while sum(vecnorm(pa.positions-exclude)<pa.fixationSize+pa.paddleHalfWidth)>0 %at least one position overlaps with fixation
                    idx = find(vecnorm(pa.positions-exclude)<pa.fixationSize+pa.paddleHalfWidth);
                    pa.positions(:,idx) = -pa.floorWidth/2+2*pa.floorWidth/2*rand(2,length(idx));
                end
                
                rotangle = [0, 90, 180, 270];
                axs = [1 0 0; 0 1 0; 0 0 1];
                pa.rotations = [rotangle(randi(length(rotangle),pa.nball,1))' axs(randi(3,pa.nball,1),:)];
                
                
                if ds.dotfield
                    pa.dotpositions = -pa.floorWidth/2+2*pa.floorWidth/2*rand(3,pa.nball*5); %uniform random positions across floor, random height within range of y position of moving target
                    pa.dotpositions(3,:) = pa.dotpositions(3,:)-pa.floorWidth/2;
                    pa.dotpositions(2,:) = pa.floorHeight;
                    
                    if ds.control
                        surround = [-1.5, pa.floorHeight, -(pa.objectdist+2.25); 1.5, pa.floorHeight, -(pa.objectdist+2.25)]';
                        surround_idx = find(vecnorm(pa.dotpositions - surround(:,1))<2 | vecnorm(pa.dotpositions - surround(:,2))<2);
                        pa.surroundpositions = pa.dotpositions(:,surround_idx);
                        pa.nball = length(surround_idx);
                    end
                end
                
      
        pa.timeToPaddle = (pa.paddleOrbitShift-pa.paddleHalfWidth*3.8) ./ sqrt(pa.xSpeed.^2 + pa.zSpeed.^2);

                if pa.timeToPaddle > 2
                    pa.speedUpFlag = pa.timeToPaddle / 2; % just make the feedback no longer than 2 s long
                else
                    pa.speedUpFlag = 1;
                end
        
        pa.timeToPaddle = (pa.paddleOrbitShift-pa.paddleHalfWidth*3.8) ./ sqrt((pa.speedUpFlag*pa.xSpeed).^2 + (pa.speedUpFlag*pa.zSpeed).^2);
        
        
        pa.paddleAngle(pa.thisTrial) = pa.paddleAngle; % deg
        pa.paddleAngleInitial = pa.paddleAngle(pa.thisTrial);
        
        pa.targetContrast = pa.fullFactorial(pa.thisTrial,1);
        
        kb.responseGiven = 0;
        pa.feedbackGiven = 0;

        pa.trialOnset = ds.vbl; % the trial starts....NOW
        oc.trialStart = pa.trialOnset;
        oc.trial_startflag = 0;
        oc.trialendflag = 0;
%         oc.UTCtrialStart = datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss:ms');

    end
end

