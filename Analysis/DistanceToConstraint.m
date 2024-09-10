function  [dconstraint, dsurround] = DistanceToConstraint(ds, pa, depth_range, theta, trial)

% simulate VR flow scene to generate distance to constraint for multiple velocities
% takes saved variables from VR experiment
% plots center, surround velocities as well as constraint
visualize = 0;
% seed=2;
% rng(seed) % to have random dots that appear in the same "random" place each time
ns = pa.targetMotionDuration; % number of seconds
world_speed = pa.translation; % m/s speed of the observer in a straight line
if ~isfield(ds,'frameRate')
    fps = 90;
else
    fps = ds.frameRate; %Screen(screenNumber,'FrameRate'); % should be 144 to match experiment
end

 
height = -pa.floorHeight;

speeds =  pa.speed;
directions = pa.direction ; 


if nargin < 4 || isempty(trial)
    theta = [];
    conditions = fullfact([numel(speeds), numel(directions)]);
    condition_indices = conditions;
    conditions = [speeds(conditions(:,1)).*cos(directions(conditions(:,2))); speeds(conditions(:,1)).*sin(directions(conditions(:,2)))]';
else
    condition_speeddir = [pa.fullFactorial(trial,3), pa.fullFactorial(trial,4)];
    conditions = [pa.fullFactorial(trial,1), pa.fullFactorial(trial,2)];
end

if ~isfield(pa,'objectdist')
    object_dist = 2;
    fixation = 3;
else
    object_dist = pa.objectdist;
    fixation = pa.fixationdist;
end

if depth_range == 0
    calculate_segment = 0; % calculate to segment or point (0)
else
    calculate_segment = 1;
end

dim = [pa.floorWidth,0,pa.floorWidth]; % extent of where dots can be in m: X, Y, Z. Depth is more than how far you're travelling (ns *speed) + a little extra 
% 5 m across
nClusters = 1000; % specify number of clusters
nDotsPerCluster = 1;% number of dots per cluster
nObjects = pa.nball;
 
% *** find these based on oculus display
view_dist = .5; %m how far the screen is from the observer
viewingdepths = [.01,   5]; % nearest and furthest dots that can show up, m
 
windowRect =  [0  0  2560 1600]; % [0  0  2560 1600] [0  0  2160 1200 ]switch with ds.windowRect screen size in pixels (origin, width of screen, height of screen) [0  0  2560 1600]
pixels = windowRect(3:4); % pixel width and height of screen
screensize = [.712 .312]; % screen size in m
 
ppcm = [39,39]; % pixel per cm of the screen
xyrat = windowRect(3)/windowRect(4);
 
clusters = rand(nClusters,3).*dim - [.5*dim(1), .5*dim(2)-height, 0]; % randomize position of dots, centered around x = 0, and ground pushed down
% [.5*dim(1), -height, 0];
 
%object positions
positions = -dim(1)/2+2*dim(1)/2*rand(nObjects,2); %uniform random positions across floor
positions(:,2) = positions(:,2)+fixation;
% positions = -[pa.positions(1,:)' pa.positions(3,:)'];
% positions(3,:) = -positions(3,:);
 
dots = repmat(clusters,nDotsPerCluster,1); % ground plane
% for no floor
% dots = [];
objectwidth = pa.paddleHalfWidth;
object = [objectwidth, objectwidth, objectwidth]; %length, width, height
dotsperobj = 15;
a = -object(1);
b = object(1);
aboveground = -pa.aboveground;
 
if ~isempty(nObjects)
    for obj = 1:nObjects
        r = (b-a).*rand(dotsperobj,3) + a; % in line below added +pa.positions(2,obj)'
        newpositions = [r(:,1)+positions(obj,1), r(:,2)+(height-b)-pa.positions(2,obj)', r(:,3)+positions(obj,2)];
        dots = [dots; newpositions];
    end
    
%     y = pa.positions(2,:)+height-pa.paddleHalfHeight;


end
            
fixation_dot = [0, height, fixation];
dots(end+1,:) = fixation_dot;
fixation_idx = length(dots);


% stationary and target positions
% stationary on left and moving on right to demonstrate simulation with
% correct direction degree labels - left hand side object always gets
% mirror added motion
stationary_target = [-0.5, aboveground+height-b, object_dist; 0.5, aboveground+height-b, object_dist];
% stationary then target
 
for obj = 1:2 %stationary obj and moving obj
    r = (b-a).*rand(dotsperobj,3) + a;
    newpositions = [r(:,1)+stationary_target(obj,1), r(:,2)+stationary_target(obj,2), r(:,3)+stationary_target(obj,3)];
    dots = [dots; newpositions];
end
 
nDots = length(dots); % total numbber of dots
stationary_idx = (length(dots)+1-2*dotsperobj):length(dots)-dotsperobj;
target_idx = (length(dots)+1-dotsperobj):length(dots);

start_dots = dots;
 
 
%visualize dots, orient so that Z axis extends from observer to direction
%of gaze
if visualize
    figure(1), scatter3(dots(:,3), dots(:,1), -dots(:,2), 'filled')
    hold on, scatter3(dots(end-2,3), dots(end-2,1), -dots(end-2,2), 'filled', 'r')
 
    xlabel('Z')
    ylabel('X')
    zlabel('Y')
    axis equal
end
 
 
% Observer trajectory
 
% create matrix where each row is a velocity vector, specify velocity
% between each frame
 
% first specify position over time
% straight line, don't change X or Y position, just Z coordinate
% start at 0, end at some distance (defined by ns*speed)
% speed/fps defines how far to go over each frame 
trajectory = [zeros(ns*fps+1,1), zeros(ns*fps+1,1),(0:(world_speed/fps):(ns*world_speed))'];
 
% loop over speeds and directions for object
% dconstraint = NaN(numel(speeds), numel(directions));
% dsurround = NaN(numel(speeds), numel(directions));
dconstraint = NaN(size(conditions,1),1);
dsurround = NaN(size(conditions,1),1);

% figure('Position', [10 10 1200 600])

for cond = 1:size(conditions, 1)
    
    % trajectory = [zeros(ns*fps+1,1), 0.2*sin(0:(speed/fps):(ns*speed))', (0:(speed/fps):(ns*speed))'];
    %x-z plane
    %target_trajectory = [sign(stationary_target(2,1))*speeds(conditions(cond,1))*cos(directions(conditions(cond,2)))*(0:1/fps:ns)', zeros(ns*fps+1,1), speeds(conditions(cond,1))*sin(directions(conditions(cond,2)))*(0:1/fps:ns)'];
    target_trajectory = [sign(stationary_target(2,1))*conditions(cond,1)*(0:1/fps:ns)', zeros(ns*fps+1,1), conditions(cond,2)*(0:1/fps:ns)'];
    
    target_trajectory = target_trajectory + stationary_target(2,:);
    
    % rotation
    secs = 1/fps:1/fps:ns;
    if isempty(theta)
        theta = atan(height./(fixation-world_speed*secs)); % update theta for observer fixating at a point at the ground in front of them, fixation m away
    end
    
    % visualize observer trajectory within dots
    if visualize
        figure(1), hold on, plot3(trajectory(:,3), -trajectory(:,1), -trajectory(:,2), 'LineWidth', 3)
        hold on, plot3(target_trajectory(:,3), -target_trajectory(:,1), -target_trajectory(:,2), 'LineWidth', 2)
        legend('dots', 'trajectory', 'moving object')
        title('environment and observer trajectory')
    end
    
    % calculate velocity between frames
    v = diff(trajectory);
    v_target = diff(target_trajectory);
    T = NaN(size(v));
    v_constraint = nan(2,nDots,ns*fps);
    v_constraint_far = nan(2,nDots,ns*fps);
    v_constraint_close = nan(2,nDots,ns*fps);
    
    % holder matrices for screen positions
    x = nan(nDots,ns*fps);
    y = nan(nDots,ns*fps);
    Z = nan(nDots,ns*fps);
    I = true(nDots,ns*fps);
    
    drawndots = NaN([size(dots) ns*fps]);
    for ii=1:ns*fps %
        velocity = v(ii,:); % how much did the observer move
        t_vel = v_target(ii,:);
        % moving observer = moving world relative to observer in equal and opposite way
        % recalculating world coordinates in terms of observer reference frame,
        % where observer is always at the origin
        if ii == 1
            dots = start_dots - velocity;
        else
            dots = dots - velocity; %shift dots in world coordinates
        end
        dots(target_idx,:) = dots(target_idx:end,:)+t_vel; %add velocity to moving object
        
        % if the observer rotates, rotate the world based on 3D rotation matrix
        observerRotation = [1, 0, 0; 0, cos(theta(ii)), -sin(theta(ii)); 0, sin(theta(ii)), cos(theta(ii))];
        drawndots(:,:,ii) = (observerRotation*dots')'; %observer coordinate dots
        
        %for calculating velocity based on constraint equation
        Z(:,ii) = drawndots(:,3,ii); %Z value of points
        T(ii,:) = (observerRotation*velocity')'; %rotate velocity vector as well
        
        % Using projective geometry (similar triangles) to calculate where on
        % the screen the dots should appear
        
        % x = x coordinate on screen, y = y coordinate on screen, convert to cm
        x(:,ii) = 100*view_dist*(drawndots(:,1,ii))./(drawndots(:,3,ii));
        y(:,ii) = 100*view_dist*(drawndots(:,2,ii))./(drawndots(:,3,ii));
        
        % calculate velocity based on constraint eq, make sure x,y in m
        v_constraint(:,:,ii) = constraint_velocity_screen(Z(:,ii), [x(:,ii)';y(:,ii)']./100, T(ii,:), view_dist,Z(fixation_idx,ii));
        v_constraint(:,:,ii) = v_constraint(:,:,ii).*100;
        
        % smaller range for far/close on constraint line
        v_constraint_far(:,:,ii) = constraint_velocity_screen(ones(size(Z(:,ii))).*(Z(:,ii)+depth_range), [x(:,ii)';y(:,ii)']./100, T(ii,:), view_dist,Z(fixation_idx,ii));
        v_constraint_far(:,:,ii) = v_constraint_far(:,:,ii)*100;
        v_constraint_close(:,:,ii) = constraint_velocity_screen(ones(size(Z(:,ii))).*(Z(:,ii))-depth_range, [x(:,ii)';y(:,ii)']./100, T(ii,:), view_dist,Z(fixation_idx,ii));
        v_constraint_close(:,:,ii) = v_constraint_close(:,:,ii)*100;
        
        % Indices of dots to show based on how close/far the dots in the real world are (viewing depths)
        I(:,ii) = drawndots(:,3,ii) > viewingdepths(1)...
            & drawndots(:,3,ii)< viewingdepths(2);
        % and the screensize
        I(:,ii) = I(:,ii) & abs(x(:,ii)*ppcm(1))<windowRect(3)/2 & abs(y(:,ii)*ppcm(2))<windowRect(4)/2;
        
        
    end
    
    rvelocityX = diff(x,1,2);
    rvelocityY = diff(y,1,2);
    
%     %visualize first frame in pixels
    if visualize
        figure, scatter(x(I(:,1),1)*ppcm(1), -y(I(:,1),1)*ppcm(2), 'filled')
        hold on, scatter(x(fixation_idx,1)*ppcm(1), -y(fixation_idx,1)*ppcm(2), 'filled', 'r') %fixation
        hold on, scatter(x(stationary_idx,1)*ppcm(1), -y(stationary_idx,1)*ppcm(2), 'filled', 'b') %stationary
        hold on, scatter(x(target_idx,1)*ppcm(1), -y(target_idx,1)*ppcm(2), 'filled', 'g') %target
        xlim([-windowRect(3)/2, windowRect(3)/2])
        ylim([-windowRect(4)/2, windowRect(4)/2])
        axis equal
        title('first frame')
    end
%     
    
    % surround velocities
    % calculate in terms of degrees
    degX = atand(drawndots(:,1,:)./drawndots(:,3,:));
    degY = atand(drawndots(:,2,:)./drawndots(:,3,:));
    
    
    % calculate velcoities in deg/s
    
    rvXdeg = atand(rvelocityX/100*view_dist./(view_dist^2+(rvelocityX/100+x(:,1:end-1).*x(:,1:end-1))));
    rvYdeg = atand(rvelocityY/100*view_dist./(view_dist^2+(rvelocityY/100+y(:,1:end-1).*y(:,1:end-1))));
    
    
    %show target vs surround velocities throughout stim
    radius = 3; %in cm
    center = target_idx; %target_idx vs stationary_idx
    xlims = [-.1, .1];
    ylims = [-.1, .1];
    
    
    
    mean_d = NaN(ns*fps-1,2);
    for ii = 1:ns*fps-1
        
        onscreen = find(I(:,ii));
        
        [ti, target_onscreen, tb] = intersect(onscreen, target_idx);
        
        if isempty(target_onscreen)
            % don't change mean_d
        else
            d = NaN(1, length(target_onscreen));
            
            for jj = 1:length(target_onscreen)
                
                % calcualte distance to constraint
                if calculate_segment
                    d(jj) = point2segment([rvelocityX(onscreen(target_onscreen(jj)),ii); rvelocityY(onscreen(target_onscreen(jj)),ii)], v_constraint_far(:,onscreen(target_onscreen(jj)),ii), v_constraint_close(:,onscreen(target_onscreen(jj)),ii));
                else
                    dif = [rvelocityX(onscreen(target_onscreen(jj)),ii); rvelocityY(onscreen(target_onscreen(jj)),ii)] - v_constraint(:,onscreen(target_onscreen(jj)),ii);
                    d(jj) = norm(dif);
                end
            end
            mean_d(ii,1) = mean(d);
            
            % calculate distance to surround
            center_point = [mean([max(x(onscreen(target_onscreen),ii)),min(x(onscreen(target_onscreen),ii))]), mean([max(y(onscreen(target_onscreen),ii)),min(y(onscreen(target_onscreen),ii))])];
            distance2center_point = vecnorm((center_point - [x(onscreen,ii),y(onscreen,ii)])');
            window_idx = find(distance2center_point<radius);
            % plot window on object
            % hold on, scatter(x(window_idx,ii), -y(window_idx,ii), 50,[0.8500 0.3250 0.0980])
            
            % find non target velocities
            surround_idx = window_idx(~ismember(window_idx, target_onscreen));
            % hold on, scatter(x(surround_idx,ii), -y(surround_idx,ii), 50,[0 0.4470 0.7410])
            
            % get target and surround velocity mean
            center_mean= mean([rvelocityX(onscreen(target_onscreen),ii), rvelocityY(onscreen(target_onscreen),ii)]);
            surround_mean = mean([rvelocityX(onscreen(surround_idx),ii), rvelocityY(onscreen(surround_idx),ii)],1);
            mean_d(ii,2) = vecnorm(center_mean-surround_mean);
            
            
        end
        % plot mean velocity object and surround
        
        if ii == 1 %round((ns*fps-1)/2)
    
            if visualize
                subplot(numel(directions),numel(speeds),cond)
                clf

                set(gcf,'position',[500, 500, 600, 400])
                set(gcf,'color','w');
            
            
            % plot suround velocities
            surround_idx = onscreen(surround_idx);
            quiver(zeros(size(rvelocityX(surround_idx,ii))),zeros(size(rvelocityX(surround_idx,ii))), rvelocityX(surround_idx,ii), -rvelocityY(surround_idx,ii), 'AutoScale', 'off', 'LineWidth', 2)
            hold on
            quiver(zeros(size(rvelocityX(center,ii))),zeros(size(rvelocityX(center,ii))), rvelocityX(center,ii), -rvelocityY(center,ii), 'AutoScale', 'off', 'LineWidth', 2)
    
%            plot mean velocity object and surround
            if dotsperobj>1
                hold on
                quiver(0,0, center_mean(1), -center_mean(2), 'r','AutoScale', 'off', 'LineWidth', 5)
                hold on
                quiver(0,0,surround_mean(1), -surround_mean(2), 'color',[0,0,0.75],'AutoScale', 'off', 'LineWidth', 5)
            end
    
            hold on,
            for jj = 1:length(center)
                plot([v_constraint_close(1,center(jj), ii) v_constraint_far(1,center(jj), ii)], -[v_constraint_close(2,center(jj),ii) v_constraint_far(2,center(jj),ii)], 'k', 'LineWidth', 2)
    
            end
            axis equal
            xlim(xlims)
            ylim(ylims)
            if exist('conditon_indices','var')
                title(['s = ', num2str(speeds(condition_indices(cond,1))), ' m/s,  dir = ', num2str(rad2deg(directions(condition_indices(cond,2)))), ' deg'])
            else
                title(['s = ', num2str(condition_speeddir(1)), ' m/s,  dir = ', num2str(rad2deg(condition_speeddir(2))), ' deg'])

            end
            pause(1)
            end
            %
        end
    end
    
    dconstraint(cond) = nanmean(mean_d(:,1));
    dsurround(cond) = nanmean(mean_d(:,2));
    
end


