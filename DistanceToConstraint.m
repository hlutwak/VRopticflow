function  D = DistanceConstraint(ds, pa, depth_range)

% simulate VR flow scene to generate distance to constraint for multiple velocities
% takes saved variables from VR experiment


seed=2;
rng(seed) % to have random dots that appear in the same "random" place each time
ns = pa.targetMotionDuration;
fps = ds.frameRate;
world_speed = pa.translation;
height = -pa.floorheight;

conditions = fullfact([numel(pa.speed), numel(pa.direction)]);
speeds = pa.speed;
directions = pa.direction;
dim = [pa.floorWidth 0 pa.floorWidth];
nObjects = 50;
view_dist = .5; %m how far the screen is from the observer

windowRect = [0           0        2560    1600]; % screen size in pixels (origin, width of screen, height of screen)
pixels = windowRect(3:4); % pixel width and height of screen
screensize = [.712 .312]; % screen size in m
 
ppcm = [39,39]; % pixel per cm of the screen ** change this to ds.pixelsPerM
xyrat = windowRect(3)/windowRect(4);

positions = -dim(1)/2+2*dim(1)/2*rand(nObjects,2); %uniform random positions across floor
positions(:,2) = positions(:,2)+dim(3)/2;
 
dots = repmat(clusters,nDotsPerCluster,1); % ground plane

object = [pa.paddleHalfWidth, pa.paddleHalfWidth, pa.paddleHalfWidth]; %length, width, height
dotsperobj = 15;
a = -object(1);
b = object(1);
aboveground = -pa.aboveground;%-.15;

if ~isempty(nObjects)
    for obj = 1:nObjects
        r = (b-a).*rand(dotsperobj,3) + a;
        newpositions = [r(:,1)+positions(obj,1), r(:,2)+(height-b), r(:,3)+positions(obj,2)];
        dots = [dots; newpositions];
    end
end
 
fixation_dot = [0, height, dim(3)/2];
dots(end+1,:) = fixation_dot;
fixation_idx = length(dots);


% old analysis where these werent preset
if ~isfield(pa, 'objectdist')
    pa.objectdist = 2;
    pa.fixationdist = 3;
end

stationary_target = [-0.5, aboveground+height-b, pa.objectdist; 0.5, aboveground+height-b, pa.objectdist];

for obj = 1:2 %stationary obj and moving obj
    r = (b-a).*rand(dotsperobj,3) + a;
    newpositions = [r(:,1)+stationary_target(obj,1), r(:,2)+stationary_target(obj,2), r(:,3)+stationary_target(obj,3)];
    dots = [dots; newpositions];
end

nDots = length(dots); % total numbber of dots
stationary_idx = (length(dots)+1-2*dotsperobj):length(dots)-dotsperobj;
target_idx = (length(dots)+1-dotsperobj):length(dots);

start_dots = dots;

trajectory = [zeros(ns*fps+1,1), zeros(ns*fps+1,1),(0:(world_speed/fps):(ns*world_speed))'];
 
% loop over speeds and directions for object
D = NaN(numel(pa.speed), numel(pa.direction));

for cond = 1:size(conditions, 1) 
    
    % trajectory = [zeros(ns*fps+1,1), 0.2*sin(0:(speed/fps):(ns*speed))', (0:(speed/fps):(ns*speed))'];
    %x-z plane
    target_trajectory = [sign(stationary_target(2,1))*speeds(conditions(cond,1))*cos(directions(conditions(cond,2)))*(0:1/fps:ns)', zeros(ns*fps+1,1), speeds(conditions(cond,1))*sin(directions(conditions(cond,2)))*(0:1/fps:ns)'];
    target_trajectory = target_trajectory + stationary_target(2,:);
    
    % rotation
    secs = 1/fps:1/fps:ns;
    theta = atan(height./(pa.fixationdist-world_speed*secs)); % update theta for observer fixating at a point at the ground in front of them, fixation m away
    
    
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
        v_constraint_close(:,:,ii) = constraint_velocity_screen(ones(size(Z(:,ii))).*(Z(:,ii)-depth_range), [x(:,ii)';y(:,ii)']./100, T(ii,:), view_dist,Z(fixation_idx,ii));
        v_constraint_close(:,:,ii) = v_constraint_close(:,:,ii)*100;
        
        % Indices of dots to show based on how close/far the dots in the real world are (viewing depths)
        I(:,ii) = drawndots(:,3,ii) > viewingdepths(1)...
            & drawndots(:,3,ii)< viewingdepths(2);
        % and the screensize
        I(:,ii) = I(:,ii) & abs(x(:,ii)*ppcm(1))<windowRect(3)/2 & abs(y(:,ii)*ppcm(2))<windowRect(4)/2;
        
        
    end
    
    rvelocityX = diff(x,1,2);
    rvelocityY = diff(y,1,2);
    
    % visualize first frame in pixels
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
    xlims = [-.01, .08];
    ylims = [-.025,.025];
    
    
%     subplot(length(speeds),length(directions), cond)
    figure(1)
    % set(gcf,'position',[500, 500, 600, 400])
    set(gcf,'color','w');
    
    ii = 1; %:ns*fps-1
        clf
        center_point = [mean([max(x(center,ii)),min(x(center,ii))]), mean([max(y(center,ii)),min(y(center,ii))])];
        distance2center_point = vecnorm((center_point - [x(:,ii),y(:,ii)])');
        window_idx = find(distance2center_point<radius);
        % plot window on object
        % hold on, scatter(x(window_idx,ii), -y(window_idx,ii), 50,[0.8500 0.3250 0.0980])
        
        % find non target velocities
        surround_idx = window_idx(~ismember(window_idx, center));
        % hold on, scatter(x(surround_idx,ii), -y(surround_idx,ii), 50,[0 0.4470 0.7410])
        
        % get target and surround velocity mean
        center_mean= mean([rvelocityX(center,ii), rvelocityY(center,ii)]);
        surround_mean = mean([rvelocityX(surround_idx,ii), rvelocityY(surround_idx,ii)],1);
        
        % plot suround velocities
        
        quiver(zeros(size(rvelocityX(surround_idx,ii))),zeros(size(rvelocityX(surround_idx,ii))), rvelocityX(surround_idx,ii), -rvelocityY(surround_idx,ii), 'AutoScale', 'off', 'LineWidth', 2)
        hold on
        quiver(zeros(size(rvelocityX(center,ii))),zeros(size(rvelocityX(center,ii))), rvelocityX(center,ii), -rvelocityY(center,ii), 'AutoScale', 'off', 'LineWidth', 2)
        
        % plot mean velocity object and surround
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
%         pause(1/fps)
    
    
    % calculating distance to constraint segment/point for just target
    %
    % figure
    % set(gcf,'position',[250, 250, 800, 400])
    % set(gcf,'color','w');
    
    mean_d = NaN(ns*fps-1,1);
    
    for ii = 1:ns*fps-1
        
        onscreen = find(I(:,ii));
        
        [ti, target_onscreen, tb] = intersect(onscreen, target_idx);
        
        if isempty(target_onscreen)
            % don't change mean_d
        else
            d = NaN(1, length(target_onscreen));
            for jj = 1:length(target_onscreen)
                if calculate_segment
                    d(jj) = point2segment([rvelocityX(onscreen(target_onscreen(jj)),ii); rvelocityY(onscreen(target_onscreen(jj)),ii)], v_constraint_far(:,onscreen(target_onscreen(jj)),ii), v_constraint_close(:,onscreen(target_onscreen(jj)),ii));
                else
                    dif = [rvelocityX(onscreen(target_onscreen(jj)),ii); rvelocityY(onscreen(target_onscreen(jj)),ii)] - v_constraint(:,onscreen(target_onscreen(jj)),ii);
                    d(jj) = norm(dif);
                end
            end
            mean_d(ii) = mean(d);

        end
        %
    end
    
    D(conditions(cond,1), conditions(cond,2)) = nanmean(mean_d);

end

