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

figure, scatter(pa.fullFactorial(:,1), pa.fullFactorial(:,2)), axis equal


%% screencap
s = size(screenCap);

figure
set(gcf,'position',[0, 0, 1280, 1600])

for ii = 1:s(3)
    imagesc(screenCap(:,:,ii)), colormap(gray)
    axis equal
    pause(1/10)
end


%% response data
% 90 is forward, 270 is backward
pcorrect = NaN(length(pa.speed),length(pa.direction));
for speed = 1:length(pa.speed)
    idx_speed = find(pa.fullFactorial(:,3) == pa.speed(speed));
    for direction = 1:length(pa.direction)
        idx_dir = find(pa.fullFactorial(:,4) == pa.direction(direction));
        idx = intersect(idx_speed, idx_dir);
        pcorrect(cond) = sum(eq(pa.LR(idx), pa.LRresponse(idx)))/pa.nRepeats;
    end
end
