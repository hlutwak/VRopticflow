%% analysis
addpath(genpath('Data'))

%% plotting HMD data
% load oc.modelViewDataLeft
load('test-fixed-20230322T142931-0')
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

figure, plot3(pa.zPosition, pa.xPosition,ones(size(pa.xPosition))*pa.floorHeight+pa.paddleHalfHeight)
axis equal
xlabel('z')
ylabel('x')
zlabel('y')

% how much did you travel
vecnorm(oc.modelViewDataLeft(end-3:end-1,4)-oc.modelViewDataLeft(1:3,4))