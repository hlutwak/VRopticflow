%% analysis
addpath(genpath('Data'))

%% plotting HMD data
% load oc.modelViewDataLeft
load('test-fixed-20230322T142931-0')
figure
plot3(oc.HMD(1:4:end,4), oc.HMD(2:4:end,4),...
    oc.HMD(3:4:end,4))