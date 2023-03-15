# 3D Pong for VR and AR

If this code is useful to you, and leads to some form of publication, please cite this paper (JM Fulvio, B Rokers. 2017. Use of Cues in Virtual Reality Depends on Visual Feedback. Nature Scientific Reports 7, Article number: 16009) and acknowledge us for sharing the code.

Also we ask that if you receive interest in the code from anyone else, please direct them to us. Lastly, we are happy for you to use or modify the code for all research, educational and private purposes, but we would like to note that it is not intended for commercial use.

To run in debug mode – without HMD connected, in SetupDisplay on line 10, set:

ds.oculusConnected = 0; % Is the HMD connected

For production mode, in RunExperiment, set DEBUG_FLAG = 0;


Best practices:

-If running the code with the Oculus Rift CV1, make sure to set the interpupillary-distance (IPD) slider on the headset to the appropriate setting and only trigger the code to start once the device is on the individual’s head. If triggered prior to then, syncing to the display will not be correct and the frame rate will be drastically-reduced from the 90 fps of the device.

-If running the code with the Oculus Rift DK2, make sure to set the IPD slider as well as the individual’s height in the Oculus Rift Configuration Utility.  View the demo scene and ensure that everything is properly centered with respect to the individual


Running through the experiment:

To begin the experiment, run RunExperiment.m
Press the spacebar once the opening message appears to begin the first trial.
To adjust the paddle, use the left and right arrow keys. To lock in the paddle setting, press the spacebar to initiate feedback. To advance to the next trial after the feedback, press the up arrow key.  Press Escape to quit out early.
