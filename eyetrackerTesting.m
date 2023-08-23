%% run python through matlab

% addpath(genpath('C:\Users\hlutw\OneDrive\Documents\MATLAB'))

% make sure eyetracker and computer connected to the same network
% following code comes from:
% https://docs.pupil-labs.com/neon/real-time-api/introduction/
% look at % https://pupil-labs-realtime-api.readthedocs.io/en/stable/examples/index.html
% trouble shooting https://docs.pupil-labs.com/neon/troubleshooting/#i-can-not-connect-to-devices-using-the-real-time-api


% make sure it's connected
py.nest_asyncio.apply()

import py.pupil_labs.realtime_api.simple.Device
ip = "10.17.201.214";
device = Device(ip, "8080");
% companion device name, are the glasses connected to companion
disp(device.phone_name)
disp(device.module_serial)

%% start and stop a recording
recording_duration = 2; %recording time in seconds

recording_id = device.recording_start();
disp(['Started recording with id', string(recording_id)]);

py.time.sleep(recording_duration)

device.recording_stop_and_save()

%%
recording_id = device.recording_start();
disp(['Started recording with id', string(recording_id)]);

py.time.sleep(2)

py.print(device.send_event("test event 1"))

py.time.sleep(2)

% send event with current timestamp
py.print(device.send_event("test event 2", py.time.time_ns()))

py.time.sleep(2)

device.recording_stop_and_save()

%%
