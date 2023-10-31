function  device = SetupEyetracker()
% running python through matlab

% make sure eyetracker and computer connected to the same network
% following code comes from:
% https://docs.pupil-labs.com/neon/real-time-api/introduction/
% look at % https://pupil-labs-realtime-api.readthedocs.io/en/stable/examples/index.html
% trouble shooting https://docs.pupil-labs.com/neon/troubleshooting/#i-can-not-connect-to-devices-using-the-real-time-api


% make sure it's connected
py.nest_asyncio.apply()

import py.pupil_labs.realtime_api.simple.Device
ip = "10.17.225.131";
device = Device(ip, "8080");
% companion device name, are the glasses connected to companion
% disp(device.phone_name)
% disp(device.module_serial)

