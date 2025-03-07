function [velocity_field] = constraint_velocity_screen(Z, xy, translate, view_dist, omegax)
% The world is a cloud of enveloped pattern motion, the observer is translating and and fixating past it.
% The observer is always at the origin, and is looking down the Z axis.
% This program takes a set of gabors and calculates the velocity of each
% one based on given translational components (for phase translating).
%
% XYZ is the XYZ position (in m) of the locations where we are measuring
% the speed.
%
% xy is the x,y positions of the dots on the screen (in m).
%
% translate tells you the observer's translation speed (Vx, Vy, Vz) in
% m/sec.
%
% view_dist is how far the observer is the from the screen in meters.
%
% velocity_field is a matrix (2x[number of locations]) reporting
% the x and y velocity in m/s on the screen.
%
% HL
% Updated 2/25/2020
% Updated 9/2/2020 takes in just depth for first argument
% Updated 2/22/2022 calculate velocity on imaging plane

% if nargin==0
%     [XYZ, centers_deg] = make_dot_plane(.1, 12.5, [40 32], [-20 -3 20 3]);
%     translate = [0.14 0 1.9];                                               % the observer is walking forward at 1.9 m/s, a brisk walk
%     view_dist = .57;
% end 

% initialize empty matrix to hold velocities
velocity_field = zeros(size(xy));                                 % number of patches, X and Y vel

% centers_m = XYZ(3,:).*tand(centers_deg);

%solve for the inverse depth at the X,Y,Z position that projects to the
%desired x,y position, then use that to compute the flow vector at each
%location (same Z value for all of the gabors).
% 
% for ii=1:length(Z)
%         inv_depth = 1/(Z(ii));
%         
%         % Heeger + Jepson, Lutwak + Bonnen + Simoncelli
%         x = xy(1,ii); % coordinates of screen
%         y = xy(2,ii);
%         f = view_dist;
%         Pxy = inv_depth;
%         
%         A = [-f 0 x; 0 -f y];
%         B = [f+x^2/f x*y/f 0; x*y/f f+y^2/f 0];
%         velocity_field(1:2,ii) = (z0*Pxy*A+B)*translate'/z0;
% end

%velocities in m/s

% % adding in eye rotations, instead of z0 put theta (omegax)
for ii=1:length(Z)
        inv_depth = 1/(Z(ii));
        
        % Heeger + Jepson
        x = xy(1,ii); % coordinates of screen
        y = xy(2,ii);
        f = view_dist;
        Pxy = inv_depth;
        
        A = [-f 0 x; 0 -f y];
        B = [(x*y)/f -(f+x^2/f) y; f+y^2/f -(x*y)/f -x];
%         rotate = [1/z0*translate(2), 0, 0];
        rotate = [omegax, 0, 0];
        velocity_field(1:2,ii) = Pxy*A*translate'+B*rotate';
end
