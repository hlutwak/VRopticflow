function  d = point2segment(point, starting, ending)
% from http://www.fundza.com/vectors/point2line/index.html

% # Given a line segment with coordinates 'starting' (far distance on constraint)  and 'ending' (near distance on constraint) and the
% # coordinates of a vel 'point' the function returns the shortest 
% # distance from point to the segment and the coordinates of the 
% # nearest point on the line.

% 1 change coordinates so start of segment is at the origin
segment_vec = starting- ending;
point_vec = point - ending;

% 2 get projection of point_vec onto segment_vec
proj = dot(point_vec, segment_vec)/norm(segment_vec);

% 3 if it's negative, closest point is start (far)
% if it's larger than length of segment, closest point is end (near)
% if it's less than or equal to segment, distance is the perpendicular 

if proj < 0
    d = norm(point_vec);
elseif proj > norm(segment_vec)
    d = norm(point_vec - segment_vec);
else 
    theta = acos(dot(point_vec, segment_vec)/(norm(point_vec)*norm(segment_vec)));
    d = norm(point_vec)*sin(theta);
end


%
% test point 2 segment
% points = -5+10*rand(2, 50);
% starting = [0,1]';
% ending = [2, -3]';
% 
% d = zeros(1, length(points));
% for ii = 1:length(points)
%     d(ii) = point2segment(points(:,ii), starting, ending);
% end
% 
% figure, scatter(points(1,:), points(2,:))
% hold on, plot([starting(1), ending(1)], [starting(2), ending(2)], 'k')
