% Given two angles, computes angular difference

% accounts for phase wrapping at 0/360
% input in degrees
%
function angular_diff = compute_relative_angle(ang_A,ang_B)

ang_A = ang_A(:)';
ang_B = ang_B(:)';

vec_A = [sind(ang_A); cosd(ang_A); 0*ang_A];
vec_B = [sind(ang_B); cosd(ang_B); 0*ang_B];

A_cross_B = cross(vec_A,vec_B);

angular_diff = 180/pi*atan2(A_cross_B(3,:),dot(vec_A,vec_B));
