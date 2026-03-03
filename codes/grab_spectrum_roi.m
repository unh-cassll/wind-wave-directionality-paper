
% Given handle of axis containing wavenumber/frequency-directional spectrum
% returns directional extrema of ROI
% Direction should be given in radians
%
% N. Laxague 2026
%
function roi_directional_extrema = grab_spectrum_roi(axis_handle,scale_vector)

h_roi = drawfreehand(axis_handle);

position_from_roi = h_roi.Position;

coordinate_direction = [position_from_roi(:,1); position_from_roi(:,1)];
coordinate_scale = [position_from_roi(:,2); position_from_roi(:,2)];

ind_max_left = find(coordinate_scale==max(coordinate_scale),1,'first');
ind_max_right = find(coordinate_scale==max(coordinate_scale),1,'last');

inds = ind_max_left:ind_max_right;

coordinate_scale = coordinate_scale(inds);
coordinate_direction = coordinate_direction(inds);

halflength = floor(length(coordinate_scale)/2);

scale_firsthalf = coordinate_scale(1:halflength);
direction_firsthalf = coordinate_direction(1:halflength);

ind_max_last = find(scale_firsthalf==max(scale_firsthalf),1,'last');
ind_min_first = find(scale_firsthalf==min(scale_firsthalf),1,'first');

scale_left = coordinate_scale(ind_max_last:ind_min_first);
direction_left = coordinate_direction(ind_max_last:ind_min_first);

coordinate_scale(1:ind_min_first) = [];
coordinate_direction(1:ind_min_first) = [];

ind_min_last = find(coordinate_scale==min(coordinate_scale),1,'last');
ind_max_first = find(coordinate_scale==max(coordinate_scale),1,'first');

scale_right = coordinate_scale(ind_min_last:ind_max_first);
direction_right = coordinate_direction(ind_min_last:ind_max_first);

[scale_left,inds_unique,~] = unique(scale_left);
direction_left = direction_left(inds_unique);

[scale_right,inds_unique,~] = unique(scale_right);
direction_right = direction_right(inds_unique);

direction_left_interp = interp1(scale_left,direction_left,scale_vector,'linear');
direction_right_interp = interp1(scale_right,direction_right,scale_vector,'linear');

if max(direction_right_interp) > max(direction_left_interp)
    roi_directional_extrema = [direction_left_interp(:) direction_right_interp(:)];
else
    roi_directional_extrema = [direction_right_interp(:) direction_left_interp(:)];
end
