% Given wavenumber/frequency/inverse phase speed-directional spectrum and
% wind direction, arranges about wind direction +/- pi radians
%
% N. Laxague 2025
%
function [out_spect,out_dir_rad] = convert_dirspect_to_downwind(in_spect,in_dir_rad,wind_dir_going_to_rad)

in_dir_rad = in_dir_rad(:)';

[s1,~] = size(in_spect);

if length(in_dir_rad) == s1

    in_spect = in_spect';

end

super_spect = [in_spect in_spect in_spect in_spect in_spect];
super_dir = [in_dir_rad-4*pi in_dir_rad-2*pi in_dir_rad in_dir_rad+2*pi in_dir_rad+4*pi];

super_dir_rel_wind = super_dir - wind_dir_going_to_rad;

ind_left = find(abs(super_dir_rel_wind)<pi,1,'first');
ind_right = ind_left + length(in_dir_rad) - 1;


out_spect = super_spect(:,ind_left:ind_right);
out_dir_rad = super_dir_rel_wind(ind_left:ind_right);