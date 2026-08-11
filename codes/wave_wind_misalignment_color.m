% Color-axis definition for the wave-wind misalignment, theta_m - theta_wind:
% the same quantity that forms the abscissa of every panel in
% wind_wave_subrange_directions, so a marker colored here reads directly against
% that figure.
%
% The misalignment wraps, and it fills the whole circle in these data, so it
% takes the cyclic twilight map over the full +/-180 deg rather than a diverging
% map over a narrow window: black on the wind, blue to the left of it, red to the
% right, white directly opposed.
%
% Returns:
% *value      - per-run misalignment [deg], theta_m - theta_wind
% *cmap       - cyclic colormap, for colormap()
% *clims      - color-axis limits, [-180 180]
% *ticks      - tick positions on the color axis
% *ticklabels - tick labels, in degrees
% *label      - colorbar title (LaTeX)
%
% Usage:
%   mc = wave_wind_misalignment_color();
%   scatter(x,y,sz,mc.value,'filled')
%   colormap(gca,mc.cmap); clim(mc.clims)
%
% N. Laxague 2026
%
function mc = wave_wind_misalignment_color()

m = subrange_wind_misalignment();

mc.value = m.dtheta_m(:);

mc.cmap = twilight(256);

mc.clims = [-180 180];
mc.ticks = -180:90:180;
mc.ticklabels = arrayfun(@(v) sprintf('%g',v),mc.ticks,'UniformOutput',false);
mc.label = '$\mathrm{\theta_m-\theta_{wind}\ [^\circ]}$';
