% Wave-age color-axis definition matching the active binning scheme
%
% The color axis is mapped into log (or bin-index) space so that the n discrete
% colors span the n bins exactly. Tick labels stay in wave age
%
% Returns:
% * transform  - function handle mapping wave age onto the color axis
% * clims      - color axis limits, for clim()
% * ticks      - tick positions on the color axis, at the bin edges
% * ticklabels - tick labels, as wave-age values
% * axis_scale - scale for a wave-age plot axis ('log' or 'linear')
% * axis_ticks - round ticks for a wave-age plot axis
%
% Usage:
%   wc = wave_age_color();
%   scatter(x,y,sz,wc.transform(wave_age),'filled')
%   clim(wc.clims)
%   cbar.Ticks = wc.ticks; cbar.TickLabels = wc.ticklabels;
%
% Nathan Laxague 2026
%
function wc = wave_age_color()

s = load(repo_data_path('global_figure_settings.mat'));
edges = wave_age_bin_edges();

switch lower(s.wave_age_bin_spacing)

    case 'log'

        wc.transform = @log10;

    case 'quantile'

        % Equal-count edges are unequal in both linear and log space; bin-index
        % coordinates make each bin exactly one unit wide
        wc.transform = @(x) interp1(edges,0:numel(edges)-1,x,'linear','extrap');

    otherwise

        wc.transform = @(x) x;

end

wc.clims = wc.transform(edges([1 end]));
wc.ticks = wc.transform(edges);

wc.ticklabels = cell(1,length(edges));
for n = 1:length(edges)
    wc.ticklabels{n} = sprintf('%.0f',edges(n));
end

% Ticks for a wave-age plot axis (as opposed to a color axis): bin edges land on
% values like 14, 21, 43 and read as arbitrary, so round ticks are used instead.
% The axis stays logarithmic wherever the binning acknowledges the skew
if strcmpi(s.wave_age_bin_spacing,'linear')
    wc.axis_scale = 'linear';
    wc.axis_ticks = unique(round(edges));
else
    wc.axis_scale = 'log';
    wc.axis_ticks = [10 20 40 80];
end
