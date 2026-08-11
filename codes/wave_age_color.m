% Wave-age color-axis definition matching the active binning scheme
%
% A linear color axis misrepresents non-uniform bins: the discrete colors spread
% evenly while the edges do not. Mapping the axis into log (or bin-index) space
% makes the n colors span the n bins exactly. Labels stay in wave age
%
% Returns:
% *transform  - function handle mapping wave age onto the color axis
% *clims      - color axis limits, for clim()
% *ticks      - tick positions on the color axis, at the bin edges
% *ticklabels - tick labels, as wave-age values
%
% Usage:
%   wc = wave_age_color();
%   scatter(x,y,sz,wc.transform(wave_age),'filled')
%   clim(wc.clims)
%   cbar.Ticks = wc.ticks; cbar.TickLabels = wc.ticklabels;
%
% N. Laxague 2026
%
function wc = wave_age_color()

s = load(repo_data_path('global_figure_settings.mat'));
edges = wave_age_bin_edges();

switch lower(s.wave_age_bin_spacing)
    case 'log'
        wc.transform = @log10;
    case 'quantile'
        % Equal-count edges are unequal in both linear and log space, so map
        % wave age onto bin-index coordinates: each bin is then exactly one
        % unit wide and the N discrete colors span the N bins
        wc.transform = @(x) interp1(edges,0:numel(edges)-1,x,'linear','extrap');
    otherwise
        wc.transform = @(x) x;
end

wc.clims = wc.transform(edges([1 end]));
wc.ticks = wc.transform(edges);
wc.ticklabels = arrayfun(@(v) sprintf('%.0f',v),edges,'UniformOutput',false);

% For a wave-age axis (as opposed to a color axis): the scale that matches the
% binning, with round integer ticks. Bin edges make poor axis ticks -- they land
% on values like 14, 21, 43 and read as arbitrary.
% The wave-age axis itself stays logarithmic whenever the binning acknowledges
% the skew, i.e. for both log and quantile spacing
if strcmpi(s.wave_age_bin_spacing,'linear')
    wc.axis_scale = 'linear';
    wc.axis_ticks = unique(round(edges));
else
    wc.axis_scale = 'log';
    wc.axis_ticks = nice_log_ticks(edges([1 end]));
end
