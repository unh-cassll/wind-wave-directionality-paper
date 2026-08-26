% Wave-age bin edges and centers, per the scheme set in the aa_step01 preamble
%
% Wave age scales roughly as 1/U10, so the sample is skewed and uniform bins
% leave the oldest bin nearly empty; log spacing spreads the counts, equal-count
% spacing balances them exactly
%
% Settings consumed from data/global_figure_settings.mat:
% *wave_age_lims        - [min max] of the binned range
% *wave_age_bin_spacing - 'log', 'linear' or 'quantile'
% *n_wave_age_bins      - bin count ('log' and 'quantile')
% *d_wave_age           - bin width ('linear' only)
%
% Returns:
% *edges   - 1 x (n+1) bin boundaries
% *centers - 1 x n bin centers: arithmetic midpoints if linear, geometric if
%             log, and the reference sample's per-bin median if quantile
%
% Nathan Laxague 2026
%
function [edges,centers] = wave_age_bin_edges()

s = load(repo_data_path('global_figure_settings.mat'));

lims = s.wave_age_lims;

switch lower(s.wave_age_bin_spacing)

    case 'log'

        if lims(1) <= 0
            error('wave_age_lims(1) must be positive for logarithmic bins.')
        end
        n = s.n_wave_age_bins;
        edges = lims(1)*(lims(2)/lims(1)).^((0:n)/n);
        centers = sqrt(edges(1:end-1).*edges(2:end));

    case 'linear'

        edges = lims(1):s.d_wave_age:lims(2);
        centers = edges(1:end-1) + diff(edges)/2;

    case 'quantile'

        % One canonical reference sample for every figure, outer edges pinned
        % to wave_age_lims so the bins tile the stated range
        n = s.n_wave_age_bins;
        wa = reference_wave_age();
        wa = wa(wa >= lims(1) & wa <= lims(2));

        if numel(wa) < 2*n
            error('Only %d reference runs inside wave_age_lims for %d equal-count bins.',numel(wa),n)
        end

        edges = prctile(wa,linspace(0,100,n+1));
        edges([1 end]) = lims;              % tile the stated range exactly

        if any(diff(edges) <= 0)
            error('Equal-count bins collapsed (tied wave ages); reduce n_wave_age_bins.')
        end

        centers = NaN(1,n);
        for i = 1:n
            if i < n
                in_bin = wa >= edges(i) & wa < edges(i+1);
            else
                in_bin = wa >= edges(i) & wa <= edges(i+1);
            end
            centers(i) = median(wa(in_bin));
        end

    otherwise

        error('wave_age_bin_spacing must be ''log'', ''linear'' or ''quantile''.')

end

% Canonical wave-age sample: gated, finite, on the 190-run production grid
function wa = reference_wave_age()

supporting_nc_name = repo_data_path('ASIT2019_supporting_environmental_observations.nc');

wind = wind_forcing(supporting_nc_name);
pk = load(repo_data_path('ASIT2019_peak_wave_phase_speed.mat'));

wa = pk.C_p(:)./wind.ustar;
wa = wa(isfinite(wa));
