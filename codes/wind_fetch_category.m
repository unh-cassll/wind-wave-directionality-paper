% Discrete wind-fetch categories for the 190-run grid, shared by the figures
% that split runs by fetch. Fetch comes from the gated run-grid wind direction
% through asit_fetch_vs_azimuth
%
% Returns:
% * cat      - per-run category: 1 offshore, 2 cross-shore, 3 onshore, 0 no wind
% * colors   - category colors from the figure_style sweep (crimson, teal, violet)
% * labels   - display names, for legends
% * fetch_km - per-run fetch [km]
%
% Nathan Laxague 2026
%
function fc = wind_fetch_category()

% Category boundaries [km]
fetch_offshore_km = 10;
fetch_onshore_km = 100;

m = subrange_wind_misalignment();

fetch_km = asit_fetch_vs_azimuth(mod(m.wind_dir_deg_going_to+180,360));
fc.fetch_km = fetch_km(:);

fc.cat = 0*fc.fetch_km;
fc.cat(fc.fetch_km < fetch_offshore_km) = 1;
fc.cat(fc.fetch_km >= fetch_offshore_km & fc.fetch_km < fetch_onshore_km) = 2;
fc.cat(fc.fetch_km >= fetch_onshore_km) = 3;
fc.cat(~isfinite(fc.fetch_km)) = 0;

color_order = get(groot,'DefaultAxesColorOrder');
fc.colors = color_order([4 2 1],:);
fc.labels = {'offshore','cross-shore','onshore'};
