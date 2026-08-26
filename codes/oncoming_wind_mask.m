% Oncoming-wind gate: 1 for runs whose wind approaches the instrument boom from
% ahead, NaN where it blows over the platform from behind and the measurements
% are flow-distorted
%
% Multiply per-run quantities by the returned mask to drop the gated runs. The
% heading comes from the aa_step01 preamble, as does the half-angle, which is
% per source: the boom-mounted sonic needs the tight ASIT_EC sector, the ASIT
% bulk record tolerates the wider ASIT_COARE one. The tighter of the two is the
% default, so an unqualified call can never be too permissive
%
% wind_forcing applies this itself, per record, using its own provenance; call
% it directly only for quantities that do not come through wind_forcing, such
% as the breaking-crest records
%
% You provide:
% *wind_dir_deg_coming_from - per-run wind direction the wind blows FROM [deg true]
% *gate_deg                 - half-angle [deg, default ASIT_EC_wind_gate_deg]
%
% Returns:
% *mask - 1 or NaN, same shape as the input
%
% Nathan Laxague 2026
%
function mask = oncoming_wind_mask(wind_dir_deg_coming_from,gate_deg)

s = load(repo_data_path('global_figure_settings.mat'));

if nargin < 2 || isempty(gate_deg); gate_deg = s.ASIT_EC_wind_gate_deg; end

diving_board_relative_angle = compute_relative_angle(wind_dir_deg_coming_from, ...
    wind_dir_deg_coming_from*0 + s.diving_board_heading_deg);

mask = 0*wind_dir_deg_coming_from + 1;
mask(abs(diving_board_relative_angle(:)) > gate_deg) = NaN;
