% Wind speed and friction velocity from the source selected in the aa_step01
% preamble, so every figure draws its forcing from one place
%
% Sources (wind_source_choice in data/global_figure_settings.mat):
%   'EC_ASIT'      - eddy covariance at ASIT; U10 from the log law using the EC
%                    measurement height and z0 = z*exp(-kappa*U/u_*)
%   'COARE_ASIT'   - COARE 3.6 at ASIT: COARE_U10 and COARE_Ust
%   'EC_COARE_ASIT'- eddy covariance with ASIT COARE 3.6 backfilling the runs
%                    where the sonic has no record; w.provenance says which
%
% Everything here is measured at ASIT. The U10_best/ustar_best variables the
% supporting netCDF carries are deliberately not offered: they backfill from an
% off-platform buoy, which this analysis no longer draws on
%
% Two gates are applied here, so every figure inherits them through one code
% path and none has to reapply them:
%
%   oncoming wind - each source on its own half-angle, ASIT_EC_wind_gate_deg
%                   for the sonic and ASIT_COARE_wind_gate_deg for the bulk
%                   record, applied BEFORE the blended source decides where to
%                   backfill. A run the sonic loses to flow distortion is
%                   therefore a gap the bulk record can fill on its own wider
%                   sector, rather than a run discarded outright
%   ustar_cut     - minimum friction velocity; set it to 0 to disable
%
% Pass apply_cut = false to skip both, for figures documenting the conditions
% themselves rather than analyzing them
%
% You provide:
% *nc_name   - supporting netCDF to read from [default the 190-run production file]
% *apply_cut - apply the ustar_cut gate [default true]
%
% Returns:
% *U10    - 10 m wind speed [m/s]
% *ustar  - air-side friction velocity [m/s]
% *tau    - wind stress rho_a*ustar^2 [N/m^2]
% *source - the source actually used
% *label  - human-readable description, for figure annotation
% *n_cut  - runs removed by the ustar gate
% *provenance - per record, which source supplied the value: 1 primary, 2 fill,
%               NaN where nothing did (blended sources only, else all ones)
%
% Nathan Laxague 2026
%
function w = wind_forcing(nc_name,apply_cut)

if nargin < 1 || isempty(nc_name)
    nc_name = repo_data_path('ASIT2019_supporting_environmental_observations.nc');
end
if nargin < 2 || isempty(apply_cut); apply_cut = true; end

s = load(repo_data_path('global_figure_settings.mat'));
choice = s.wind_source_choice;

vars = {ncinfo(nc_name).Variables.Name};
has = @(v) any(strcmp(vars,v));

rho_a = 1.25;
kappa = 0.4;

% Oncoming-wind masks, built up front so the blended source can gate the sonic
% BEFORE deciding where to backfill: a run the sonic loses to flow distortion
% is a genuine gap, and the bulk record is entitled to fill it on its own wider
% sector. Gating after the merge would instead discard that run outright
if has('COARE_Wdir')
    wind_dir_deg_coming_from = ncread(nc_name,'COARE_Wdir');
else
    wind_dir_deg_coming_from = ncread(nc_name,'COARE_winddir');
end

if apply_cut
    mask_EC = oncoming_wind_mask(wind_dir_deg_coming_from,s.ASIT_EC_wind_gate_deg);
    mask_COARE = oncoming_wind_mask(wind_dir_deg_coming_from,s.ASIT_COARE_wind_gate_deg);
else
    mask_EC = ones(size(wind_dir_deg_coming_from));
    mask_COARE = mask_EC;
end

mask_EC = mask_EC(:);
mask_COARE = mask_COARE(:);

switch lower(choice)

    case 'ec_asit'

        [U10,ustar] = read_ec_asit(nc_name,has,kappa);
        available = isfinite(ustar(:));
        U10 = mask_EC.*U10(:);
        ustar = mask_EC.*ustar(:);
        w.source = 'EC_ASIT';
        w.label = 'eddy covariance, ASIT';

    case 'coare_asit'

        [U10,ustar] = read_coare_asit(nc_name,has);
        available = isfinite(ustar(:));
        U10 = mask_COARE.*U10(:);
        ustar = mask_COARE.*ustar(:);
        w.source = 'COARE_ASIT';
        w.label = 'COARE 3.6, ASIT';

    case 'ec_coare_asit'

        [U10,ustar] = read_ec_asit(nc_name,has,kappa);
        [U10_fill,ustar_fill] = read_coare_asit(nc_name,has);

        available = isfinite(ustar(:)) | isfinite(ustar_fill(:));

        % Gate each source on its own sector first
        U10 = mask_EC.*U10(:);
        ustar = mask_EC.*ustar(:);
        U10_fill = mask_COARE.*U10_fill(:);
        ustar_fill = mask_COARE.*ustar_fill(:);

        % Then backfill whatever the sonic no longer supplies, whether it never
        % had the run or lost it to the gate. The ustar_cut is deliberately not
        % part of this: a sonic reporting weak wind is a measurement, not a gap,
        % so it keeps its value and meets the cut on its own further down
        fill = ~isfinite(ustar) & isfinite(ustar_fill);
        U10(fill) = U10_fill(fill);
        ustar(fill) = ustar_fill(fill);

        provenance = NaN(size(ustar));
        provenance(isfinite(ustar)) = 1;
        provenance(fill) = 2;

        w.source = 'EC_COARE_ASIT';
        w.label = 'eddy covariance, ASIT with ASIT COARE backfill';

    otherwise

        error(['wind_source_choice must be ''EC_ASIT'', ''COARE_ASIT'' ' ...
            'or ''EC_COARE_ASIT''.'])

end

w.n_gated = sum(available & ~isfinite(ustar(:)));

if ~exist('provenance','var')
    provenance = ones(size(U10));
    provenance(~isfinite(U10)) = NaN;
end

w.U10 = U10(:);
w.ustar = ustar(:);
w.provenance = provenance(:);

% Minimum-friction-velocity gate; below u_* ~ 0.1 m/s the short waves decouple
% from the local wind and the wind direction itself is poorly defined
w.n_cut = 0;
if apply_cut && isfield(s,'ustar_cut') && s.ustar_cut > 0
    below = w.ustar < s.ustar_cut;
    w.n_cut = sum(below);
    w.U10(below) = NaN;
    w.ustar(below) = NaN;
    w.provenance(below) = NaN;
    w.label = sprintf('%s, u_* > %g m/s',w.label,s.ustar_cut);
end

w.tau = rho_a*w.ustar.^2;


% Eddy covariance at ASIT; U10 from the log law using the measurement height
function [U10,ustar] = read_ec_asit(nc_name,has,kappa)

if has('EC_ustar_m_s')
    ustar = ncread(nc_name,'EC_ustar_m_s');
    U = ncread(nc_name,'EC_U_m_s');
    z = ncread(nc_name,'EC_z_m_above_water');
else
    ustar = ncread(nc_name,'EC_ustar');
    U = ncread(nc_name,'EC_U');
    z = ncread(nc_name,'EC_z');
end

z0 = z.*exp(-kappa*U./ustar);
U10 = ustar/kappa.*log(10./z0);


% COARE 3.6 driven by ASIT observations
function [U10,ustar] = read_coare_asit(nc_name,has)

U10 = ncread(nc_name,'COARE_U10');
if has('COARE_Ust')
    ustar = ncread(nc_name,'COARE_Ust');
else
    ustar = ncread(nc_name,'COARE_ustar');
end
