% Mean propagation direction of each spectral subrange, and its misalignment
% with the wind, per run
%
% Subrange directions are energy-weighted circular means of the directional
% spectrum over each band:
% * equilibrium and saturation, from the combined frequency/wavenumber spectra
% * short gravity, between k_sat_end and 117 rad/m
% * gravity-capillary, over 117-371 rad/m and within +/-120 deg of the short
%   gravity direction, matching the window used for the 50% spreading limits
%
% Runs with the wind blowing over the platform from behind are dropped
%
% Returns (all per run, in degrees):
% * theta_m/eq/sat/sg/gc  - subrange mean directions [deg true, going to]
% * wind_dir_deg_going_to - wind direction
% * dtheta_*              - subrange direction minus wind direction
% * wave_age              - c_p/u_*
% * breakers, current     - matching pairs for the two non-spectral panels
%
% Nathan Laxague 2026
%
function m = subrange_wind_misalignment()

load('data/frequency_spect_range_limits.mat','f_eq_start','f_eq_end');
k_limits = load('data/wavenumber_spect_range_limits.mat');
k_sat_end = k_limits.k_sat_end;

short_wave_nc_name = 'data/ASIT2019_wave_spectra_stats_timeseries_empirical_gain.nc';
long_wave_nc_name = 'data/ASIT2019_EPSS_directional_spectra.nc';
supporting_nc_name = 'data/ASIT2019_supporting_environmental_observations.nc';

integrated_wave_breaking_quantities = load('data/integrated_wave_breaking_quantities.mat');

wind_driven_current_data = load('data/ASIT2019_skin_Doppler_subsurface_velocities.mat');

current_Winddir = wind_driven_current_data.depth_01cm_Winddir;
current_U_dir = wind_driven_current_data.depth_01cm_Udir;
current_Wavedir = wind_driven_current_data.depth_01cm_Wavedir;

f_Hz_Pyxis = double(ncread(short_wave_nc_name,'f_Hz'));
k_rad_m = double(ncread(short_wave_nc_name,'k_rad_m'));
theta_rad_Pyxis = double(ncread(short_wave_nc_name,'theta_rad'));

wind_dir_deg_coming_from = ncread(supporting_nc_name,'COARE_Wdir');
wind_dir_deg_going_to = mod(wind_dir_deg_coming_from+180,360);
% Wind forcing per the source selected in the aa_step01 preamble
wind = wind_forcing(supporting_nc_name);
EC_ustar_m_s = wind.ustar;

SFTHETA = double(ncread(short_wave_nc_name,'S_f_theta'));
SKTHETA = double(ncread(short_wave_nc_name,'S_k_theta'));

f_Hz_EPSS = double(ncread(long_wave_nc_name,'frequency'));

% Auto-detect direction units (legacy: deg/per-deg, current: rad/per-rad)
dir_EPSS = double(ncread(long_wave_nc_name,'direction'));
FFTHETA = double(ncread(long_wave_nc_name,'F_f_d'));
if max(abs(dir_EPSS)) > 2*pi
    theta_rad_EPSS = dir_EPSS*pi/180;
    FFTHETA = FFTHETA*180/pi;
else
    theta_rad_EPSS = dir_EPSS;
end

f_Hz_EPSS = f_Hz_EPSS(2:end);
FFTHETA = FFTHETA(:,:,2:end);

FFTHETA = permute(FFTHETA,[3 2 1]);

SFTHETA(SFTHETA==0) = NaN;
FFTHETA(FFTHETA==0) = NaN;

f_cut_low = 0.75;
f_cut_high = 1.25;

FFTHETA = FFTHETA(f_Hz_EPSS<f_cut_low,:,:);
f_Hz_EPSS = f_Hz_EPSS(f_Hz_EPSS<f_cut_low);

SFTHETA = SFTHETA(f_Hz_Pyxis>f_cut_high,:,:);
f_Hz_Pyxis = f_Hz_Pyxis(f_Hz_Pyxis>f_cut_high);

water_depth_m = 15;
[c_disp,~] = lindisp_with_current(2*pi*f_Hz_Pyxis,water_depth_m,0);
k_disp_Pyxis = 2*pi*f_Hz_Pyxis./c_disp;

FFTHETA_direct_big = k_disp_Pyxis.^-2.*[SFTHETA SFTHETA SFTHETA];

theta_rad_direct_big = [theta_rad_Pyxis(:)'-2*pi theta_rad_Pyxis(:)' theta_rad_Pyxis(:)'+2*pi];
inds_keep_direct = theta_rad_direct_big >= -pi & theta_rad_direct_big < pi;

% Shared Pyxis direction grid ([-pi,pi)) for the combined spectrum
theta_rad = theta_rad_direct_big(inds_keep_direct);

% Resample E-PSS spectrum from its native direction grid onto the Pyxis grid
FFTHETA_EPSS_resamp = NaN*ones(size(FFTHETA,1),length(theta_rad),size(FFTHETA,3));
for run_ind = 1:size(FFTHETA,3)
    FFTHETA_EPSS_resamp(:,:,run_ind) = regrid_directional_spectrum(squeeze(FFTHETA(:,:,run_ind)),theta_rad_EPSS,theta_rad);
end

FFTHETA_combined = [FFTHETA_EPSS_resamp; FFTHETA_direct_big(:,inds_keep_direct,:)];

f_Hz_combined = [f_Hz_EPSS; f_Hz_Pyxis];

frequency_spect = load('data/ASIT2019_combined_frequency_slope_spectra.mat');

f_p = sum(frequency_spect.f_Hz_combined.*frequency_spect.F_f_block.^4,1,'omitnan')./sum(frequency_spect.F_f_block.^4,1,'omitnan');
f_p = f_p(:);

[c_p,~] = lindisp_with_current(2*pi*f_p,water_depth_m,0);

wave_age = c_p./EC_ustar_m_s;

Sm = squeeze(mean(sin(theta_rad(:)').*FFTHETA_combined,[1 2],'omitnan'));
Cm = squeeze(mean(cos(theta_rad(:)').*FFTHETA_combined,[1 2],'omitnan'));

Vm = atan2(Sm,Cm);

theta_m = 180/pi*Vm;

theta_eq = NaN*theta_m;
theta_sat = theta_eq;
theta_sg = theta_eq;
theta_gc = theta_eq;

inds_gc = k_rad_m > 117 & k_rad_m < 371;

for n = 1:length(f_eq_start)

    Fftheta_slice = squeeze(FFTHETA_combined(:,:,n));

    inds_eq = f_Hz_combined >= f_eq_start(n) & f_Hz_combined <= f_eq_end(n);

    Sm = mean(sin(theta_rad(:)').*Fftheta_slice(inds_eq,:),'all','omitnan');
    Cm = mean(cos(theta_rad(:)').*Fftheta_slice(inds_eq,:),'all','omitnan');
    Vm = atan2(Sm,Cm);
    theta_eq(n) = 180/pi*Vm;

    Sktheta_slice = squeeze(SKTHETA(:,:,n));

    inds_sat = k_rad_m <= k_sat_end(n);

    Sm = mean(sin(theta_rad_Pyxis(:)').*k_rad_m(inds_sat).^-2.*Sktheta_slice(inds_sat,:),'all','omitnan');
    Cm = mean(cos(theta_rad_Pyxis(:)').*k_rad_m(inds_sat).^-2.*Sktheta_slice(inds_sat,:),'all','omitnan');
    Vm = atan2(Sm,Cm);
    theta_sat(n) = 180/pi*Vm;

    inds_sg = k_rad_m > k_sat_end(n) & k_rad_m < 117;

    Sm = mean(sin(theta_rad_Pyxis(:)').*k_rad_m(inds_sg).^-2.*Sktheta_slice(inds_sg,:),'all','omitnan');
    Cm = mean(cos(theta_rad_Pyxis(:)').*k_rad_m(inds_sg).^-2.*Sktheta_slice(inds_sg,:),'all','omitnan');
    Vm = atan2(Sm,Cm);
    theta_sg(n) = 180/pi*Vm;

    Sktheta_big = [Sktheta_slice Sktheta_slice Sktheta_slice];
    theta_rad_big = [theta_rad_Pyxis(:)'-2*pi theta_rad_Pyxis(:)' theta_rad_Pyxis(:)'+2*pi] - Vm;
    mask = NaN*Sktheta_big;
    mask(:,theta_rad_big>=-2*pi/3&theta_rad_big<2*pi/3) = 1;

    Sm = mean(mask(inds_gc,:).*sin(theta_rad_big(:)').*k_rad_m(inds_gc).^-2.*Sktheta_big(inds_gc,:),'all','omitnan');
    Cm = mean(mask(inds_gc,:).*cos(theta_rad_big(:)').*k_rad_m(inds_gc).^-2.*Sktheta_big(inds_gc,:),'all','omitnan');
    Vm = atan2(Sm,Cm);
    theta_gc(n) = 180/pi*Vm + theta_sg(n);

end

theta_m = mod(theta_m+360,360);
theta_eq = mod(theta_eq+360,360);
theta_sat = mod(theta_sat+360,360);
theta_sg = mod(theta_sg+360,360);
theta_gc = mod(theta_gc+360,360);

% The run-grid oncoming-wind gate is applied inside wind_forcing; carry it onto
% the directions so gated runs drop out here too
wind_dir_deg_going_to(~isfinite(EC_ustar_m_s(:))) = NaN;

% Breaking-crest records store the direction the wind blows towards, and carry
% bulk COARE wind rather than the sonic, so they take the COARE half-angle
gset = load(repo_data_path('global_figure_settings.mat'));
integrated_wave_breaking_quantities.wdir_deg = ...
    oncoming_wind_mask(mod(integrated_wave_breaking_quantities.wdir_deg+180,360), ...
    gset.ASIT_COARE_wind_gate_deg).*integrated_wave_breaking_quantities.wdir_deg;

% Offshore-propagating crest gate (+/-90 degrees), currently disabled
% integrated_wave_breaking_quantities.theta_br(integrated_wave_breaking_quantities.theta_br<270 & integrated_wave_breaking_quantities.theta_br > 90) = NaN;

% Contiguous acquisition blocks, for tests that must tolerate the hour-to-hour
% correlation between runs. Breaking crests sit on their own acquisition list
% (camera 02, 212 records), not the 190-run production grid
t_run = ncread(supporting_nc_name,'time');
L = load('data/ASIT2019_Lambda_c_theta_all.mat','DTime_02');
t_br = 86400*datenum(L.DTime_02);

m.block_id = cumsum([1; diff(t_run)/3600 > 3]);

m.theta_m = theta_m;
m.theta_eq = theta_eq;
m.theta_sat = theta_sat;
m.theta_sg = theta_sg;
m.theta_gc = theta_gc;
m.wind_dir_deg_going_to = wind_dir_deg_going_to;
m.wave_age = wave_age;

m.dtheta_m = compute_relative_angle(theta_m,wind_dir_deg_going_to)';
m.dtheta_eq = compute_relative_angle(theta_eq,wind_dir_deg_going_to)';
m.dtheta_sat = compute_relative_angle(theta_sat,wind_dir_deg_going_to)';
m.dtheta_sg = compute_relative_angle(theta_sg,wind_dir_deg_going_to)';
m.dtheta_gc = compute_relative_angle(theta_gc,wind_dir_deg_going_to)';

% Breaking crests and wind-driven current: reference/particular pairs
m.breakers.D_ref = compute_relative_angle(integrated_wave_breaking_quantities.D_m_deg,integrated_wave_breaking_quantities.wdir_deg)';
m.breakers.D_particular = compute_relative_angle(integrated_wave_breaking_quantities.theta_br,integrated_wave_breaking_quantities.wdir_deg)';
m.breakers.wave_age = integrated_wave_breaking_quantities.c_p_m_s./integrated_wave_breaking_quantities.ustar_m_s;
m.breakers.block_id = cumsum([1; diff(t_br)/3600 > 3]);

m.current.D_ref = compute_relative_angle(current_Wavedir,current_Winddir)';
m.current.D_particular = compute_relative_angle(current_U_dir,current_Winddir)';
m.current.wave_age = wave_age;
