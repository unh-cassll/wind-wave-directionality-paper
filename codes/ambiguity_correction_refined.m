% Sources the ambiguity-correction constants used by
% compute_theta_halfwidth_and_delta:
% * R0_upwind, the genuine upwind energy fraction, taken as the median
%   mirror-lobe ratio over the directionally resolved band (k = 20-60 rad/m)
% * the ambiguity onset wavenumber, where the mirror lobe first reaches 90% of
%   the true lobe
%
% Also prints Delta uncorrected, with a flat factor of two, and with the
% continuous leakage factor, at a set of representative wavenumbers
%
% Diagnostic only: writes nothing. Run from codes/
%
% Nathan Laxague 2026
%
close all;clear;clc

g = 9.806;
sigma = 0.072;
rho_w = 1020;
water_depth_m = 15;

% Pyxis acquisition rate, for the samples-per-wave-period check
f_sample_Hz = 30;
f_nyquist_Hz = f_sample_Hz/2;

short_wave_nc_name = '../data/ASIT2019_wave_spectra_stats_timeseries_empirical_gain.nc';

S_k_theta_block = double(ncread(short_wave_nc_name,'S_k_theta'));
k_rad_m = double(ncread(short_wave_nc_name,'k_rad_m'));
theta_rad = double(ncread(short_wave_nc_name,'theta_rad'));

% Order the block as (wavenumber, direction, run) whichever way it was stored
sz = size(S_k_theta_block);
dim_k = find(sz==numel(k_rad_m),1);
dim_theta = find(sz==numel(theta_rad),1);
dim_run = setdiff(1:3,[dim_k dim_theta]);

S_k_theta_block = permute(S_k_theta_block,[dim_k dim_theta dim_run]);

n_k = size(S_k_theta_block,1);
n_theta = size(S_k_theta_block,2);
nrun = size(S_k_theta_block,3);

n_theta_half = n_theta/2;

% Mirror-lobe ratio: energy 180 deg from the peak, over the peak itself
R_mirror = NaN*ones(n_k,nrun);

for run_num = 1:nrun

    S_slice = S_k_theta_block(:,:,run_num);

    for i = 1:n_k

        S_row = S_slice(i,:);

        if all(~isfinite(S_row)) || max(S_row) <= 0
            continue
        end

        [S_peak,ind_peak] = max(S_row);
        ind_mirror = mod(ind_peak-1+n_theta_half,n_theta)+1;

        if isfinite(S_row(ind_mirror))
            R_mirror(i,run_num) = S_row(ind_mirror)/S_peak;
        end

    end

end

R_mirror_med = median(R_mirror,2,'omitnan');

% Frequency of each wavenumber, full gravity-capillary dispersion
omega = sqrt((g*k_rad_m+sigma/rho_w*k_rad_m.^3).*tanh(k_rad_m*water_depth_m));
f_of_k = omega/(2*pi);

% Resolved-band baseline: the smallest sustained ratio, taken as genuine
% upwind energy rather than ambiguity leakage
inds_band = k_rad_m > 20 & k_rad_m < 60 & isfinite(R_mirror_med);
R0 = median(R_mirror_med(inds_band));

fprintf('resolved-band baseline R0 = %.3f (k = 20-60 rad/m)\n',R0);

% Leakage runs 0 where the lobes are resolved to 1 where they are not, so the
% downwind correction factor runs from 1 to 2
leakage = max(R_mirror_med-R0,0)./max(1-R0,eps);
F_leak = 1 + leakage;

delta_data = load('../data/delta_with_wind_speed_data.mat');
delta_k = delta_data.delta_k;

delta_k_mean = mean(delta_k,2,'omitnan');
n_k_delta = min(numel(delta_k_mean),n_k);

% Delta implies a crosswind/downwind ratio, which the correction factor then
% rescales
S_ratio_90_0 = (1-delta_k_mean(1:n_k_delta))./(1+delta_k_mean(1:n_k_delta));
F_leak_k = F_leak(1:n_k_delta);

delta_k_flat = (2-S_ratio_90_0)./(2+S_ratio_90_0);
delta_k_continuous = (F_leak_k-S_ratio_90_0)./(F_leak_k+S_ratio_90_0);

fprintf('\n     k       f      R     factor   Delta_raw  Delta_x2  Delta_cont\n');

k_query = [5 10 20 50 80 100 130 160 181 200 250 300 371];

for n = 1:numel(k_query)

    [~,ind_k] = min(abs(k_rad_m(1:n_k_delta)-k_query(n)));

    if isfinite(delta_k_mean(ind_k))
        fprintf('%6.1f  %6.2f  %5.3f  %6.3f  %9.3f %9.3f %11.3f\n', ...
            k_rad_m(ind_k),f_of_k(ind_k),R_mirror_med(ind_k),F_leak_k(ind_k), ...
            delta_k_mean(ind_k),delta_k_flat(ind_k),delta_k_continuous(ind_k));
    end

end

% Ambiguity onset: the first wavenumber at which the mirror lobe reaches 90%
% of the true lobe
inds_finite = find(isfinite(R_mirror_med));
ind_onset = inds_finite(find(R_mirror_med(inds_finite)>0.90,1,'first'));

fprintf('\nambiguity onset (R>0.90): k = %.1f rad/m, f = %.2f Hz = f_nyq/%.1f\n', ...
    k_rad_m(ind_onset),f_of_k(ind_onset),f_nyquist_Hz/f_of_k(ind_onset));
fprintf('samples per wave period at that frequency: %.1f\n',f_sample_Hz/f_of_k(ind_onset));
