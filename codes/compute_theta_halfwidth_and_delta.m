% Regenerates, per run:
% * the 50% directional spreading limits (theta_halfwidth edges, via
%   fifty_percent_edges on the downwave-frame spreading function)
% * the spreading ratio Delta, with wavenumber/frequency-dependent centering
% * the wavenumber at which the upwind and downwind lobes stop being separable
%
% Writes ../data/directional_spreading_quantification_data.mat and
% ../data/delta_with_wind_speed_data.mat
%
% Nathan Laxague 2026
%
close all;clear;clc

load('../data/global_figure_settings.mat')

% Skip the expensive recompute unless explicitly requested
if ~exist('recompute_derived_products','var') || ~recompute_derived_products
    return
end

smooth_type = 'movmean';

% Genuine upwind energy fraction: the median mirror-lobe ratio over the
% directionally resolved band, k = 20-60 rad/m. Separates real upwind energy
% from 180-degree ambiguity leakage. Derived by ambiguity_correction_refined
R0_upwind = 0.165;

f_smooth = 5;
k_smooth = 5;
ang_smooth = 3;

short_wave_nc_name = '../data/ASIT2019_wave_spectra_stats_timeseries_empirical_gain.nc';

S_k_theta_block = double(ncread(short_wave_nc_name,'S_k_theta'));
S_f_theta_block = double(ncread(short_wave_nc_name,'S_f_theta'));

f_Hz_Pyxis = ncread(short_wave_nc_name,'f_Hz');
k_rad_m_Pyxis = ncread(short_wave_nc_name,'k_rad_m');
theta_rad_Pyxis = ncread(short_wave_nc_name,'theta_rad');

nrun = size(S_f_theta_block,3);
ntheta = length(theta_rad_Pyxis);

big_S_k_theta = [S_k_theta_block S_k_theta_block S_k_theta_block];
big_S_f_theta = [S_f_theta_block S_f_theta_block S_f_theta_block];

for i = 1:nrun
    big_S_k_theta(:,:,i) = smoothdata2(squeeze(big_S_k_theta(:,:,i)),smooth_type,{k_smooth,ang_smooth},'omitnan');
    big_S_f_theta(:,:,i) = smoothdata2(squeeze(big_S_f_theta(:,:,i)),smooth_type,{f_smooth,ang_smooth},'omitnan');
end

f_low = 0.6;

big_S_k_theta(k_rad_m_Pyxis<k_low,:,:) = NaN;
big_S_f_theta(f_Hz_Pyxis<f_low,:,:) = NaN;

% Per-wavenumber min-max normalization: the offset field, used for the 50%
% spreading limits and for locating the run-level peak direction
S_k_theta_star = (big_S_k_theta-min(big_S_k_theta,[],2,'omitnan'))./(max(big_S_k_theta,[],2,'omitnan')-min(big_S_k_theta,[],2,'omitnan'));
S_f_theta_star = (big_S_f_theta-min(big_S_f_theta,[],2,'omitnan'))./(max(big_S_f_theta,[],2,'omitnan')-min(big_S_f_theta,[],2,'omitnan'));

S_theta = squeeze(sum(S_f_theta_star,1,'omitnan'));

bigtheta = [theta_rad_Pyxis-2*pi; theta_rad_Pyxis; theta_rad_Pyxis+2*pi];

D_spread_holder_struc = struct();

delta_f = NaN*ones(length(f_Hz_Pyxis),nrun);
delta_k = NaN*ones(length(k_rad_m_Pyxis),nrun);

% Per-run wavenumber beyond which the upwind and downwind lobes are no longer
% separable. Currents and bound-wave content shift this scale run to run, so it
% is measured rather than assumed
k_ambiguity_run = NaN*ones(nrun,1);

k_ind_high = find(k_rad_m_Pyxis<371,1,'last');
mid_cols = ntheta + (1:ntheta);

for example_run_ind = 1:nrun

    % Energy-weighted mean wave direction about the frequency-integrated peak
    ind_p = find(S_theta(:,example_run_ind)==max(S_theta(:,example_run_ind),[],'all','omitnan'),1,'first');
    if ind_p < 37
        ind_p = ind_p + 72;
    end
    if ind_p > 180
        ind_p = ind_p - 2*72;
    end

    big_S_theta = [S_theta(:,example_run_ind); S_theta(:,example_run_ind); S_theta(:,example_run_ind)];

    inds_consider = ind_p + ntheta + (-35:1:36);

    Sm = mean(sin(bigtheta(inds_consider)).*big_S_theta(inds_consider),'omitnan');
    Cm = mean(cos(bigtheta(inds_consider)).*big_S_theta(inds_consider),'omitnan');
    wavedir_rad = atan2(Sm,Cm);

    % Offset field: the 50% spreading limits, which need the Eq. (3) offset to
    % put the crosswind level below half the peak
    S_k_theta_slice = squeeze(S_k_theta_star(:,:,example_run_ind));
    S_f_theta_slice = squeeze(S_f_theta_star(:,:,example_run_ind));

    % Raw field: Delta, whose ratio does not cancel the offset
    S_k_theta_raw_slice = squeeze(big_S_k_theta(:,:,example_run_ind));
    S_f_theta_raw_slice = squeeze(big_S_f_theta(:,:,example_run_ind));

    % Rotate to the downwave frame and smooth over a periodic-extended axis
    [S_f_theta_raw,theta_rad] = convert_dirspect_to_downwind(S_f_theta_raw_slice(:,inds_consider),bigtheta(inds_consider),wavedir_rad);
    [S_k_theta_raw,~] = convert_dirspect_to_downwind(S_k_theta_raw_slice(:,inds_consider),bigtheta(inds_consider),wavedir_rad);

    [S_f_theta,~] = convert_dirspect_to_downwind(S_f_theta_slice(:,inds_consider),bigtheta(inds_consider),wavedir_rad);
    [S_k_theta,~] = convert_dirspect_to_downwind(S_k_theta_slice(:,inds_consider),bigtheta(inds_consider),wavedir_rad);

    big_S_k_theta_downwave = smoothdata2([S_k_theta S_k_theta S_k_theta],'movmean',{3,3},'includenan');
    big_S_f_theta_downwave = smoothdata2([S_f_theta S_f_theta S_f_theta],'movmean',{5,3},'includenan');

    big_S_k_theta_downwave_raw = smoothdata2([S_k_theta_raw S_k_theta_raw S_k_theta_raw],'movmean',{3,3},'includenan');
    big_S_f_theta_downwave_raw = smoothdata2([S_f_theta_raw S_f_theta_raw S_f_theta_raw],'movmean',{5,3},'includenan');

    % 50% spreading limits from the peak-normalized spreading function. Search
    % restricted to +/-120 deg of the downwave direction: the Pyxis spectra
    % carry an inherent 180-deg ambiguity, so the mirror lobe also exceeds 0.5
    wcols = theta_rad >= -2*pi/3 & theta_rad <= 2*pi/3;
    D_k = big_S_k_theta_downwave(:,mid_cols);
    D_f = big_S_f_theta_downwave(:,mid_cols);
    D_k = D_k(:,wcols)./max(D_k(:,wcols),[],2,'omitnan');
    D_f = D_f(:,wcols)./max(D_f(:,wcols),[],2,'omitnan');

    [Dk_left,Dk_right] = fifty_percent_edges(D_k,theta_rad(wcols));
    [Df_left,Df_right] = fifty_percent_edges(D_f,theta_rad(wcols));

    D_spread_holder_struc(example_run_ind).D_k_limits = [Dk_left Dk_right];
    D_spread_holder_struc(example_run_ind).D_f_limits = [Df_left Df_right];

    % Delta centering: energy-weighted circular mean of the raw field at each
    % scale, smoothed across the scale axis
    dir_peak_k = smoothdata(weighted_center_direction(big_S_k_theta_downwave_raw,theta_rad,ntheta),'movmean',3,'omitnan');
    dir_peak_f = smoothdata(weighted_center_direction(big_S_f_theta_downwave_raw,theta_rad,ntheta),'movmean',3,'omitnan');

    % Periodic-extended direction axis, sorted once so Delta can be sampled at
    % exact angles rather than snapped to the nearest 5 deg bin
    theta_big = [theta_rad(:)-2*pi; theta_rad(:); theta_rad(:)+2*pi];
    [theta_big_sorted,inds_sort] = sort(theta_big);

    S_k_theta_raw_sorted = big_S_k_theta_downwave_raw(:,inds_sort);
    S_f_theta_raw_sorted = big_S_f_theta_downwave_raw(:,inds_sort);

    % Three-point averaging width, one direction bin either side of center
    dtheta_rad = median(diff(theta_rad));
    offsets = [-dtheta_rad 0 dtheta_rad];

    R_mirror = NaN*ones(k_ind_high,1);

    for i = 1:k_ind_high

        angles_0 = dir_peak_k(i) + offsets;
        angles_90 = dir_peak_k(i) + [offsets-pi/2, offsets+pi/2];
        angles_180 = dir_peak_k(i) + pi + offsets;

        S_0 = sample_direction(S_k_theta_raw_sorted(i,:),theta_big_sorted,angles_0);
        S_90 = sample_direction(S_k_theta_raw_sorted(i,:),theta_big_sorted,angles_90);
        S_180 = sample_direction(S_k_theta_raw_sorted(i,:),theta_big_sorted,angles_180);

        % Below four samples per wave period the retrieval splits energy evenly
        % between the true lobe and its mirror, halving the downwind sample.
        % The leakage factor runs from 1 where the lobes are resolved to 2
        % where they are fully ambiguous
        F_leak = ambiguity_leakage_factor(S_180,S_0,R0_upwind);

        delta_k(i,example_run_ind) = (F_leak*S_0-S_90)./(F_leak*S_0+S_90);

        if isfinite(S_0) && S_0 > 0
            R_mirror(i) = S_180/S_0;
        end

    end

    % First wavenumber at which the mirror lobe reaches 90% of the true lobe
    R_mirror_smooth = smoothdata(R_mirror,'movmean',5,'omitnan');
    i_ambiguity = find(R_mirror_smooth > 0.9,1,'first');

    if ~isempty(i_ambiguity)
        k_ambiguity_run(example_run_ind) = k_rad_m_Pyxis(i_ambiguity);
    end

    for i = 1:length(f_Hz_Pyxis)

        angles_0 = dir_peak_f(i) + offsets;
        angles_90 = dir_peak_f(i) + [offsets-pi/2, offsets+pi/2];
        angles_180 = dir_peak_f(i) + pi + offsets;

        S_0 = sample_direction(S_f_theta_raw_sorted(i,:),theta_big_sorted,angles_0);
        S_90 = sample_direction(S_f_theta_raw_sorted(i,:),theta_big_sorted,angles_90);
        S_180 = sample_direction(S_f_theta_raw_sorted(i,:),theta_big_sorted,angles_180);

        F_leak = ambiguity_leakage_factor(S_180,S_0,R0_upwind);

        delta_f(i,example_run_ind) = (F_leak*S_0-S_90)./(F_leak*S_0+S_90);

    end

end

save('../data/directional_spreading_quantification_data.mat','D_spread_holder_struc','-v7.3')

save('../data/delta_with_wind_speed_data.mat','delta_k','delta_f','k_ambiguity_run')

% Ambiguity leakage factor from the measured upwind/downwind ratio R. Returns 1
% where the lobes are resolved and 2 where the mirror matches the true lobe
%
% R0 is the genuine upwind energy fraction, so only the excess above it counts
% as leakage
function F = ambiguity_leakage_factor(S_180,S_0,R0)

if ~isfinite(S_0) || S_0 <= 0 || ~isfinite(S_180)
    F = 1;
    return
end

R = S_180/S_0;

F = 1 + max(R-R0,0)/(1-R0);
F = min(max(F,1),2);

end

% Energy-weighted circular mean direction per scale, over the +/-90 deg window
% about the downwave direction. Continuous in angle, so the Delta sampling
% window does not step by a whole bin as noise moves the maximum
%
% S_tiled is the periodic-extended field; theta_rad is the single-period axis
function dir_rad = weighted_center_direction(S_tiled,theta_rad,ntheta)

n_scale = size(S_tiled,1);

% The field is already rotated so that the run wave direction sits at zero
mid_cols = ntheta + (1:ntheta);
inds_win = abs(theta_rad) <= pi/2;
theta_win = theta_rad(inds_win);
cols_win = mid_cols(inds_win);

dir_rad = NaN*ones(n_scale,1);

for i = 1:n_scale

    S_win = S_tiled(i,cols_win);
    S_win(~isfinite(S_win)) = 0;

    if max(S_win) <= 0
        continue
    end

    Sm = sum(S_win(:).*sin(theta_win(:)));
    Cm = sum(S_win(:).*cos(theta_win(:)));

    dir_rad(i) = atan2(Sm,Cm);

end

end

% Directional spectrum evaluated at exact angles by linear interpolation over
% the periodic-extended axis, then averaged. theta_axis must be ascending, with
% S_row given in the same column order
function S_at = sample_direction(S_row,theta_axis,angles_rad)

inds_good = isfinite(S_row);

if sum(inds_good) < 2
    S_at = NaN;
    return
end

S_at = mean(interp1(theta_axis(inds_good),S_row(inds_good),angles_rad,'linear'),'omitnan');

end
