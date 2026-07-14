
close all;clear;clc

load('../data/global_figure_settings.mat')

% Skip expensive recompute unless explicitly requested
if ~exist('recompute_derived_products','var') || ~recompute_derived_products
    return
end

% Regenerates, per run, the 50% directional spreading limits (theta_halfwidth
% edges, via fifty_percent_edges on the downwave-frame spreading function) and
% the spreading ratio Delta with wavenumber/frequency-dependent centering.
% Replaces the legacy interactive-ROI product with a deterministic method
% consistent with the fig-10 overlay.

smooth_type = 'movmean';

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

S_k_theta_star = (big_S_k_theta-min(big_S_k_theta,[],2,'omitnan'))./(max(big_S_k_theta,[],2,'omitnan')-min(big_S_k_theta,[],2,'omitnan'));
S_f_theta_star = (big_S_f_theta-min(big_S_f_theta,[],2,'omitnan'))./(max(big_S_f_theta,[],2,'omitnan')-min(big_S_f_theta,[],2,'omitnan'));

S_theta = squeeze(sum(S_f_theta_star,1,'omitnan'));

bigtheta = [theta_rad_Pyxis-2*pi; theta_rad_Pyxis; theta_rad_Pyxis+2*pi];

D_spread_holder_struc = struct();

delta_f = NaN*ones(length(f_Hz_Pyxis),nrun);
delta_k = NaN*ones(length(k_rad_m_Pyxis),nrun);

k_ind_high = find(k_rad_m_Pyxis<371,1,'last');
mid_cols = ntheta + (1:ntheta);

for example_run_ind = 1:nrun

    S_f_theta_slice = squeeze(S_f_theta_star(:,:,example_run_ind));
    S_k_theta_slice = squeeze(S_k_theta_star(:,:,example_run_ind));

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

    % Rotate to the downwave frame and smooth over a periodic-extended axis
    [S_f_theta,theta_rad] = convert_dirspect_to_downwind(S_f_theta_slice(:,inds_consider),bigtheta(inds_consider),wavedir_rad);
    [S_k_theta,~] = convert_dirspect_to_downwind(S_k_theta_slice(:,inds_consider),bigtheta(inds_consider),wavedir_rad);

    big_S_k_theta_downwave = smoothdata2([S_k_theta S_k_theta S_k_theta],'movmean',{3,3},'includenan');
    big_S_f_theta_downwave = smoothdata2([S_f_theta S_f_theta S_f_theta],'movmean',{5,3},'includenan');

    % 50% spreading limits from the peak-normalized spreading function.
    % Search restricted to +/-90 deg of the downwave direction: the Pyxis
    % spectra carry an inherent 180-deg ambiguity, so the mirror lobe also
    % exceeds 0.5 — the half-plane window is the deterministic equivalent of
    % the legacy interactive-ROI ambiguity resolution.
    wcols = theta_rad >= -pi/2 & theta_rad <= pi/2;
    D_k = big_S_k_theta_downwave(:,mid_cols);
    D_f = big_S_f_theta_downwave(:,mid_cols);
    D_k = D_k(:,wcols)./max(D_k(:,wcols),[],2,'omitnan');
    D_f = D_f(:,wcols)./max(D_f(:,wcols),[],2,'omitnan');

    [Dk_left,Dk_right] = fifty_percent_edges(D_k,theta_rad(wcols));
    [Df_left,Df_right] = fifty_percent_edges(D_f,theta_rad(wcols));

    D_spread_holder_struc(example_run_ind).D_k_limits = [Dk_left Dk_right];
    D_spread_holder_struc(example_run_ind).D_f_limits = [Df_left Df_right];

    % Delta about the scale-dependent center direction (mean of the 50% limits)
    dir_mean_k = smoothdata(mean([Dk_left Dk_right],2,'omitnan'),'movmean',3,'omitnan');
    dir_mean_f = smoothdata(mean([Df_left Df_right],2,'omitnan'),'movmean',3,'omitnan');

    for i = 1:k_ind_high

        dir_diff = abs(theta_rad-dir_mean_k(i));
        ind_center = find(dir_diff==min(dir_diff,[],'all','omitnan'),1,'first') + ntheta;
        if isempty(ind_center)
            ind_center = ntheta + 36;
        end
        inds_0_deg = ind_center + [-1 0 1];
        inds_90_deg = ind_center + [[-1 0 1]-18,[-1 0 1]+18];
        delta_k(i,example_run_ind) = (mean(big_S_k_theta_downwave(i,inds_0_deg),2,'omitnan')-mean(big_S_k_theta_downwave(i,inds_90_deg),2,'omitnan'))./(mean(big_S_k_theta_downwave(i,inds_0_deg),2,'omitnan')+mean(big_S_k_theta_downwave(i,inds_90_deg),2,'omitnan'));

    end

    for i = 1:length(f_Hz_Pyxis)

        dir_diff = abs(theta_rad-dir_mean_f(i));
        ind_center = find(dir_diff==min(dir_diff,[],'all','omitnan'),1,'first') + ntheta;
        if isempty(ind_center)
            ind_center = ntheta + 36;
        end
        inds_0_deg = ind_center + [-1 0 1];
        inds_90_deg = ind_center + [[-1 0 1]-18,[-1 0 1]+18];
        delta_f(i,example_run_ind) = (mean(big_S_f_theta_downwave(i,inds_0_deg),2,'omitnan')-mean(big_S_f_theta_downwave(i,inds_90_deg),2,'omitnan'))./(mean(big_S_f_theta_downwave(i,inds_0_deg),2,'omitnan')+mean(big_S_f_theta_downwave(i,inds_90_deg),2,'omitnan'));

    end

end

save('../data/directional_spreading_quantification_data.mat','D_spread_holder_struc','-v7.3')

save('../data/delta_with_wind_speed_data.mat','delta_k','delta_f')
