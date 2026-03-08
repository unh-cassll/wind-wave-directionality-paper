
close all;clear;clc

load('../data/global_figure_settings.mat')

smooth_type = 'movmean';

f_smooth = 5;
k_smooth = 5;
ang_smooth = 3;

short_wave_nc_name = '../data/ASIT2019_wave_spectra_stats_timeseries_empirical_gain.nc';
long_wave_nc_name = '../data/ASIT2019_EPSS_directional_spectra.nc';
supporting_nc_name = '../data/ASIT2019_supporting_environmental_observations.nc';

g = 9.81;

S_k_theta_block = double(ncread(short_wave_nc_name,'S_k_theta'));
S_f_theta_block = double(ncread(short_wave_nc_name,'S_f_theta'));

f_Hz_Pyxis = ncread(short_wave_nc_name,'f_Hz');
k_rad_m_Pyxis = ncread(short_wave_nc_name,'k_rad_m');
theta_rad_Pyxis = ncread(short_wave_nc_name,'theta_rad');

big_S_k_theta = [S_k_theta_block S_k_theta_block S_k_theta_block];
big_S_f_theta = [S_f_theta_block S_f_theta_block S_f_theta_block];

for i = 1:190
    big_S_k_theta(:,:,i) = smoothdata2(squeeze(big_S_k_theta(:,:,i)),smooth_type,{k_smooth,ang_smooth},'omitnan');
    big_S_f_theta(:,:,i) = smoothdata2(squeeze(big_S_f_theta(:,:,i)),smooth_type,{f_smooth,ang_smooth},'omitnan');
end

f_low = 0.6;

big_S_k_theta(k_rad_m_Pyxis<k_low,:,:) = NaN;
big_S_f_theta(f_Hz_Pyxis<f_low,:,:) = NaN;

S_k_theta_star = (big_S_k_theta-min(big_S_k_theta,[],2,'omitnan'))./(max(big_S_k_theta,[],2,'omitnan')-min(big_S_k_theta,[],2,'omitnan'));
S_f_theta_star = (big_S_f_theta-min(big_S_f_theta,[],2,'omitnan'))./(max(big_S_f_theta,[],2,'omitnan')-min(big_S_f_theta,[],2,'omitnan'));

ntheta = length(theta_rad_Pyxis);

S_theta = squeeze(sum(S_f_theta_star,1,'omitnan'));

bigtheta = [theta_rad_Pyxis-2*pi; theta_rad_Pyxis; theta_rad_Pyxis+2*pi];

kappa = 0.4;

EC_U_m_s = ncread(supporting_nc_name,'EC_U_m_s');
EC_ustar_m_s = ncread(supporting_nc_name,'EC_ustar_m_s');
EC_z_m_above_water = ncread(supporting_nc_name,'EC_z_m_above_water');
EC_z0_m = EC_z_m_above_water.*exp(-kappa*EC_U_m_s./EC_ustar_m_s);
EC_U10_m_s = EC_ustar_m_s/kappa.*log(10./EC_z0_m);

D_spread_holder_struc = struct();

delta_f = NaN*ones(length(f_Hz_Pyxis),length(EC_ustar_m_s));
delta_k = NaN*ones(length(k_rad_m_Pyxis),length(EC_ustar_m_s));

inds_0_deg = 35:37;
inds_90_deg = [17:19 53:55];

%%

% run_start = 1;
%
% for example_run_ind = run_start:size(S_f_theta_block,3)
%
%     S_f_theta_slice = squeeze(S_f_theta_star(:,:,example_run_ind));
%     S_k_theta_slice = squeeze(S_k_theta_star(:,:,example_run_ind));
%
%     down_up_vals_k = NaN*ones(length(k_rad_m_Pyxis),2);
%     ind_p_k = NaN*k_rad_m_Pyxis;
%     i_start = find(~isnan(sum(big_S_k_theta,2)),1,'first');
%     for i = i_start:length(k_rad_m_Pyxis)
%         S_theta_slice = S_k_theta_slice(i,:)';
%         ind_p = find(S_theta_slice==max(S_theta_slice,[],'all','omitnan'),1,'first');
%         if ind_p < 37
%             ind_p = ind_p + 72;
%         end
%         if ind_p > 180
%             ind_p = ind_p - 72;
%         end
%         inds_consider = ind_p + (-36:1:36);
%         Sm = mean(sin(bigtheta(inds_consider)).*S_theta_slice(inds_consider),'omitnan');
%         Cm = mean(cos(bigtheta(inds_consider)).*S_theta_slice(inds_consider),'omitnan');
%         wavedir_rad = atan2(Sm,Cm);
%         ind_p_new = find(bigtheta>wavedir_rad,1,'first') + 72;
%         ind_p_k(i) = ind_p_new;
%         down_up_vals_k(i,:) = [mean(S_theta_slice(ind_p_new + (-1:1:1))) mean(S_theta_slice(ind_p_new + (-1:1:1) - 36))];
%     end
%
%     ind_p = find(S_theta(:,example_run_ind)==max(S_theta(:,example_run_ind),[],'all','omitnan'),1,'first');
%     if ind_p < 37
%         ind_p = ind_p + 72;
%     end
%     if ind_p > 180
%         ind_p = ind_p - 2*72;
%     end
%
%     big_S_theta = [S_theta(:,example_run_ind); S_theta(:,example_run_ind); S_theta(:,example_run_ind)];
%
%     inds_consider = ind_p + ntheta + (-35:1:36);
%
%     Sm = mean(sin(bigtheta(inds_consider)).*big_S_theta(inds_consider),'omitnan');
%     Cm = mean(cos(bigtheta(inds_consider)).*big_S_theta(inds_consider),'omitnan');
%     wavedir_rad = atan2(Sm,Cm);
%
%     [S_f_theta,theta_rad] = convert_dirspect_to_downwind(S_f_theta_slice(:,inds_consider),bigtheta(inds_consider),wavedir_rad);
%     [S_k_theta,~] = convert_dirspect_to_downwind(S_k_theta_slice(:,inds_consider),bigtheta(inds_consider),wavedir_rad);
%
%     big_S_k_theta_downwave = [S_k_theta S_k_theta S_k_theta];
%     big_S_f_theta_downwave = [S_f_theta S_f_theta S_f_theta];
%     big_theta_downwind = [theta_rad-2*pi theta_rad theta_rad+2*pi];
%
%     big_S_k_theta_downwave_copy = big_S_k_theta_downwave;
%     big_S_f_theta_downwave_copy = big_S_f_theta_downwave;
%
%     % pause()
%
%     figure(11);clf
%     set(gcf,'Position',[1500 600 1100 550])
%     imagesc(big_theta_downwind,f_Hz_Pyxis,big_S_f_theta_downwave_copy)
%     colorbar
%     clim([0.1 1])
%     box on
%     grid off
%     view([0 -90])
%     title(num2str(example_run_ind))
%
%     roi_directional_extrema_core = grab_spectrum_roi(gca,f_Hz_Pyxis);
%
%     % roi_directional_extrema_ambiguous = grab_spectrum_roi(gca,f_Hz_Pyxis);
%
%     [core_spect,ambiguous_spect] = combine_core_and_ambiguous_spectra_portions(big_S_f_theta_downwave_copy,big_theta_downwind,f_Hz_Pyxis,roi_directional_extrema_core,roi_directional_extrema_core);
%
%     % core_spect(core_spect < 0.5) = 0;
%     % ambiguous_spect(ambiguous_spect < 0.5) = 0;
%
%     ambiguous_mask = 0*core_spect + 1;
%     ambiguous_mask(ambiguous_spect>0) = 0.5;
%
%     combined_spect = (core_spect + ambiguous_spect).*ambiguous_mask;
%
%     max_ratio = smoothdata(max(core_spect,[],2,'omitnan'),'movmean',3)./smoothdata(max(combined_spect,[],2,'omitnan'),'movmean',3);
%     max_ratio(sum(ambiguous_spect,2,'omitnan')==0) = 1;
%
%     combined_spect = combined_spect.*max_ratio;
%
%     figure(22);clf
%     set(gcf,'Position',[50 600 1100 550])
%     imagesc(big_theta_downwind,f_Hz_Pyxis,core_spect)
%     colorbar
%     xlim([-180 180]*pi/180)
%     clim([0.5 1])
%     box on
%     grid off
%     view([0 -90])
%     title(num2str(example_run_ind))
%
%     figure(33);clf
%     set(gcf,'Position',[50 1150 1100 550])
%     imagesc(big_theta_downwind,f_Hz_Pyxis,ambiguous_spect)
%     colorbar
%     xlim([-180 180]*pi/180)
%     clim([0.5 1])
%     box on
%     grid off
%     view([0 -90])
%     title(num2str(example_run_ind))
%
%     combined_spect = core_spect;
%     big_S_f_theta_downwave_copy = smoothdata2(combined_spect,'movmean',{3,3},'omitnan');
%
%     figure(44);clf
%     set(gcf,'Position',[50 50 1100 550])
%     imagesc(big_theta_downwind,f_Hz_Pyxis,big_S_f_theta_downwave_copy)
%     colorbar
%     xlim([-180 180]*pi/180)
%     clim([0.5 1])
%     box on
%     grid off
%     view([0 -90])
%     title(num2str(example_run_ind))
%
%     center_inds = size(big_S_f_theta_downwave_copy,2)/2 + (-35:1:36);
%     S_f_theta_downwind = big_S_f_theta_downwave_copy(:,center_inds);
%
%     figure(10);clf
%     set(gcf,'Position',[2600 600 550 550])
%     semilogx(down_up_vals_k(:,1)./down_up_vals_k(:,2),k_rad_m_Pyxis,'linewidth',2)
%     xlim([1 nanmax(down_up_vals_k(:,1)./down_up_vals_k(:,2))])
%     ylim([0 550])
%
%     figure(11);clf
%     set(gcf,'Position',[1500 600 1100 550])
%     imagesc(big_theta_downwind,k_rad_m_Pyxis,big_S_k_theta_downwave_copy)
%     ylim([0 550])
%     colorbar
%     clim([0.1 1])
%     box on
%     grid off
%     view([0 -90])
%     title(num2str(example_run_ind))
%
%     roi_directional_extrema_core = grab_spectrum_roi(gca,k_rad_m_Pyxis);
%
%     roi_directional_extrema_ambiguous = grab_spectrum_roi(gca,k_rad_m_Pyxis);
%
%     [core_spect,ambiguous_spect] = combine_core_and_ambiguous_spectra_portions(big_S_k_theta_downwave_copy,big_theta_downwind,k_rad_m_Pyxis,roi_directional_extrema_core,roi_directional_extrema_ambiguous);
%
%     % core_spect(core_spect < 0.5) = 0;
%     % ambiguous_spect(ambiguous_spect < 0.5) = 0;
%     %
%     % ambiguous_mask = 0*core_spect + 1;
%     % ambiguous_mask(ambiguous_spect>0) = 0.5;
%
%     % combined_spect = (core_spect + ambiguous_spect).*ambiguous_mask;
%
%     max_core = smoothdata(max(core_spect,[],2,'omitnan'),'movmean',3);
%     max_ambiguous = smoothdata(max(ambiguous_spect,[],2,'omitnan'),'movmean',3);
%     max_core(max_ambiguous<0.1) = 0.5;
%     max_ambiguous(max_ambiguous<0.1) = 1;
%
%     combined_spect = smoothdata2((core_spect./max_core + ambiguous_spect./max_ambiguous)/2,'movmean',{5,5},'omitnan');
%
%     figure(22);clf
%     set(gcf,'Position',[50 600 1100 550])
%     imagesc(big_theta_downwind,k_rad_m_Pyxis,core_spect)
%     colorbar
%     xlim([-180 180]*pi/180)
%     clim([0.5 1])
%     box on
%     grid off
%     view([0 -90])
%     title(num2str(example_run_ind))
%
%     figure(33);clf
%     set(gcf,'Position',[50 1150 1100 550])
%     imagesc(big_theta_downwind,k_rad_m_Pyxis,ambiguous_spect)
%     colorbar
%     xlim([-180 180]*pi/180)
%     clim([0.5 1])
%     box on
%     grid off
%     view([0 -90])
%     title(num2str(example_run_ind))
%
%     % combined_spect = core_spect;
%     big_S_k_theta_downwave_copy = smoothdata2(combined_spect,'movmean',{3,3},'omitnan');
%
%     figure(44);clf
%     set(gcf,'Position',[50 50 1100 550])
%     imagesc(big_theta_downwind,k_rad_m_Pyxis,big_S_k_theta_downwave_copy)
%     colorbar
%     xlim([-180 180]*pi/180)
%     clim([0.5 1])
%     box on
%     grid off
%     view([0 -90])
%     title(num2str(example_run_ind))
%
%     center_inds = size(big_S_k_theta_downwave_copy,2)/2 + (-35:1:36);
%     S_k_theta_downwind = big_S_k_theta_downwave_copy(:,center_inds);
%
%     % pause()
%
%     D_f_limits = NaN*ones(length(f_Hz_Pyxis),2);
%     D_f_p = NaN*f_Hz_Pyxis;
%     for i = 1:length(f_Hz_Pyxis)
%         D_slice = S_f_theta_downwind(i,:);
%         ind_p = find(D_slice==max(D_slice,[],'all','omitnan'),1,'first');
%
%         inds_left = 1:ind_p;
%         inds_right = ind_p:72;
%         D_slice_left = fliplr(D_slice(inds_left));
%         D_slice_right = (D_slice(inds_right));
%         ind_left = ind_p - find(D_slice_left>0.5,1,'last');
%         ind_right = ind_p + find(D_slice_right>0.5,1,'last');
%         try
%             D_f_limits(i,:) = theta_rad([ind_left ind_right]);
%             D_f_p(i) = theta_rad(ind_p);
%         end
%     end
%
%     D_k_limits = NaN*ones(length(k_rad_m_Pyxis),2);
%     D_k_p = NaN*k_rad_m_Pyxis;
%     for i = 1:length(k_rad_m_Pyxis)
%         D_slice = S_k_theta_downwind(i,:);
%         ind_p = find(D_slice==max(D_slice,[],'all','omitnan'),1,'first');
%
%         inds_left = 1:ind_p;
%         inds_right = ind_p:72;
%         D_slice_left = fliplr(D_slice(inds_left));
%         D_slice_right = (D_slice(inds_right));
%         ind_left = ind_p - find(D_slice_left>0.5,1,'last');
%         ind_right = ind_p + find(D_slice_right>0.5,1,'last');
%         try
%             D_k_limits(i,:) = theta_rad([ind_left ind_right]);
%             D_k_p(i) = theta_rad(ind_p);
%         end
%     end
%
%     D_spread_holder_struc(example_run_ind).D_f_limits = D_f_limits;
%     D_spread_holder_struc(example_run_ind).D_k_limits = D_k_limits;
%
%     S_k_theta = smoothdata2(S_k_theta_downwind,'movmean',{3,3},'includenan');
%     S_f_theta = smoothdata2(S_f_theta_downwind,'movmean',{5,3},'includenan');
%     % S_k_theta = S_k_theta_downwind;
%     % S_f_theta = S_f_theta_downwind;
%
%     delta_k(:,example_run_ind) = (mean(S_k_theta(:,inds_0_deg),2,'omitnan')-mean(S_k_theta(:,inds_90_deg),2,'omitnan'))./(mean(S_k_theta(:,inds_0_deg),2,'omitnan')+mean(S_k_theta(:,inds_90_deg),2,'omitnan'));
%     delta_f(:,example_run_ind) = (mean(S_f_theta(:,inds_0_deg),2,'omitnan')-mean(S_f_theta(:,inds_90_deg),2,'omitnan'))./(mean(S_f_theta(:,inds_0_deg),2,'omitnan')+mean(S_f_theta(:,inds_90_deg),2,'omitnan'));
%
%
% end
%
% save('../data/directional_spreading_quantification_data.mat','D_spread_holder_struc','-v7.3')
%
% save('../data/delta_with_wind_speed_data.mat','delta_k','delta_f')

%% Update: make Delta computation dynamic around center direction

delta_f = NaN*ones(length(f_Hz_Pyxis),length(EC_ustar_m_s));
delta_k = NaN*ones(length(k_rad_m_Pyxis),length(EC_ustar_m_s));

k_ind_high = find(k_rad_m_Pyxis<371,1,'last');

load('../data/directional_spreading_quantification_data.mat')

for example_run_ind = 1:190

    D_k_limits = D_spread_holder_struc(example_run_ind).D_k_limits;
    D_f_limits = D_spread_holder_struc(example_run_ind).D_f_limits;


    S_f_theta_slice = squeeze(S_f_theta_star(:,:,example_run_ind));
    S_k_theta_slice = squeeze(S_k_theta_star(:,:,example_run_ind));

    down_up_vals_k = NaN*ones(length(k_rad_m_Pyxis),2);
    ind_p_k = NaN*k_rad_m_Pyxis;
    i_start = find(~isnan(sum(big_S_k_theta,2)),1,'first');
    for i = i_start:length(k_rad_m_Pyxis)
        S_theta_slice = S_k_theta_slice(i,:)';
        ind_p = find(S_theta_slice==max(S_theta_slice,[],'all','omitnan'),1,'first');
        if ind_p < 37
            ind_p = ind_p + 72;
        end
        if ind_p > 180
            ind_p = ind_p - 72;
        end
        inds_consider = ind_p + (-36:1:36);
        Sm = mean(sin(bigtheta(inds_consider)).*S_theta_slice(inds_consider),'omitnan');
        Cm = mean(cos(bigtheta(inds_consider)).*S_theta_slice(inds_consider),'omitnan');
        wavedir_rad = atan2(Sm,Cm);
        ind_p_new = find(bigtheta>wavedir_rad,1,'first') + 72;
        ind_p_k(i) = ind_p_new;
        down_up_vals_k(i,:) = [mean(S_theta_slice(ind_p_new + (-1:1:1))) mean(S_theta_slice(ind_p_new + (-1:1:1) - 36))];
    end

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

    [S_f_theta,theta_rad] = convert_dirspect_to_downwind(S_f_theta_slice(:,inds_consider),bigtheta(inds_consider),wavedir_rad);
    [S_k_theta,~] = convert_dirspect_to_downwind(S_k_theta_slice(:,inds_consider),bigtheta(inds_consider),wavedir_rad);

    big_S_k_theta_downwave = smoothdata2([S_k_theta S_k_theta S_k_theta],'movmean',{3,3},'includenan');
    big_S_f_theta_downwave = smoothdata2([S_f_theta S_f_theta S_f_theta],'movmean',{5,3},'includenan');
    big_theta_downwind = [theta_rad-2*pi theta_rad theta_rad+2*pi];

    S_k_theta = big_S_k_theta_downwave;
    S_f_theta = big_S_f_theta_downwave;

    dir_mean_k = smoothdata(mean(D_k_limits,2,'omitnan'),'movmean',3,'omitnan');
    dir_mean_f = smoothdata(mean(D_f_limits,2,'omitnan'),'movmean',3,'omitnan');

    for i = 1:k_ind_high

        dir_diff = abs(theta_rad-dir_mean_k(i));
        ind_center = find(dir_diff==min(dir_diff,[],'all','omitnan'),1,'first') + 72;
        if isempty(ind_center)
            ind_center = 72 + 36;
        end
        inds_0_deg = ind_center + [-1 0 1];
        inds_90_deg = ind_center + [[-1 0 1]-18,[-1 0 1]+18];
        delta_k(i,example_run_ind) = (mean(S_k_theta(i,inds_0_deg),2,'omitnan')-mean(S_k_theta(i,inds_90_deg),2,'omitnan'))./(mean(S_k_theta(i,inds_0_deg),2,'omitnan')+mean(S_k_theta(i,inds_90_deg),2,'omitnan'));

    end

    for i = 1:length(f_Hz_Pyxis)

        dir_diff = abs(theta_rad-dir_mean_f(i));
        ind_center = find(dir_diff==min(dir_diff,[],'all','omitnan'),1,'first') + 72;
        if isempty(ind_center)
            ind_center = 72 + 36;
        end
        inds_0_deg = ind_center + [-1 0 1];
        inds_90_deg = ind_center + [[-1 0 1]-18,[-1 0 1]+18];
        delta_f(i,example_run_ind) = (mean(S_f_theta(i,inds_0_deg),2,'omitnan')-mean(S_f_theta(i,inds_90_deg),2,'omitnan'))./(mean(S_f_theta(i,inds_0_deg),2,'omitnan')+mean(S_f_theta(i,inds_90_deg),2,'omitnan'));

    end

end

% figure(1);pcolor(big_theta_downwind*180/pi,k_rad_m_Pyxis,log10(S_k_theta));shading('flat');colorbar;ax=gca;ax.YScale='log';xlim([-180 180])
% figure(2);semilogy(delta_k,k_rad_m_Pyxis,'linewidth',2)
% 
% figure(11);pcolor(big_theta_downwind*180/pi,f_Hz_Pyxis,log10(S_f_theta));shading('flat');colorbar;ax=gca;ax.YScale='log';xlim([-180 180])
% figure(12);semilogy(delta_f,f_Hz_Pyxis,'linewidth',2)

save('../data/delta_with_wind_speed_data.mat','delta_k','delta_f')