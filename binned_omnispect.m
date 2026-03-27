function binned_omnispect(fignum,fsize)

supporting_data_nc_name = 'data/ASIT2019_supporting_environmental_observations.nc';
short_wave_spectra_nc_name = 'data/ASIT2019_wave_spectra_stats_timeseries_empirical_gain.nc';
long_wave_spectra_nc_name = 'data/ASIT2019_EPSS_directional_spectra.nc';

g = 9.81;

f_Hz_EPSS = ncread(long_wave_spectra_nc_name,'frequency');
theta_rad_EPSS = double(ncread(long_wave_spectra_nc_name,'direction'))*pi/180;

f_Hz_Pyxis = ncread(short_wave_spectra_nc_name,'f_Hz');
k_rad_m_Pyxis = ncread(short_wave_spectra_nc_name,'k_rad_m');
theta_rad_Pyxis = ncread(short_wave_spectra_nc_name,'theta_rad');

kappa = 0.4;

EC_U_m_s = ncread(supporting_data_nc_name,'EC_U_m_s');
EC_ustar_m_s = ncread(supporting_data_nc_name,'EC_ustar_m_s');
EC_z_m_above_water = ncread(supporting_data_nc_name,'EC_z_m_above_water');
EC_z0_m = EC_z_m_above_water.*exp(-kappa*EC_U_m_s./EC_ustar_m_s);
EC_U10_m_s = EC_ustar_m_s/kappa.*log(10./EC_z0_m);

SKTHETA = double(ncread(short_wave_spectra_nc_name,'S_k_theta'));
SFTHETA = double(ncread(short_wave_spectra_nc_name,'S_f_theta'));
FFTHETA = double(ncread(long_wave_spectra_nc_name,'F_f_d'))*180/pi;

SKTHETA(SKTHETA==0) = NaN;
SFTHETA(SFTHETA==0) = NaN;
FFTHETA(FFTHETA==0) = NaN;

U_E_m_s = double(ncread(supporting_data_nc_name,'U_E_m_s'));
U_N_m_s = double(ncread(supporting_data_nc_name,'U_N_m_s'));
z_m = double(ncread(supporting_data_nc_name,'z_m'));
tail_flag = false;
k_max = 10;

f_cut = 0.35;

ind_trim = find(f_Hz_Pyxis>f_cut,1,'first');

FFTHETA = FFTHETA(:,:,f_Hz_EPSS<f_cut);
f_Hz_EPSS = f_Hz_EPSS(f_Hz_EPSS<f_cut);

ind_cut = length(f_Hz_EPSS);

water_depth_m = 15;

N = length(EC_U10_m_s);

F_f_block = NaN*ones(length(f_Hz_EPSS)+length(f_Hz_Pyxis(ind_trim:end)),N);
S_f_block = F_f_block;
F_k_block = NaN*ones(1024+length(k_rad_m_Pyxis),N);

for example_run_ind = 1:N

    try

        wdir_deg = double(ncread(supporting_data_nc_name,'COARE_Wdir'));

        wdir_deg = mod(wdir_deg(example_run_ind)+180,360);
        wdir_rad = pi/180*wdir_deg;

        S_k_theta = squeeze(SKTHETA(:,:,example_run_ind));

        S_f_theta = squeeze(SFTHETA(:,:,example_run_ind));
        F_f_theta_m2_Hz_rad = squeeze(FFTHETA(example_run_ind,:,:))';

        [S_f_theta,theta_rad] = convert_dirspect_to_downwind(S_f_theta,theta_rad_Pyxis,wdir_rad);
        [F_f_theta_m2_Hz_rad,~] = convert_dirspect_to_downwind(F_f_theta_m2_Hz_rad,theta_rad_EPSS,wdir_rad);

        [c_disp_EPSS,~] = lindisp_with_current(2*pi*f_Hz_EPSS,water_depth_m,0);
        [c_disp_Pyxis,~] = lindisp_with_current(2*pi*f_Hz_Pyxis,water_depth_m,0);

        k_disp_EPSS = 2*pi*f_Hz_EPSS./c_disp_EPSS;
        k_disp_Pyxis = 2*pi*f_Hz_Pyxis./c_disp_Pyxis;

        dtheta = median(diff(theta_rad_Pyxis));

        S_f_Pyxis = sum(S_f_theta,2,'omitnan')*dtheta;
        F_f_Pyxis = k_disp_Pyxis.^-2.*S_f_Pyxis;
        F_f_EPSS = sum(F_f_theta_m2_Hz_rad,2,'omitnan')*dtheta;

        f_Hz_combined = [f_Hz_EPSS; f_Hz_Pyxis(ind_trim:end)];
        F_f_combined = [F_f_EPSS; F_f_Pyxis(ind_trim:end)];
        S_f_combined = [k_disp_EPSS.^2.*F_f_EPSS; S_f_Pyxis(ind_trim:end)];

        [S_k_theta,~] = convert_dirspect_to_downwind(S_k_theta,theta_rad_Pyxis,wdir_rad);

        Umag_m_s = sqrt(U_E_m_s(example_run_ind,:).^2+U_N_m_s(example_run_ind,:).^2);
        Udir_deg = mod(atan2(U_E_m_s(example_run_ind,:),U_N_m_s(example_run_ind,:))*180/pi+360,360);

        [wave_F_k_theta,k_rad_m,~,~] = directional_Doppler_shift_spectrum(Umag_m_s,Udir_deg,z_m(example_run_ind,:),water_depth_m,f_Hz_EPSS,F_f_theta_m2_Hz_rad,theta_rad,tail_flag,k_max);

        k_rad_m_combined = [k_rad_m; k_rad_m_Pyxis(:)];
        F_k_theta_combined = [wave_F_k_theta; k_rad_m_Pyxis.^-2.*S_k_theta];

        F_k_combined = sum(k_rad_m_combined.*F_k_theta_combined,2,'omitnan')*dtheta;

        F_f_block(:,example_run_ind) = F_f_combined(:);
        S_f_block(:,example_run_ind) = S_f_combined(:);
        F_k_block(:,example_run_ind) = F_k_combined(:);

    end

end

F_f_block(F_f_block==0) = NaN;
S_f_block(S_f_block==0) = NaN;
F_k_block(F_k_block==0) = NaN;

[c_m_s_disp,~] = lindisp_with_current(2*pi*f_Hz_combined,water_depth_m,0);

k_rad_m_disp = 2*pi*f_Hz_combined./c_m_s_disp;

save('data/ASIT2019_combined_frequency_slope_spectra.mat','F_f_block','S_f_block','f_Hz_combined','-v7.3')
save('data/ASIT2019_combined_wavenumber_elevation_spectra.mat','F_k_block','k_rad_m_combined','-v7.3')

%% FIGURES

text_x = 0.05;
text_y = 0.95;
labels = {'(a)','(b)','(c)','(d)','(e)','(f)'};

f_p = sum(f_Hz_combined.*F_f_block.^4,1,'omitnan')./sum(F_f_block.^4,1,'omitnan');
[c_p,~] = lindisp_with_current(2*pi*f_p,water_depth_m,0);

waveage = c_p./EC_ustar_m_s(:);

d_waveage_norm = 10;
waveage_centers = 15:d_waveage_norm:55;
nU = length(waveage_centers);

Ulims = [waveage_centers(1) waveage_centers(end)] + [-1 1]*d_waveage_norm/2;

cmap_binned = (magma(nU));

F_f_binned = NaN*ones(size(F_f_combined,1),nU);
S_f_binned = F_f_binned;
F_k_binned = NaN*ones(size(F_k_combined,1),nU);

for n = 1:nU

    inds_consider = waveage >= waveage_centers(n)-d_waveage_norm/2 & waveage < waveage_centers(n)+d_waveage_norm/2;
    F_f_binned(:,n) = mean(F_f_block(:,inds_consider),2,'omitnan');
    S_f_binned(:,n) = mean(S_f_block(:,inds_consider),2,'omitnan');
    F_k_binned(:,n) = mean(F_k_block(:,inds_consider),2,'omitnan');

end

F_f_binned = fliplr(F_f_binned);
F_k_binned = fliplr(F_k_binned);

f_ind_cut = ind_cut;
k_ind_cut = 1024;

F_f_binned(f_ind_cut,:) = NaN;
F_k_binned(k_ind_cut,:) = NaN;

f_eq = [2e-1 5e-1];
f_sat = [5e-1 1.5e0];

k_eq = [1e-1 1e0];
k_sat = [1e0 2e1];

figure(fignum);clf
tlayout = tiledlayout(2,2);
ax_struc = struct();

nexttile(1)
hold on
loglog(f_Hz_combined,F_f_binned,'k-','linewidth',3)
loglog(f_Hz_combined,F_f_binned,'linewidth',2)
plot(f_eq,1e-2*f_eq.^-4,'k--','linewidth',2.5)
plot(f_sat,5e-3*f_sat.^-5,'k:','linewidth',3)
hold off
colororder(cmap_binned)
colormap(cmap_binned)
clim(Ulims)
xlim([1e-2 2e1])
ylim([1e-10 1e2])
ylabel('F(f) [m^2Hz^{-1}]')
text(mean(f_eq),10.^mean(log10(1e-2*f_eq.^-4))*7,'f^{-4}','FontSize',fsize,'HorizontalAlignment','center')
text(mean(f_sat),10.^mean(log10(5e-3*f_sat.^-5))*7,'f^{-5}','FontSize',fsize,'HorizontalAlignment','center')

nexttile(2)
hold on
loglog(k_rad_m_combined,F_k_binned,'k-','linewidth',3)
loglog(k_rad_m_combined,F_k_binned,'linewidth',2)
plot(k_eq,4e-2*k_eq.^-2.5,'k--','linewidth',2.5)
plot(k_sat,4e-2*k_sat.^-3,'k:','linewidth',3)
hold off
colororder(cmap_binned)
xlim([1e-2 1e3])
ylim([1e-12 1e2])
ylabel('F(k) [m^3]')
text(1.5*10^mean(log10(k_eq)),10.^mean(log10(4e-2*k_eq.^-2.5))*20,'k^{-2.5}','FontSize',fsize,'HorizontalAlignment','center')
text(1.5*10^mean(log10(k_sat)),10.^mean(log10(4e-2*k_sat.^-3))*20,'k^{-3}','FontSize',fsize,'HorizontalAlignment','center')

nexttile(3)
hold on
loglog(f_Hz_combined,k_rad_m_disp.^2.*f_Hz_combined.*F_f_binned,'k-','linewidth',3)
loglog(f_Hz_combined,k_rad_m_disp.^2.*f_Hz_combined.*F_f_binned,'linewidth',2)
plot(f_eq,1.2*3e-2*f_eq.^1,'k--','linewidth',2.5)
plot(f_sat,1.2*1.5e-2*f_sat.^0,'k:','linewidth',3)
hold off
colororder(cmap_binned)
xlim([1e-2 2e1])
ylim([1e-4 5e-2])
xlabel('f [Hz]')
ylabel('(2\pif)^5/g^2F(f) [rad]')
text(mean(f_eq),1.3*10.^mean(log10(3e-2*f_eq.^1))*1.5,'f^{1}','FontSize',fsize,'HorizontalAlignment','center')
text(mean(f_sat),1.3*10.^mean(log10(1.5e-2*f_sat.^0))*1.5,'f^{0}','FontSize',fsize,'HorizontalAlignment','center')

nexttile(4)
hold on
loglog(k_rad_m_combined,k_rad_m_combined.^3.*F_k_binned,'k-','linewidth',3)
loglog(k_rad_m_combined,k_rad_m_combined.^3.*F_k_binned,'linewidth',2)
plot(k_eq,1.8e-2*k_eq.^0.5,'k--','linewidth',2.5)
plot(k_sat,1.7e-2*k_sat.^0,'k:','linewidth',3)
hold off
colororder(cmap_binned)
xlim([1e-2 1e3])
ylim([1e-4 5e-2])
xlabel('k [rad m^{-1}]')
ylabel('k^3F(k) [rad]')
text(10^mean(log10(k_eq)),1.1*10.^mean(log10(1.8e-2*k_eq.^0.5))*1.5,'k^{0.5}','FontSize',fsize,'HorizontalAlignment','center')
text(10^mean(log10(k_sat)),1.1*10.^mean(log10(1.7e-2*k_sat.^0))*1.5,'k^{0}','FontSize',fsize,'HorizontalAlignment','center')

for n = 1:4
    ax_struc(n).ax = nexttile(n);
    ax_struc(n).ax.XScale = 'log';
    ax_struc(n).ax.YScale = 'log';
    ax_struc(n).ax.XTick = 10.^(-3:1:3);
    box on
end

ax_struc(2).ax.YAxisLocation = 'right';
ax_struc(4).ax.YAxisLocation = 'right';

for n = 1:4
    nexttile(n)
    text(text_x,text_y,labels{n},'Color',[0 0 0],'Units','normalized','FontSize',fsize,'HorizontalAlignment','center')
    clim(Ulims)
end

tlayout.TileSpacing = 'tight';

for n = 1:2
    ax_struc(n).ax.XTickLabel = '';
    ax_struc(n).ax.YTick = 10.^(-20:5:5);
end

cbar = colorbar(ax_struc(2).ax);

cbar.Layout.Tile = 'north';
cbar.Layout.TileSpan = [2 2];
cbar.Ticks = Ulims(1):d_waveage_norm:Ulims(2);
cbar.TickLabels = flipud(cbar.TickLabels);
cbar.Direction = 'reverse';
set(get(cbar,'Label'),'String','c_p/u_*')
