
function wavenumber_directional_and_omni_spect(fignum,fsize)


supporting_data_nc_name = 'data/ASIT2019_supporting_environmental_observations.nc';
short_wave_spectra_nc_name = 'data/ASIT2019_wave_spectra_stats_timeseries_empirical_gain.nc';
long_wave_spectra_nc_name = 'data/ASIT2019_EPSS_directional_spectra.nc';

s = load('data/global_figure_settings.mat');
example_run_ind = s.example_run_ind;

cmap = flipud(spectral(7));
violet = cmap(1,:);
teal = cmap(2,:);
goldenrod = cmap(5,:);
crimson = cmap(7,:);

wdir_deg = ncread(supporting_data_nc_name,'COARE_Wdir');
wdir_deg = wdir_deg(example_run_ind);
wdir_deg = mod(wdir_deg+180,360);
wdir_rad = pi/180*wdir_deg;

k_rad_m_Pyxis = ncread(short_wave_spectra_nc_name,'k_rad_m');
theta_rad_Pyxis = ncread(short_wave_spectra_nc_name,'theta_rad');
S_k_theta = double(ncread(short_wave_spectra_nc_name,'S_k_theta'));

f_Hz_EPSS = ncread(long_wave_spectra_nc_name,'frequency');
theta_deg_EPSS = double(ncread(long_wave_spectra_nc_name,'direction'));
F_f_theta_m2_Hz_deg = ncread(long_wave_spectra_nc_name,'F_f_d');
F_f_theta_m2_Hz_rad = F_f_theta_m2_Hz_deg*180/pi;

S_k_theta = squeeze(S_k_theta(:,:,example_run_ind));
F_f_theta_m2_Hz_rad = squeeze(F_f_theta_m2_Hz_rad(example_run_ind,:,:))';

[S_k_theta,theta_rad] = convert_dirspect_to_downwind(S_k_theta,theta_rad_Pyxis,wdir_rad);
[F_f_theta_m2_Hz_rad,~] = convert_dirspect_to_downwind(F_f_theta_m2_Hz_rad,theta_deg_EPSS*pi/180,wdir_rad-pi);
F_f_theta_m2_Hz_rad(isnan(F_f_theta_m2_Hz_rad)) = 0;

inds_trim = f_Hz_EPSS < 0.35;
f_Hz_EPSS = f_Hz_EPSS(inds_trim);
F_f_theta_m2_Hz_rad = F_f_theta_m2_Hz_rad(inds_trim,:);

dtheta = median(diff(theta_rad_Pyxis));

downwind_mask = 0*theta_rad+1;
upwind_mask = 0*theta_rad+1;
downwind_mask(abs(theta_rad)>pi/2) = NaN;
upwind_mask(abs(theta_rad)<pi/2) = NaN;

U_E_m_s = ncread(supporting_data_nc_name,'U_E_m_s');
U_N_m_s = ncread(supporting_data_nc_name,'U_N_m_s');
z_m = ncread(supporting_data_nc_name,'z_m');

U_E_m_s = U_E_m_s(example_run_ind,:);
U_N_m_s = U_N_m_s(example_run_ind,:);
z_m = z_m(example_run_ind,:);

water_depth_m = max(abs(z_m));

Umag_m_s = sqrt(U_E_m_s.^2+U_N_m_s.^2);
Udir_deg = mod(atan2(U_E_m_s,U_N_m_s)*180/pi+360,360);

tail_flag = false;
k_max = 10;

klims = [1e-2 1e3];
theta_lims = 180*[-1 1];

k_ticks = 10.^(-2:1:3);

dir_ticks = 180*(-1:0.25:1);

[wave_F_k_theta,k_rad_m,~,~] = directional_Doppler_shift_spectrum(Umag_m_s,Udir_deg,z_m,water_depth_m,f_Hz_EPSS,F_f_theta_m2_Hz_rad,theta_rad,tail_flag,k_max);

text_x = 0.02;
text_y = 0.94;

figure(fignum);clf
tlayout = tiledlayout(2,1);

nexttile(1)
hold on
pcolor(k_rad_m_Pyxis,(theta_rad-dtheta/2)*180/pi,log10(k_rad_m_Pyxis.^2.*S_k_theta)')
pcolor(k_rad_m,(theta_rad-dtheta/2)*180/pi,log10(k_rad_m.^4.*wave_F_k_theta)')
plot([1e-10 1e10],0*[1 1],'--','Color',0.5*[1 1 1],'linewidth',2)
hold off
box on
shading('flat')
cbar = colorbar;
set(get(cbar,'Label'),'String','$\mathrm{log_{10}}\{B(k,\theta)\}$\ ','Interpreter','LaTeX')
cbar.Location = 'northoutside';
ylim(theta_lims)
xlim(klims)
clim([-5.5 -2.5])
ax_struc(1).ax = gca;
ax_struc(1).ax.YDir = 'reverse';
ax_struc(1).ax.XScale = 'log';
ax_struc(1).ax.XTick = k_ticks;
ax_struc(1).ax.YTick = dir_ticks;
ylabel('\theta [\circ]')
xlabel('k [rad m^{-1}]')
text(text_x,text_y,'(a)','Color','w','FontSize',fsize,'Units','normalized')

nexttile(2)
hold on

% omni
plot(k_rad_m,sum(k_rad_m.^4.*wave_F_k_theta,2,'omitnan')*dtheta,'--','Color',violet,'linewidth',3);
h_omni = plot(k_rad_m_Pyxis,k_rad_m_Pyxis.^2.*sum(S_k_theta,2,'omitnan')*dtheta,'-','Color',violet,'linewidth',3);

% downwind
plot(k_rad_m,sum(cos(theta_rad).^2.*downwind_mask.*k_rad_m.^4.*wave_F_k_theta,2,'omitnan')*dtheta,'--','Color',teal,'linewidth',3)
h_down = plot(k_rad_m_Pyxis,k_rad_m_Pyxis.^2.*sum(cos(theta_rad).^2.*downwind_mask.*S_k_theta,2,'omitnan')*dtheta,'-','Color',teal,'linewidth',3);

% crosswind
plot(k_rad_m,sum(sin(theta_rad).^2.*k_rad_m.^4.*wave_F_k_theta,2,'omitnan')*dtheta,'--','Color',goldenrod,'linewidth',3)
h_cross = plot(k_rad_m_Pyxis,k_rad_m_Pyxis.^2.*sum(sin(theta_rad).^2.*S_k_theta,2,'omitnan')*dtheta,'-','Color',goldenrod,'linewidth',3);

% upwind
plot(k_rad_m,sum(cos(theta_rad).^2.*upwind_mask.*k_rad_m.^4.*wave_F_k_theta,2,'omitnan')*dtheta,'--','Color',crimson,'linewidth',3)
h_up = plot(k_rad_m_Pyxis,k_rad_m_Pyxis.^2.*sum(cos(theta_rad).^2.*upwind_mask.*S_k_theta,2,'omitnan')*dtheta,'-','Color',crimson,'linewidth',3);

hold off
box on
xlim(klims)
ylim([1e-5 3e-2])
ax_struc(2).ax = gca;
ax_struc(2).ax.XTick = k_ticks;
ax_struc(2).ax.YTick = 10.^(-6:1:0);
ax_struc(2).ax.XScale = 'log';
ax_struc(2).ax.YScale = 'log';
xlabel('k [rad m^{-1}]')
ylabel('B(k) [rad]')
H = [h_omni h_down h_cross h_up];
L = {'omni','downwind','crosswind','upwind'};
legend(H,L,'Location','southeast')
text(text_x,text_y,'(b)','FontSize',fsize,'Units','normalized')

text(0.25,0.92,'E-PSS','FontSize',fsize,'Units','normalized','HorizontalAlignment','center')
text(0.75,0.92,'direct','FontSize',fsize,'Units','normalized','HorizontalAlignment','center')

tile_cleaner(ax_struc,tlayout)

tlayout.TileSpacing = 'tight';