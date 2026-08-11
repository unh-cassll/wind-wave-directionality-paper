
function wavenumber_directional_and_omni_spect(fignum,fsize)


supporting_nc_name = 'data/ASIT2019_supporting_environmental_observations.nc';
short_wave_spectra_nc_name = 'data/ASIT2019_wave_spectra_stats_timeseries_empirical_gain.nc';
long_wave_spectra_nc_name = 'data/ASIT2019_EPSS_directional_spectra.nc';

s = load('data/global_figure_settings.mat');
example_run_ind = s.example_run_ind;

cmap = flipud(spectral(7));
violet = cmap(1,:);
teal = cmap(2,:);
goldenrod = cmap(5,:);
crimson = cmap(7,:);

wdir_deg = ncread(supporting_nc_name,'COARE_Wdir');
wdir_deg = wdir_deg(example_run_ind);
wdir_deg = mod(wdir_deg+180,360);
wdir_rad = pi/180*wdir_deg;

k_rad_m_Pyxis = ncread(short_wave_spectra_nc_name,'k_rad_m');
theta_rad_Pyxis = ncread(short_wave_spectra_nc_name,'theta_rad');
S_k_theta = double(ncread(short_wave_spectra_nc_name,'S_k_theta'));

dir_EPSS = double(ncread(long_wave_spectra_nc_name,'direction'));
k_EWDM = ncread(long_wave_spectra_nc_name,'wavenumber');
F_k_d = ncread(long_wave_spectra_nc_name,'F_k_d');

% Auto-detect direction units (legacy: deg/per-deg, current: rad/per-rad)
if max(abs(dir_EPSS)) > 2*pi
    theta_rad_EPSS = dir_EPSS*pi/180;
    F_k_d = F_k_d*180/pi;
else
    theta_rad_EPSS = dir_EPSS;
end

S_k_theta = squeeze(S_k_theta(:,:,example_run_ind));
F_k_d = squeeze(F_k_d(example_run_ind,:,:)).';   % (wavenumber, direction)

[S_k_theta,theta_rad] = convert_dirspect_to_downwind(S_k_theta,theta_rad_Pyxis,wdir_rad);
% E-PSS/EWDM directional convention matches Pyxis (no -pi offset); the legacy
% F_f_d handling used wdir_rad-pi, which flips the new EWDM spectrum 180 deg
[F_k_d,theta_rad_EPSS_dw] = convert_dirspect_to_downwind(F_k_d,theta_rad_EPSS,wdir_rad);
% Resample E-PSS onto the Pyxis downwind direction grid used downstream
F_k_d = regrid_directional_spectrum(F_k_d,theta_rad_EPSS_dw,theta_rad);
F_k_d(isnan(F_k_d)) = 0;

% EWDM long waves shown only below the k = 2 rad/m handoff to Pyxis direct
keepE = k_EWDM < 2;

dtheta = median(diff(theta_rad_Pyxis));

downwind_mask = 0*theta_rad+1;
upwind_mask = 0*theta_rad+1;
downwind_mask(abs(theta_rad)>pi/2) = NaN;
upwind_mask(abs(theta_rad)<pi/2) = NaN;

klims = [1e-2 1e3];
theta_lims = 180*[-1 1];

k_ticks = 10.^(-2:1:3);

dir_ticks = 180*(-1:0.25:1);

text_x = 0.02;
text_y = 0.94;

figure(fignum);clf
tlayout = tiledlayout(2,1);

nexttile(1)
hold on
pcolor(k_rad_m_Pyxis,(theta_rad-dtheta/2)*180/pi,log10(k_rad_m_Pyxis.^2.*S_k_theta)')
pcolor(k_EWDM(keepE),(theta_rad-dtheta/2)*180/pi,log10(k_EWDM(keepE).^4.*F_k_d(keepE,:))')
plot([1e-10 1e10],0*[1 1],'--','Color',0.5*[1 1 1],'linewidth',2)
hold off
box on
shading('flat')
cbar = colorbar;
set(get(cbar,'Label'),'String','$\mathrm{log_{10}\{B(k,\theta)\}}$\ ','Interpreter','LaTeX')
cbar.Location = 'northoutside';
ylim(theta_lims)
xlim(klims)
clim([-5.5 -2.5])
ax_struc(1).ax = gca;
% ax_struc(1).ax.YDir = 'reverse';
ax_struc(1).ax.XScale = 'log';
ax_struc(1).ax.XTick = k_ticks;
ax_struc(1).ax.YTick = dir_ticks;
ylabel('$\mathrm{\theta-\theta_{\mathrm{wind}}\ [^\circ]}$','Interpreter','latex')
xlabel('$\mathrm{k\ [rad\ m^{-1}]}$','Interpreter','latex')
text(text_x,text_y,'(a)','Color','w','FontSize',fsize,'Units','normalized')

nexttile(2)
hold on

% omni
plot(k_EWDM(keepE),sum(k_EWDM(keepE).^4.*F_k_d(keepE,:),2,'omitnan')*dtheta,'--','Color',violet,'linewidth',3);
h_omni = plot(k_rad_m_Pyxis,k_rad_m_Pyxis.^2.*sum(S_k_theta,2,'omitnan')*dtheta,'-','Color',violet,'linewidth',3);

% downwind
plot(k_EWDM(keepE),sum(cos(theta_rad).^2.*downwind_mask.*k_EWDM(keepE).^4.*F_k_d(keepE,:),2,'omitnan')*dtheta,'--','Color',teal,'linewidth',3)
h_down = plot(k_rad_m_Pyxis,k_rad_m_Pyxis.^2.*sum(cos(theta_rad).^2.*downwind_mask.*S_k_theta,2,'omitnan')*dtheta,'-','Color',teal,'linewidth',3);

% crosswind
plot(k_EWDM(keepE),sum(sin(theta_rad).^2.*k_EWDM(keepE).^4.*F_k_d(keepE,:),2,'omitnan')*dtheta,'--','Color',goldenrod,'linewidth',3)
h_cross = plot(k_rad_m_Pyxis,k_rad_m_Pyxis.^2.*sum(sin(theta_rad).^2.*S_k_theta,2,'omitnan')*dtheta,'-','Color',goldenrod,'linewidth',3);

% upwind
plot(k_EWDM(keepE),sum(cos(theta_rad).^2.*upwind_mask.*k_EWDM(keepE).^4.*F_k_d(keepE,:),2,'omitnan')*dtheta,'--','Color',crimson,'linewidth',3)
h_up = plot(k_rad_m_Pyxis,k_rad_m_Pyxis.^2.*sum(cos(theta_rad).^2.*upwind_mask.*S_k_theta,2,'omitnan')*dtheta,'-','Color',crimson,'linewidth',3);

hold off
box on
xlim(klims)
ylim([1e-6 1e-2])
ax_struc(2).ax = gca;
ax_struc(2).ax.XTick = k_ticks;
ax_struc(2).ax.YTick = 10.^(-6:1:0);
ax_struc(2).ax.XScale = 'log';
ax_struc(2).ax.YScale = 'log';
xlabel('$\mathrm{k\ [rad\ m^{-1}]}$','Interpreter','latex')
ylabel('$\mathrm{B(k)\ [rad]}$','Interpreter','latex')
H = [h_omni h_down h_cross h_up];
L = {'omni','downwind','crosswind','upwind'};
legend(H,L,'Location','southeast')
text(text_x,text_y,'(b)','FontSize',fsize,'Units','normalized')

text(0.25,0.92,'E-PSS','FontSize',fsize,'Units','normalized','HorizontalAlignment','center')
text(0.85,0.92,'direct','FontSize',fsize,'Units','normalized','HorizontalAlignment','center')

tile_cleaner(ax_struc,tlayout)

tlayout.TileSpacing = 'tight';
