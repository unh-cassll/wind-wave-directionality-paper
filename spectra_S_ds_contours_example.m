%
function spectra_S_ds_contours_example(fignum,fsize)

supporting_data_nc_name = 'data/ASIT2019_supporting_environmental_observations.nc';
short_wave_spectra_nc_name = 'data/ASIT2019_wave_spectra_stats_timeseries_empirical_gain.nc';
long_wave_spectra_nc_name = 'data/ASIT2019_EPSS_directional_spectra.nc';

s = load('data/global_figure_settings.mat');
example_run_ind = s.example_run_ind;

wdir_deg = ncread(supporting_data_nc_name,'COARE_Wdir');
wdir_deg = wdir_deg(example_run_ind);
wdir_deg = mod(wdir_deg+180,360);

EC_ustar_m_s = ncread(supporting_data_nc_name,'EC_ustar_m_s');
EC_ustar_m_s_particular = EC_ustar_m_s(example_run_ind);

load('data/ASIT2019_Lambda_c_theta_Pyxis_runs.mat')

cmap = flipud(spectral(7));
crimson = cmap(7,:);

cmap = viridis(7);
cerulean = cmap(3,:);

IPX_offsets_deg = [38 208];

g = 9.806;
sigma = 0.072;
rho_w = 1020;

k = 0.01:0.01:300;
water_depth_m = 15;

omega_disp = sqrt((g*k+sigma/rho_w*k.^3).*tanh(k*water_depth_m));

c_disp = omega_disp./k;

k_of_c = fit(c_disp(:),k(:),'spline');

dir_ticks = 180*(-1:0.25:1);

S_k_theta = ncread(short_wave_spectra_nc_name,'S_k_theta');
S_f_theta = ncread(short_wave_spectra_nc_name,'S_f_theta');
k_rad_m = ncread(short_wave_spectra_nc_name,'k_rad_m');
theta_rad = ncread(short_wave_spectra_nc_name,'theta_rad');
f_Hz = ncread(short_wave_spectra_nc_name,'f_Hz');

F_f_d = double(ncread(long_wave_spectra_nc_name,'F_f_d'));
frequency = double(ncread(long_wave_spectra_nc_name,'frequency'));
direction = double(ncread(long_wave_spectra_nc_name,'direction'));

nc_name_supporting = 'data/ASIT2019_supporting_environmental_observations.nc';

U_sfc_mag_m_s = ncread(nc_name_supporting,'U_sfc_mag_m_s');

freq_spect_range_limits = load('data/frequency_spect_range_limits.mat');
wavenumber_spect_range_limits = load('data/wavenumber_spect_range_limits.mat');

freq_spect_range_limits.f_eq_start = freq_spect_range_limits.f_eq_start(:);
freq_spect_range_limits.f_eq_end = freq_spect_range_limits.f_eq_end(:);
freq_spect_range_limits.f_sat_end = freq_spect_range_limits.f_sat_end(:);
wavenumber_spect_range_limits.k_sat_end = wavenumber_spect_range_limits.k_sat_end(:);

S_f_particular = squeeze(S_f_theta(:,:,example_run_ind));

Ff = trapz(direction,squeeze(F_f_d(example_run_ind,:,:)));
Sf = trapz(theta_rad,S_f_particular')';

f_combined = [frequency(frequency>frequency(1) & frequency <= 0.35); f_Hz(f_Hz>0.4 & f_Hz <=5)];
Ff_combined = [Ff(frequency>frequency(1) & frequency <= 0.35)'; ((2*pi*f_Hz(f_Hz>0.4 & f_Hz <=5)).^2/g).^-2.*Sf(f_Hz>0.4 & f_Hz <=5)];
T_E = trapz(f_combined,f_combined.^-1.*Ff_combined)/trapz(f_combined,Ff_combined);
f_E = T_E^-1;

[c_E,~] = lindisp_with_current(2*pi*f_E,water_depth_m,0);

[c_direct_disp,cg_direct_disp] = lindisp_with_current(2*pi*f_Hz,water_depth_m,U_sfc_mag_m_s(example_run_ind));

k_direct_disp = 2*pi*f_Hz./c_direct_disp;
S_k_particular_disp = 1./(2*pi*k_direct_disp).*cg_direct_disp.*S_f_particular;

transition_wavenumber = 5;

inds_direct = k_rad_m > transition_wavenumber;
inds_disp = k_direct_disp < transition_wavenumber;

k_low = 0.1;
k_high = 200;

k_interp = (k_low:k_low:k_high)';

k_rad_m_combined = [k_direct_disp(inds_disp); k_rad_m(inds_direct)];
S_k_theta_combined = [S_k_particular_disp(inds_disp,:); squeeze(S_k_theta(inds_direct,:,example_run_ind))];

S_k_theta_interp = interp1(k_rad_m_combined,S_k_theta_combined,k_interp,'pchip');

s = load('data/global_figure_settings.mat');
klow = s.k_low;

S_k_theta_interp(k_interp<klow,:) = NaN;

camera_choice = 2;

c_phase = c_br/(1/0.8);

inds_retain = c_phase > 0.23 & c_phase < sqrt(g*water_depth_m);

c_phase = c_phase(inds_retain);

Lambda_C_theta_particular = squeeze(Lambda_c_theta_02(inds_retain,:,example_run_ind))*180/pi;


c_phase_interp = logspace(log10(0.25),log10(12),100);

Lambda_C_theta_particular = interp1(c_phase,Lambda_C_theta_particular,c_phase_interp,'pchip');
c_phase = c_phase_interp;

k_disp = k_of_c(c_phase);

Lambda_C_theta_particular = smoothdata2(Lambda_C_theta_particular,'movmedian',{5,15});

dc_vec = diff(c_phase);
dk_vec = diff(k_disp);

dc_dk = dc_vec(:)./dk_vec(:);
dc_dk = [dc_dk(1); dc_dk];

Lambda_k_theta_particular = -dc_dk.*c_phase(:)./k_disp(:).*Lambda_C_theta_particular;

% S_ds, using breaking strength 'b' from Zappa et al. [2016]
% A = 3.482e-3;
% B = -4.691e-5;
% Updated values of A and B for this field campaign from Hogan et al. [2025]
A = 2.027e-3;
B = -2.166e-5;
b = A+B*(c_E/EC_ustar_m_s_particular);
S_ds_particular = b/g^2*c_phase(:).^5.*Lambda_k_theta_particular;

bigtheta_deg = [theta_br-360 theta_br theta_br+360] + IPX_offsets_deg(camera_choice) - wdir_deg;
big_S_ds = [S_ds_particular S_ds_particular S_ds_particular];

bigtheta_rad = [theta_rad-2*pi; theta_rad; theta_rad+2*pi]' - wdir_deg*pi/180;
big_S_k_particular = [S_k_theta_interp S_k_theta_interp S_k_theta_interp];
big_B_k_particular =  k_interp.^2.*big_S_k_particular;

big_S_k_particular_disp_interp = interp1(bigtheta_rad*180/pi,interp1(k_interp,big_S_k_particular,k_disp)',bigtheta_deg)';

beta_Plant = 0.04*(EC_ustar_m_s_particular./c_phase(:)).^2.*(c_phase(:).*k_disp(:)).*cosd(bigtheta_deg(:)');
S_in_Plant = beta_Plant.*big_S_k_particular_disp_interp/(rho_w*g);

big_S_in_Plant_downwind = S_in_Plant.*cosd(bigtheta_deg).*k_disp;
big_S_in_Plant_downwind(isnan(big_S_in_Plant_downwind)) = 0;
big_S_in_Plant_downwind(:,abs(bigtheta_deg)>90) = 0;
big_S_in_Plant_downwind(k_disp>70,:) = 0;

inds_consider = k_disp > 2e-1;

S_ds_levels = -9:0.5:-6;

big_S_ds_contours = smoothdata2(k_disp.^4.*big_S_ds,'movmean',{3,11});

big_S_ds_for_summing = big_S_ds;
big_S_ds_for_summing(isnan(big_S_ds)) = 0;
big_S_ds_for_summing = smoothdata2(big_S_ds_for_summing,'movmean',{3,11},'omitnan');

contour_color = flipud(reds(length(S_ds_levels)));

text_x = 0.05;
text_y = 0.95;
labels = {'(a)','(b)','(c)','(d)','(e)','(f)'};

YTicks = linspace(-5e-8,5e-8,5);

figure(fignum);clf
tlayout = tiledlayout(3,1);
ax_struc = struct();

nexttile(1)
hold on
ax_struc(1).ax = gca;
f = fill([-180 180 180 -180],[1e-1 1e-1 1e2 1e2],'k');
f.FaceAlpha = 0.5;
pcolor(bigtheta_rad*180/pi,k_interp,log10(big_B_k_particular));colormap(viridis)
for n = 1:length(S_ds_levels)
    [contour_c,contour_h] = contour(bigtheta_deg + 0.5,k_disp,log10(big_S_ds_contours),[S_ds_levels(n) 10],'-','Color',contour_color(n,:),'linewidth',2);
    clabel(contour_c,contour_h,'FontSize',fsize*0.6,'Color',contour_color(n,:),'FontName','Liberation Serif')
end
plot([0 0],[1e-10 1e10],'--','Color',0.5*[1 1 1],'linewidth',2)
hold off
box on
shading('flat')
cbar = colorbar;
cbar.Location = 'northoutside';
set(get(cbar,'Label'),'String','$\mathrm{log_{10}\{B(k,\theta)\}}$','Interpreter','LaTeX')
text(-175,0.65,'$\mathrm{log_{10}\{k^4S_{ds}(k,\theta)\}}$','Interpreter','LaTeX','FontSize',fsize,'Color','w','HorizontalAlignment','left')
ax_struc(1).ax.YScale = 'log';
shading('flat')
clim([-4 -2.5])
xlim([-180 180])
ylim([5e-1 1e2])
ylabel('$\mathrm{k\ [rad\ m^{-1}]}$','Interpreter','LaTeX')

S_in_Plant = -trapz(k_disp(inds_consider),big_S_in_Plant_downwind(inds_consider,:)./c_phase(inds_consider)');

nexttile(2)
hold on
plot([-180 180],[0 0],'--','Color',0.5*[1 1 1],'linewidth',2)
plot([0 0],[-180 180],'--','Color',0.5*[1 1 1],'linewidth',2)
% plot(bigtheta_deg,1.5*S_in_Plant,':','Color',cerulean,'linewidth',2)
% plot(bigtheta_deg,0.5*S_in_Plant,':','Color',cerulean,'linewidth',2)
f_in = fill([bigtheta_deg fliplr(bigtheta_deg)],[1.5*S_in_Plant fliplr(0.5*S_in_Plant)],cerulean);
h_in = plot(bigtheta_deg,S_in_Plant,'-','Color',cerulean,'linewidth',2);
f_in.FaceAlpha = 0.25;
f_in.LineStyle = 'none';
h_ds = plot(bigtheta_deg,trapz(k_disp(inds_consider),big_S_ds_for_summing(inds_consider,:)./c_phase(inds_consider)'),'-','Color',crimson,'linewidth',2);
hold off
box on
xlim([-180 180])
ylim([YTicks(1) YTicks(end)])
ax_struc(2).ax=gca;
ax_struc(2).ax.YTick = YTicks;
ylabel('$\mathrm{\tau(\theta)\ [N\ m^{-2} rad^{-1}]}$','Interpreter','LaTeX')
H = [h_in h_ds];
L = {'$\mathrm{\tau_{w}}$','$\mathrm{\tau_{br}}$'};
legend(H,L,'Location','southeast','Interpreter','LaTeX')

residual_upper = -1.5*trapz(k_disp(inds_consider),big_S_in_Plant_downwind(inds_consider,:)./c_phase(inds_consider)')+trapz(k_disp(inds_consider),big_S_ds_for_summing(inds_consider,:)./c_phase(inds_consider)');
residual_lower = -0.5*trapz(k_disp(inds_consider),big_S_in_Plant_downwind(inds_consider,:)./c_phase(inds_consider)')+trapz(k_disp(inds_consider),big_S_ds_for_summing(inds_consider,:)./c_phase(inds_consider)');

nexttile(3)
hold on
plot([-180 180],[0 0],'--','Color',0.5*[1 1 1],'linewidth',2)
plot([0 0],[-180 180],'--','Color',0.5*[1 1 1],'linewidth',2)
% plot(bigtheta_deg,residual_upper,'k:','linewidth',2)
% plot(bigtheta_deg,residual_lower,'k:','linewidth',2)
f_residual = fill([bigtheta_deg fliplr(bigtheta_deg)],[residual_upper fliplr(residual_lower)],'k');
f_residual.FaceAlpha = 0.25;
f_residual.LineStyle = 'none';
plot(bigtheta_deg,-trapz(k_disp(inds_consider),big_S_in_Plant_downwind(inds_consider,:)./c_phase(inds_consider)')+trapz(k_disp(inds_consider),big_S_ds_for_summing(inds_consider,:)./c_phase(inds_consider)'),'k-','linewidth',2)
hold off
xlim([-180 180])
ylim([YTicks(1) YTicks(end)])
ax_struc(3).ax=gca;
ax_struc(3).ax.YTick = YTicks;
ax_struc(3).ax.YColor = 'k';
ylabel('$\mathrm{\tau_{w}(\theta)-\tau_{br}(\theta)\ [N\ m^{-2}\ rad^{-1}]}$','Interpreter','LaTeX')
box on

for n = 1:3
    nexttile(n)
    xlabel('$\mathrm{\theta\ [^\circ]}$','Interpreter','LaTeX')
    ax_struc(n).ax.XTick = dir_ticks;
    if n == 1
        labelcolor = 'w';
    else
        labelcolor = 'k';
    end
    text(text_x,text_y,labels{n},'Color',labelcolor,'Units','normalized','HorizontalAlignment','center','FontSize',fsize)
end

tile_cleaner(ax_struc,tlayout)
