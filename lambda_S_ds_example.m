%
function lambda_S_ds_example(fignum,fsize)

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

S_f_theta = ncread(short_wave_spectra_nc_name,'S_f_theta');
f_Hz = ncread(short_wave_spectra_nc_name,'f_Hz');
theta_rad = ncread(short_wave_spectra_nc_name,'theta_rad');

F_f_d = double(ncread(long_wave_spectra_nc_name,'F_f_d'));
frequency = double(ncread(long_wave_spectra_nc_name,'frequency'));
direction = double(ncread(long_wave_spectra_nc_name,'direction'));

S_f_particular = squeeze(S_f_theta(:,:,example_run_ind));

Ff = trapz(direction,squeeze(F_f_d(example_run_ind,:,:)));
Sf = trapz(theta_rad,S_f_particular')';

f_combined = [frequency(frequency>frequency(1) & frequency <= 0.35); f_Hz(f_Hz>0.4 & f_Hz <=5)];
Ff_combined = [Ff(frequency>frequency(1) & frequency <= 0.35)'; ((2*pi*f_Hz(f_Hz>0.4 & f_Hz <=5)).^2/g).^-2.*Sf(f_Hz>0.4 & f_Hz <=5)];
T_E = trapz(f_combined,f_combined.^-1.*Ff_combined)/trapz(f_combined,Ff_combined);
f_E = T_E^-1;

[c_E,~] = lindisp_with_current(2*pi*f_E,water_depth_m,0);

camera_choice = 2;

c_phase = c_br/(1/0.8);

inds_retain = c_phase > 0.23 & c_phase < sqrt(g*water_depth_m);

c_phase = c_phase(inds_retain);

k_disp = k_of_c(c_phase);

Lambda_C_theta_particular = squeeze(Lambda_c_theta_02(inds_retain,:,example_run_ind))*180/pi;

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

big_Lambda_c_theta = [Lambda_C_theta_particular Lambda_C_theta_particular Lambda_C_theta_particular];
big_Lambda_c_theta(big_Lambda_c_theta==0) = NaN;

bigtheta_deg = [theta_br-360 theta_br theta_br+360] + IPX_offsets_deg(camera_choice) - wdir_deg;
big_S_ds = [S_ds_particular S_ds_particular S_ds_particular];

inds_downwind = abs(bigtheta_deg) < 180;

Sm = mean(sind(bigtheta_deg(inds_downwind)).*big_S_ds(:,inds_downwind).*k_disp./c_phase(:),'all','omitnan');
Cm = mean(cosd(bigtheta_deg(inds_downwind)).*big_S_ds(:,inds_downwind).*k_disp./c_phase(:),'all','omitnan');
theta_mean_br = 180/pi*atan2(Sm,Cm);

k_mean_br = sum(k_disp.*big_S_ds(:,inds_downwind).*k_disp./c_phase(:),'all','omitnan')/sum(big_S_ds(:,inds_downwind).*k_disp./c_phase(:),'all','omitnan');

text_x = 0.04;
text_y = 1.04;
labels = {'(a)','(b)'};

fA = 0.8;

Lambda_contour_levels = -5.25:0.25:-4;
S_ds_contour_levels = -14:1:-8;

figure(fignum);clf
tlayout = tiledlayout(2,1);
ax_struc = struct();

nexttile(1)
hold on
contourf(bigtheta_deg+0.5,c_phase,log10(big_Lambda_c_theta),Lambda_contour_levels)
plot([0 0],[1e-10 1e10],'--','Color',0.5*[1 1 1],'linewidth',2)
f = fill([-175 -175 175 175],[4.5 5.9 5.9 4.5],'w');
hold off
box on
f.FaceAlpha = fA;
f.LineStyle = 'none';
cbar = colorbar;
set(get(cbar,'Label'),'String','$\mathrm{log_{10}\{\Lambda(c,\theta)\}\ [s^2 m^{-1} rad^{-1}]}$','Interpreter','LaTeX')
cbar.Location = 'north';
xlim([-180 180])
ylim([0 6])
clim([Lambda_contour_levels(2) Lambda_contour_levels(end)]+[-1 1]*median(diff(Lambda_contour_levels)))
ax_struc(1).ax = gca;
xlabel('$\mathrm{\theta\ [^\circ]}$','Interpreter','LaTeX')
ylabel('$\mathrm{c\ [m\ s^{-1}]}$','Interpreter','LaTeX')
cbar.Ticks = Lambda_contour_levels(2:end);

nexttile(2)
hold on
plot(theta_mean_br*[1 1],[1e-10 k_mean_br],':','Color',0.5*[1 1 1],'linewidth',2.5)
contourf(bigtheta_deg + 0.5,k_disp,log10(big_S_ds),S_ds_contour_levels)
plot([0 0],[1e-10 1e10],'--','Color',0.5*[1 1 1],'linewidth',2)
f = fill([-175 -175 175 175],[22 90 90 22],'w');
hold off
box on
f.FaceAlpha = fA;
f.LineStyle = 'none';
xlim([-180 180])
ylim([1e-1 1e2])
clim([S_ds_contour_levels(2) S_ds_contour_levels(end)]+[-1 1]*median(diff(S_ds_contour_levels)))
xlabel('$\mathrm{\theta\ [^\circ]}$','Interpreter','LaTeX')
ylabel('$\mathrm{k\ [rad\ m^{-1}]}$','Interpreter','LaTeX')
text(80,3e-1,['$\mathrm{\bar{\theta}_{br}}$ = ' sprintf('%0.2f',theta_mean_br) '$^\circ$'],'HorizontalAlignment','Center','FontSize',fsize,'Interpreter','LaTeX')
cbar = colorbar;
set(get(cbar,'Label'),'String','$\mathrm{log_{10}\{S_{ds}(k,\theta)\}\ [m^4 rad^{-1} s^{-1}]}$','Interpreter','LaTeX')
cbar.Ticks = S_ds_contour_levels(2:end);
cbar.Location = 'north';
ax_struc(2).ax = gca;

for i = 1:2
    nexttile(i)
    ax_struc(i).ax.XTick = dir_ticks;
    text(text_x,text_y,labels{i},'FontSize',fsize,'Units','Normalized','HorizontalAlignment','Center')
end

ax_struc(2).ax.YScale = 'log';

tile_cleaner(ax_struc,tlayout)

colormap(ax_struc(1).ax,deep(length(Lambda_contour_levels)))
colormap(ax_struc(2).ax,cividis(length(S_ds_contour_levels)))

tlayout.TileSpacing = 'tight';