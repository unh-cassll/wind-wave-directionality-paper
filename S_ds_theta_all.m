%
function S_ds_theta_all(fignum,fsize)

supporting_data_nc_name = 'data/ASIT2019_supporting_data_compilation.nc';

load('data/ASIT2019_Lambda_c_theta_all.mat')

g = 9.806;
sigma = 0.072;
rho_w = 1020;

k = 0.01:0.01:300;
water_depth_m = 15;

omega_disp = sqrt((g*k+sigma/rho_w*k.^3).*tanh(k*water_depth_m));

c_disp = omega_disp./k;

k_of_c = fit(c_disp(:),k(:),'spline');

dir_ticks = 180*(-1:0.25:1);

DN_start = datenum(2019,10,16,0,0,0);
time = double(ncread(supporting_data_nc_name,'time'));
DN = DN_start + time/60/24;
DTime = datetime(datevec(DN));

wdir_deg_all = mod(ncread(supporting_data_nc_name,'COARE_winddir')+180,360);
EC_ustar_m_s = ncread(supporting_data_nc_name,'EC_ustar');
EC_U_m_s = ncread(supporting_data_nc_name,'EC_U');
EC_z_m_above_water = ncread(supporting_data_nc_name,'EC_z');

EC_tau_downwind_N_m2 = ncread(supporting_data_nc_name,'EC_tau_downwind');
EC_tau_crosswind_N_m2 = ncread(supporting_data_nc_name,'EC_tau_crosswind');

EC_stress_angle_deg = 180/pi*atan2(EC_tau_crosswind_N_m2,EC_tau_downwind_N_m2);

kappa = 0.4;

EC_z0_m = EC_z_m_above_water.*exp(-kappa*EC_U_m_s./EC_ustar_m_s);
EC_U10_m_s = EC_ustar_m_s/kappa.*log(10./EC_z0_m);

f_Hz = ncread(supporting_data_nc_name,'frequency_omni');
F_f_m2_Hz = ncread(supporting_data_nc_name,'Omni_F_f');

inds_keep = f_Hz > 5e-2 & f_Hz < 1;

F_f_m2_Hz(isnan(F_f_m2_Hz)) = 0;
f_p = trapz(f_Hz(inds_keep),f_Hz(inds_keep).*F_f_m2_Hz(inds_keep,:).^4)'./trapz(f_Hz(inds_keep),F_f_m2_Hz(inds_keep,:).^4)';
f_p = f_p(:);

f_Hz_ADCP = ncread(supporting_data_nc_name,'frequency_dir')';
theta_rad_ADCP = ncread(supporting_data_nc_name,'direction')';
F_f_theta_m2_Hz_rad_ADCP = ncread(supporting_data_nc_name,'Dir_F_f_theta');

Sm = squeeze(mean(sin(theta_rad_ADCP(:)').*F_f_theta_m2_Hz_rad_ADCP,[1 2],'omitnan'));
Cm = squeeze(mean(cos(theta_rad_ADCP(:)').*F_f_theta_m2_Hz_rad_ADCP,[1 2],'omitnan'));

Vm = atan2(Sm, Cm);

D_m = 180/pi*Vm;

[c_p,~] = lindisp_with_current(2*pi*f_p,water_depth_m,0);

c_p = interp1(DTime,medfilt1(c_p,5,'includenan'),DTime,'pchip');
EC_ustar_m_s_filt = interp1(DTime,medfilt1(EC_ustar_m_s,5,'includenan'),DTime,'pchip');
EC_U10_m_s_filt = interp1(DTime,medfilt1(EC_U10_m_s,5,'includenan'),DTime,'pchip');

EC_ustar_m_s_filt(EC_ustar_m_s_filt<=0) = NaN;
EC_U10_m_s_filt(EC_U10_m_s_filt<=0) = NaN;

camera_choice = 2;

c_phase = c_br/(1/0.8);

inds_retain = c_phase > 0.23 & c_phase < sqrt(g*water_depth_m);

c_phase = c_phase(inds_retain);

k_disp = k_of_c(c_phase);

dc_vec = diff(c_phase)';
dc_vec = [dc_vec(1); dc_vec];
dk_vec = diff(k_disp);
dk_vec = [dk_vec(1); dk_vec];

dc_dk = dc_vec(:)./dk_vec(:);

dtheta_deg = median(diff(theta_br));

N = length(DTime_02);

k_ds = NaN*ones(N,1);
theta_p_ds = k_ds;
tau_theta = NaN*ones(N,length(theta_br));

ustar_m_s = k_ds;
stress_ang_deg = k_ds;
U10_m_s = k_ds;
c_p_m_s = k_ds;
wdir_deg = k_ds;
D_m_deg = k_ds;
Energy_flux_ds_int = k_ds;
Momentum_flux_ds_int = k_ds;

theta_downwind = theta_br - mean(theta_br);

quantiles = [10 25 50 75 90];
k_ds_quantiles = NaN*ones(N,length(quantiles));

for IPX_ind = 1:N

    DN_camera = datenum(DTime_02(IPX_ind));

    DN_diff = abs(DN-DN_camera);

    waves_and_wind_ind = find(DN_diff==min(DN_diff),1,'first');

    if min(DN_diff)*86400 < 3*3600

        wdir_deg_particular = wdir_deg_all(waves_and_wind_ind);
        EC_ustar_m_s_particular = EC_ustar_m_s_filt(waves_and_wind_ind);
        EC_U10_m_s_particular = EC_U10_m_s_filt(waves_and_wind_ind);
        EC_stress_angle_deg_particular = EC_stress_angle_deg(waves_and_wind_ind);

        c_p_particular = c_p(waves_and_wind_ind);
        D_m_particular = D_m(waves_and_wind_ind);

        Lambda_C_theta_particular = squeeze(Lambda_c_theta_02(inds_retain,:,IPX_ind))*180/pi;
        Lambda_C_theta_plot = smoothdata2(Lambda_C_theta_particular,'movmedian',{5,15});      
        sum_Lambda_c_theta = sum(Lambda_C_theta_particular.*c_phase(:).*dc_vec(:),'all','omitnan');

        Lambda_C_theta_particular = smoothdata2(Lambda_C_theta_particular,'movmedian',{5,5});        

        Lambda_C_theta_particular = Lambda_C_theta_particular*sum_Lambda_c_theta/sum(Lambda_C_theta_particular.*c_phase(:).*dc_vec(:),'all','omitnan');
        Lambda_C_theta_plot = Lambda_C_theta_plot*sum_Lambda_c_theta/sum(Lambda_C_theta_plot.*c_phase(:).*dc_vec(:),'all','omitnan');

        if sum_Lambda_c_theta > 1e-3 && EC_ustar_m_s_particular > 0

            Lambda_k_theta_plot = -dc_dk.*c_phase(:)./k_disp(:).*Lambda_C_theta_plot;

            % S_ds, using breaking strength 'b' from Zappa et al. [2016]
            % A = 3.482e-3;
            % B = -4.691e-5;
            % Updated values of A and B for this field campaign from Hogan et al. [2025]
            A = 2.027e-3;
            B = -2.166e-5;
            b = A+B*(c_p_particular/EC_ustar_m_s_particular);

            Momentum_flux_ds_int(IPX_ind) = rho_w/g*b*sum(c_phase(:).^5.*dc_vec(:).*Lambda_C_theta_particular,'all','omitnan')*dtheta_deg;
            Energy_flux_ds_int(IPX_ind) = rho_w/g*b*sum(c_phase(:).^6.*dc_vec(:).*Lambda_C_theta_particular,'all','omitnan')*dtheta_deg;

            S_ds_particular = b/g^2*c_phase(:).^5.*Lambda_k_theta_plot;
            S_ds_particular(S_ds_particular==0) = NaN;

            bigtheta_deg = [theta_br-720 theta_br-360 theta_br theta_br+360 theta_br+720] + IPX_offsets_deg(camera_choice);
            big_S_ds = [S_ds_particular S_ds_particular S_ds_particular S_ds_particular S_ds_particular];

            tau_theta_bit = sum(big_S_ds./c_phase(:).*k_disp.*dk_vec,1,'omitnan');
            S_ds_k_bit = sum(big_S_ds(:,abs(bigtheta_deg)<180).*k_disp,2,'omitnan');

            S_ds_k_bit_cumu_norm = 1-cumtrapz(k_disp,S_ds_k_bit)/trapz(k_disp,S_ds_k_bit);

            [S_ds_k_bit_cumu_norm_unique,ia,~] = unique(S_ds_k_bit_cumu_norm);
            k_unique = k_disp(ia);

            S_ds_frac_fit = fit(S_ds_k_bit_cumu_norm_unique(1:end-1),k_unique(1:end-1),'pchip');

            k_ds_quantiles(IPX_ind,:) = S_ds_frac_fit(quantiles/100);

            inds_max = find(abs(tau_theta_bit)==max(abs(tau_theta_bit)));

            inds_downbreak = (-179:1:180) + inds_max(3);

            theta_p_ds(IPX_ind) = bigtheta_deg(inds_max(3));

            big_S_ds_for_summing = smoothdata2(big_S_ds,'movmean',{3,3},'omitnan');

            tau_theta_bit = sum(big_S_ds_for_summing(:,inds_downbreak).*k_disp.*dk_vec,1,'omitnan');

            k_ds(IPX_ind) = sum(k_disp.^4.*big_S_ds(:,inds_downbreak),'all','omitnan')/sum(k_disp.^3.*big_S_ds(:,inds_downbreak),'all','omitnan');

            tau_theta(IPX_ind,:) = tau_theta_bit;

            ustar_m_s(IPX_ind) = EC_ustar_m_s_particular;
            U10_m_s(IPX_ind) = EC_U10_m_s_particular;
            c_p_m_s(IPX_ind) = c_p_particular;
            wdir_deg(IPX_ind) = wdir_deg_particular;
            D_m_deg(IPX_ind) = D_m_particular;
            stress_ang_deg(IPX_ind) = EC_stress_angle_deg_particular;

        end

    end

end

Sm = mean(sind(theta_br).*tau_theta,2,'omitnan');
Cm = mean(cosd(theta_br).*tau_theta,2,'omitnan');

Vm = atan2(Sm, Cm);
Vm(Vm<0) = Vm(Vm<0) + 2*pi;

theta_br = mod(theta_p_ds + 180/pi*Vm+360,360);
theta_br(theta_br-wdir_deg>180) = theta_br(theta_br-wdir_deg>180) - 360;
theta_br(theta_br-wdir_deg>120) = theta_br(theta_br-wdir_deg>120) - 180;
theta_br(wdir_deg-theta_br>120) = theta_br(wdir_deg-theta_br>120) + 180;

save('data/integrated_wave_breaking_quantities.mat','theta_br','Energy_flux_ds_int','Momentum_flux_ds_int','k_ds','ustar_m_s','U10_m_s','c_p_m_s','wdir_deg','D_m_deg','stress_ang_deg')

text_x = 0.05;
text_y = 0.95;
labels = {'(a)','(b)'};

s = load('data/global_figure_settings.mat');
wave_age_lims = s.wave_age_lims;
d_wave_age = 10;
wave_age_centers = wave_age_lims(1)+d_wave_age/2:d_wave_age:wave_age_lims(2)-d_wave_age/2;

wave_age = c_p_m_s./ustar_m_s;

U10_centers = 2:2:12;
dU10 = median(diff(U10_centers));

cmap_binned_U10 = cividis(length(U10_centers));
cmap_binned_wave_age = flipud(magma(length(wave_age_centers)));

tau_br_theta_binned_by_wave_age = NaN*ones(length(theta_downwind),length(wave_age_centers));
tau_theta_binned_by_U10 = NaN*tau_br_theta_binned_by_wave_age;

k_ds_quantiles_binned_by_wave_age = NaN*ones(length(quantiles),length(wave_age_centers));
ustar_binned = NaN*wave_age_centers;

for n = 1:length(wave_age_centers)

    wave_age_low = wave_age_centers(n) - d_wave_age/2;
    wave_age_high = wave_age_centers(n) + d_wave_age/2;
    inds_consider = wave_age >= wave_age_low & wave_age < wave_age_high;

    tau_br_theta_binned_by_wave_age(:,n) = median(tau_theta(inds_consider,:),1,'omitnan');
    k_ds_quantiles_binned_by_wave_age(:,n) = median(k_ds_quantiles(inds_consider,:),1,'omitnan');

    ustar_binned(n) = mean(ustar_m_s(inds_consider),'all','omitnan');

end

for n = 1:length(U10_centers)

    U_low = U10_centers(n) - dU10/2;
    U_high = U10_centers(n) + dU10/2;
    inds_consider = U10_m_s >= U_low & U10_m_s < U_high;

    tau_theta_binned_by_U10(:,n) = median(tau_theta(inds_consider,:),1,'omitnan');

end

option = 'waveage';
% option = 'U10';

switch option

    case 'waveage'

        tau_br_theta_binned = tau_br_theta_binned_by_wave_age;
        cmap_binned = cmap_binned_wave_age;
        clims = [wave_age_centers(1) wave_age_centers(end)] + d_wave_age/2*[-1 1];

    case 'U10'

        tau_br_theta_binned = tau_theta_binned_by_U10;
        cmap_binned = cmap_binned_U10;
        clims = [U10_centers(1) U10_centers(end)] + [-1 1]*dU10/2;

end

d_line = 0.012;
bottom_plot_vertical_levels = linspace(-1.2+d_line*2,-1-d_line*2,length(wave_age_centers));

momentum_flux_ds_theta_binned = rho_w*g*tau_br_theta_binned;

figure(fignum);clf
tlayout = tiledlayout(2,1);
ax_struc = struct();

nexttile(1)
hold on
plot([0 0],[-10 10],'--','Color',0.5*[1 1 1],'linewidth',2)
plot(theta_downwind,momentum_flux_ds_theta_binned,'k-','linewidth',3)
for i = 1:length(wave_age_centers)
    n = length(wave_age_centers)-i+1;
    plot(theta_downwind,momentum_flux_ds_theta_binned(:,i),'Color',cmap_binned(i,:),'linewidth',2)
end
hold off
cbar = colorbar;
cbar.Layout.Tile = 'north';
colormap(cmap_binned)
switch option
    case 'U10'
        set(get(cbar,'Label'),'String','$\mathrm{U_{10}\ [m\ s^{-1}]}$','Interpreter','LaTeX')
        clim(clims)
    case 'waveage'
        set(get(cbar,'Label'),'String','$\mathrm{c_p/u_*}$','Interpreter','LaTeX')
        cbar.Ticks = clims(1):d_wave_age:clims(end);
        clim(clims)
end
box on
ylim(-[1.0e-2 0])
xlim([-1 1]*180)
xlabel('$\mathrm{\theta-\theta_{br}\ [rad]}$','Interpreter','LaTeX')
ylabel('$\mathrm{\tau_{br}(\theta)\ [N\ m^{-2}rad^{-1}]}$','Interpreter','LaTeX')
ax_struc(1).ax = gca;

tau_br_norm = tau_br_theta_binned./min(tau_br_theta_binned);

tau_br_cumu_int_norm = cumtrapz(theta_downwind,tau_br_theta_binned)./trapz(theta_downwind,tau_br_theta_binned);

tau_br_quantile_angles = NaN*ones(length(quantiles),length(wave_age_centers));

nexttile(2)
hold on
plot([0 0],[-10 10],'--','Color',0.5*[1 1 1],'linewidth',2)
plot(theta_downwind,-tau_br_norm,'k-','linewidth',3)
for i = 1:length(wave_age_centers)
    plot(theta_downwind,-tau_br_norm(:,i),'Color',cmap_binned(i,:),'linewidth',2)

    i_s = find(tau_br_cumu_int_norm(:,i)>0.1,1,'first');
    i_f = find(tau_br_cumu_int_norm(:,i)>0.9,1,'first');

    x = [[theta_downwind(i_s) theta_downwind(i_f)] fliplr([theta_downwind(i_s) theta_downwind(i_f)])];
    y = bottom_plot_vertical_levels(i)+[-1 -1 1 1]*d_line;
    f = fill(x,y,cmap_binned(i,:));
    f.FaceAlpha = 0.25;

    tau_br_quantile_angles(1,i) = theta_downwind(i_s);
    tau_br_quantile_angles(5,i) = theta_downwind(i_f);

    i_m = find(tau_br_cumu_int_norm(:,i)>0.5,1,'first');

    tau_br_quantile_angles(3,i) = theta_downwind(i_m);

        i_s = find(tau_br_cumu_int_norm(:,i)>0.25,1,'first');
    i_f = find(tau_br_cumu_int_norm(:,i)>0.75,1,'first');

    x = [[theta_downwind(i_s) theta_downwind(i_f)] fliplr([theta_downwind(i_s) theta_downwind(i_f)])];
    y = bottom_plot_vertical_levels(i)+[-1 -1 1 1]*d_line;
    f = fill(x,y,cmap_binned(i,:));
    f.FaceAlpha = 0.5;

    tau_br_quantile_angles(2,i) = theta_downwind(i_s);
    tau_br_quantile_angles(4,i) = theta_downwind(i_f);

end
hold off
box on
ylim([-1.2 0])
xlim([-1 1]*180)
xlabel('$\mathrm{\theta-\theta_{br}\ [rad]}$','Interpreter','LaTeX')
ylabel('$\tau_{\mathrm{br}}(\theta)$, norm.','Interpreter','latex')
ax_struc(2).ax = gca;

for i = 1:2
    nexttile(i)
    ax_struc(i).ax.XTick = dir_ticks;
    text(text_x,text_y,labels{i},'FontSize',fsize,'Units','Normalized','HorizontalAlignment','Center')
end

tile_cleaner(ax_struc,tlayout)
tlayout.TileSpacing = 'tight';

save('data/wave_age_binned_breaking_quantities.mat','ustar_binned','tau_br_quantile_angles','k_ds_quantiles_binned_by_wave_age')
