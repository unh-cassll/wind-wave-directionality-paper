%
function wave_age_binned_directional_spreading(fignum,fsize)

short_wave_nc_name = 'data/ASIT2019_wave_spectra_stats_timeseries_empirical_gain.nc';
supporting_nc_name = 'data/ASIT2019_supporting_environmental_observations.nc';

load('data/delta_with_wind_speed_data.mat')
load('data/directional_spreading_quantification_data.mat')

wave_age_binned_breaking_quantities = load('data/wave_age_binned_breaking_quantities.mat');

f_Hz_Pyxis = double(ncread(short_wave_nc_name,'f_Hz'));
k_rad_m_Pyxis = double(ncread(short_wave_nc_name,'k_rad_m'));
EC_ustar_m_s = ncread(supporting_nc_name,'EC_ustar_m_s');
load('data/ASIT2019_peak_wave_phase_speed.mat')

load('data/LenainMelville2017_kn.mat');
load('data/LenainMelville2017_r.mat');

freq_spect_range_limits = load('data/frequency_spect_range_limits.mat');
wavenumber_spect_range_limits = load('data/wavenumber_spect_range_limits.mat');

U_sfc_mag_m_s = ncread(supporting_nc_name,'U_sfc_mag_m_s');

g = 9.81;
water_depth_m = 15;

f_p = f_p(:);
f_eq_start = freq_spect_range_limits.f_eq_start(:);
f_eq_end = freq_spect_range_limits.f_eq_end(:);
f_sat_end = freq_spect_range_limits.f_sat_end(:);
k_sat_end = wavenumber_spect_range_limits.k_sat_end(:);

k_p_disp = NaN*f_p;
k_eq_start_disp = k_p_disp;
k_eq_end_disp = k_p_disp;
k_sat_end_disp = k_p_disp;

for n = 1:length(U_sfc_mag_m_s)

    try

        [c,~] = lindisp_with_current(2*pi*f_p(n),water_depth_m,U_sfc_mag_m_s(n));
        k_p_disp(n) = 2*pi*f_p(n)./c;

        [c,~] = lindisp_with_current(2*pi*f_eq_start(n),water_depth_m,U_sfc_mag_m_s(n));
        k_eq_start_disp(n) = 2*pi*f_eq_start(n)./c;

        [c,~] = lindisp_with_current(2*pi*f_eq_end(n),water_depth_m,U_sfc_mag_m_s(n));
        k_eq_end_disp(n) = 2*pi*f_eq_end(n)./c;

        [c,~] = lindisp_with_current(2*pi*f_sat_end(n),water_depth_m,U_sfc_mag_m_s(n));
        k_sat_end_disp(n) = 2*pi*f_sat_end(n)./c;

    end

end

wave_age_lims = [10 70];
d_wave_age = 10;
wave_age_centers = wave_age_lims(1):d_wave_age:wave_age_lims(2);

clims = [wave_age_centers(1) wave_age_centers(end)] + d_wave_age/2*[-1 1];

wave_age = C_p./EC_ustar_m_s;

text_x = 0.05;
text_y = 0.95;

labels = {'(a)','(b)','(c)','(d)'};

D_k_limits_block = NaN*ones(size(D_spread_holder_struc(1).D_k_limits,1),size(D_spread_holder_struc(1).D_k_limits,2),length(D_spread_holder_struc));
D_f_limits_block = NaN*ones(size(D_spread_holder_struc(1).D_f_limits,1),size(D_spread_holder_struc(1).D_f_limits,2),length(D_spread_holder_struc));

k_disp_block = NaN*ones(length(f_Hz_Pyxis),190);

for n = 1:length(D_spread_holder_struc)
    D_k_limits_block(:,:,n) = D_spread_holder_struc(n).D_k_limits;
    D_f_limits_block(:,:,n) = D_spread_holder_struc(n).D_f_limits;

    try

        [c,~] = lindisp_with_current(2*pi*f_Hz_Pyxis,water_depth_m,U_sfc_mag_m_s(n));
        k_disp_block(:,n) = 2*pi*f_Hz_Pyxis./c;

    end

end

D_k_50th_halfwidth = 180/pi*squeeze(D_k_limits_block(:,2,:) - D_k_limits_block(:,1,:))/2;
D_f_50th_halfwidth = 180/pi*squeeze(D_f_limits_block(:,2,:) - D_f_limits_block(:,1,:))/2;

s = load('data/global_figure_settings.mat');
klow = s.k_low;

flow = sqrt(9.81*klow)/(2*pi);

% higher limits than for directional plots
k_high = 550;
f_high = 15;

halfwidth_lims = [0 180];
k_hat_lims = [1e-3 1e2];
k_hat_norm_lims = [5e-2 1e2];
Delta_lims = [-1 1]*1;

dir_ticks = 180*(-1:0.25:1);
Delta_ticks = -1:0.5:1;

delta_k(k_rad_m_Pyxis<klow,:) = NaN;
delta_f(f_Hz_Pyxis<flow,:) = NaN;
D_k_50th_halfwidth(k_rad_m_Pyxis<klow,:) = NaN;
D_f_50th_halfwidth(f_Hz_Pyxis<flow,:) = NaN;

delta_k(k_rad_m_Pyxis>k_high,:) = NaN;
delta_f(f_Hz_Pyxis>f_high,:) = NaN;
D_k_50th_halfwidth(k_rad_m_Pyxis>k_high,:) = NaN;
D_f_50th_halfwidth(f_Hz_Pyxis>f_high,:) = NaN;

k_hat_binned = logspace(log10(k_hat_lims(1)),log10(k_hat_lims(2)),100)';
delta_k_binned = NaN*ones(length(k_hat_binned),length(wave_age_centers));
D_k_50th_binned = delta_k_binned;

k_hat_disp_binned = k_hat_binned;
delta_f_binned = NaN*ones(length(k_hat_disp_binned),length(wave_age_centers));
D_f_50th_binned = delta_f_binned;
k_n_binned = delta_f_binned(1,:);
k_n_hat_binned = k_n_binned;
k_s_hat_binned = k_n_binned;
ind_match_n = k_n_binned;
ind_match_s = k_n_binned;

k_sg_norm_n_s_binned = NaN*[k_n_binned; k_n_binned];
k_gc_norm_n_s_binned = k_sg_norm_n_s_binned;

k_sg = 117;
k_gc = 371;

LM2017_r_binned = NaN*wave_age_centers;

for n = 1:length(wave_age_centers)

    wave_age_low = wave_age_centers(n) - d_wave_age/2;
    wave_age_high = wave_age_centers(n) + d_wave_age/2;
    inds_consider = wave_age >= wave_age_low & wave_age < wave_age_high;

    inds_LM2017 = LM2017_Cp_ustar >= wave_age_low & LM2017_Cp_ustar < wave_age_high;
    LM2017_r_binned(n) = mean(LM2017_r(inds_LM2017));

    delta_k_consider = delta_k(:,inds_consider);
    D_k_50th_halfwidth_consider = D_k_50th_halfwidth(:,inds_consider);
    ustar_consider = EC_ustar_m_s(inds_consider);

    k_hat = k_rad_m_Pyxis*median(ustar_consider,'all','omitnan')^2/g;

    delta_k_block = NaN*ones(length(k_hat_binned),length(ustar_consider));
    D_k_50th_block = delta_k_block;

    for m = 1:length(ustar_consider)
        delta_k_block(:,m) = interp1(k_hat,delta_k_consider(:,m),k_hat_binned);
        D_k_50th_block(:,m) = interp1(k_hat,D_k_50th_halfwidth_consider(:,m),k_hat_binned);
    end

    k_n_binned(n) = mean(k_eq_end_disp(inds_consider),'omitnan');
    k_n_hat_binned(n) = mean(k_eq_end_disp(inds_consider),'omitnan')*median(ustar_consider,'all','omitnan')^2/g;
    k_s_hat_binned(n) = mean(k_sat_end(inds_consider),'omitnan')*median(ustar_consider,'all','omitnan')^2/g;

    k_sg_norm_n_s_binned(1,n) = k_sg/mean(k_eq_end_disp(inds_consider),'omitnan');
    k_sg_norm_n_s_binned(2,n) = k_sg/mean(k_sat_end(inds_consider),'omitnan');

    k_gc_norm_n_s_binned(1,n) = k_gc/mean(k_eq_end_disp(inds_consider),'omitnan');
    k_gc_norm_n_s_binned(2,n) = k_gc/mean(k_sat_end(inds_consider),'omitnan');

    diff_k_n = abs(k_n_hat_binned(n)-k_hat_binned);
    ind_match_n(n) = find(diff_k_n==min(diff_k_n,[],'all','omitnan'));
    diff_k_s = abs(k_s_hat_binned(n)-k_hat_binned);
    ind_match_s(n) = find(diff_k_s==min(diff_k_s,[],'all','omitnan'));

    delta_k_binned(:,n) = mean(delta_k_block,2,'omitnan');
    D_k_50th_binned(:,n) = pi/180*mean(D_k_50th_block,2,'omitnan');

    delta_f_consider = delta_f(:,inds_consider);
    D_f_50th_halfwidth_consider = D_f_50th_halfwidth(:,inds_consider);
    ustar_consider = EC_ustar_m_s(inds_consider);

    k_disp = k_disp_block(:,n);
    k_hat_disp = k_disp*median(ustar_consider,'all','omitnan')^2/g;

    delta_f_block = NaN*ones(length(k_hat_disp_binned),length(ustar_consider));
    D_f_50th_block = delta_f_block;

    for m = 1:length(ustar_consider)
        delta_f_block(:,m) = interp1(k_hat_disp,delta_f_consider(:,m),k_hat_disp_binned);
        D_f_50th_block(:,m) = interp1(k_hat_disp,D_f_50th_halfwidth_consider(:,m),k_hat_disp_binned);
    end

    delta_f_binned(:,n) = mean(delta_f_block,2,'omitnan');
    D_f_50th_binned(:,n) = pi/180*mean(D_f_50th_block,2,'omitnan');

    ind_max_theta_halfwidth_f(n) = find(D_f_50th_binned(:,n)==nanmax(D_f_50th_binned(:,n)),1,'first');
    ind_min_Delta_f(n) = find(delta_f_binned(:,n)==nanmin(delta_f_binned(:,n)),1,'first');

end

cmap_binned = flipud(magma(length(wave_age_centers)));

alpha_vec = flipud(linspace(0.3,1,length(wave_age_centers))');

S_ds_id_range = wave_age_binned_breaking_quantities.tau_br_quantile_angles(5,:) - wave_age_binned_breaking_quantities.tau_br_quantile_angles(1,:);

k_ds_hat = wave_age_binned_breaking_quantities.k_ds_quantiles_binned_by_wave_age.*wave_age_binned_breaking_quantities.ustar_binned.^2/g;

figure(fignum);clf
tlayout = tiledlayout(2,2);
ax_struc = struct();

nexttile(1)
hold on
for i = 1:length(wave_age_centers)
    n = length(wave_age_centers)-i+1;
    plot([3*pi/4 D_f_50th_binned(ind_max_theta_halfwidth_f(n),n)]*180/pi,k_n_hat_binned(n)*[1 1],'Color',[cmap_binned(n,:) alpha_vec(n)],'linewidth',3)
    plot([D_f_50th_binned(ind_max_theta_halfwidth_f(n),n) 3*pi/4]*180/pi,0.96*k_n_hat_binned(n)*[1 1],'k','linewidth',0.5)
    plot([D_f_50th_binned(ind_max_theta_halfwidth_f(n),n) 3*pi/4]*180/pi,1.04*k_n_hat_binned(n)*[1 1],'k','linewidth',0.5)
    plot(3*pi/4*180/pi,k_n_hat_binned(n),'o','markersize',6,'markerfacecolor',cmap_binned(n,:),'markeredgecolor','k','linewidth',0.75)
    S_ds_fill = fill(S_ds_id_range(i)+[-1 -1 1 1]*2.5,[k_ds_hat(1,i) k_ds_hat(5,i) k_ds_hat(5,i) k_ds_hat(1,i)],cmap_binned(i,:));
    S_ds_fill.FaceAlpha = 0.5;
end
plot(D_f_50th_binned*180/pi,k_hat_disp_binned,'k-','linewidth',3)
for i = 1:length(wave_age_centers)
    n = length(wave_age_centers)-i+1;
    plot(D_f_50th_binned(:,n)*180/pi,k_hat_binned,'Color',cmap_binned(n,:),'linewidth',2)
end
hold off
colororder(cmap_binned)
box on
ylim(k_hat_lims)
xlim(halfwidth_lims)
xlabel('\theta_{halfwidth} [\circ]')
ylabel('$\mathrm{\hat{k}_{disp}\equiv k_{disp}u_*^2/g}$','Interpreter','LaTeX')
ax_struc(1).ax = gca;

nexttile(2)
hold on
for i = 1:length(wave_age_centers)
    n = length(wave_age_centers)-i+1;
    plot([-0.5 delta_f_binned(ind_min_Delta_f(n),n)],k_n_hat_binned(n)*[1 1],'Color',[cmap_binned(n,:) alpha_vec(n)],'linewidth',3)
    plot([delta_f_binned(ind_min_Delta_f(n),n) -0.5],0.96*k_n_hat_binned(n)*[1 1],'k','linewidth',0.5)
    plot([delta_f_binned(ind_min_Delta_f(n),n) -0.5],1.04*k_n_hat_binned(n)*[1 1],'k','linewidth',0.5)
    plot(-0.5,k_n_hat_binned(n),'o','markersize',6,'markerfacecolor',cmap_binned(n,:),'markeredgecolor','k','linewidth',0.75)
end

cmap_LM2017 = flipud(magma(255));
inds = floor((LM2017_Cp_ustar-clims(1))/(clims(2)-clims(1))*255);
inds(inds>255) = 255;
for i = 1:length(LM2017_Cp_ustar)
    LM2017_r_fill = fill([0.6 0.6 0.9 0.9],[0.95 1.05 1.05 0.95]*LM2017_r(i),cmap_LM2017(inds(i),:));
    LM2017_r_fill.FaceAlpha = 0.5;
    LM2017_r_fill.LineStyle = 'none';
end
text(0.75,6e-3,{'L&M 2017'},'HorizontalAlignment','center','FontSize',fsize,'Color',0.35*[1 1 1])
plot([0 0],k_hat_lims,'--','Color',0.5*[1 1 1],'linewidth',2)
plot(delta_f_binned,k_hat_disp_binned,'k-','linewidth',3)
for i = 1:length(wave_age_centers)
    n = length(wave_age_centers)-i+1;
    plot(delta_f_binned(:,n),k_hat_binned,'Color',cmap_binned(n,:),'linewidth',2)
end
hold off
box on
ylim(k_hat_lims)
xlim(Delta_lims)
xlabel('\Delta [rad]')
ylabel('$\mathrm{\hat{k}_{disp}\equiv\omega^2u_*^2/g^2}$','Interpreter','LaTeX')
ax_struc(2).ax = gca;
clim(clims)
colormap(cmap_binned)

nexttile(3)
hold on
for i = 1:length(wave_age_centers)
    n = length(wave_age_centers)-i+1;
    plot([D_k_50th_binned(ind_match_s(n),n) 3*pi/4]*180/pi,k_s_hat_binned(n)*[1 1],'Color',[cmap_binned(n,:) alpha_vec(n)],'linewidth',3)
    plot([D_k_50th_binned(ind_match_s(n),n) 3*pi/4]*180/pi,0.96*k_s_hat_binned(n)*[1 1],'k','linewidth',0.5)
    plot([D_k_50th_binned(ind_match_s(n),n) 3*pi/4]*180/pi,1.04*k_s_hat_binned(n)*[1 1],'k','linewidth',0.5)
    plot(3*pi/4*180/pi,k_s_hat_binned(n),'s','markersize',7,'markerfacecolor',cmap_binned(n,:),'markeredgecolor','k','linewidth',0.5)
    S_ds_fill = fill(S_ds_id_range(i)+[-1 -1 1 1]*2.5,[k_ds_hat(1,i) k_ds_hat(5,i) k_ds_hat(5,i) k_ds_hat(1,i)],cmap_binned(i,:));
    S_ds_fill.FaceAlpha = 0.5;
end
plot(D_k_50th_binned*180/pi,k_hat_binned,'k-','linewidth',3)
for i = 1:length(wave_age_centers)
    n = length(wave_age_centers)-i+1;
    plot(D_k_50th_binned(:,n)*180/pi,k_hat_binned,'Color',cmap_binned(n,:),'linewidth',2)
end
hold off
colororder(cmap_binned)
clim(clims)
box on
ylim(k_hat_lims)
xlim(halfwidth_lims)
xlabel('\theta_{halfwidth} [\circ]')
ylabel('$\mathrm{\hat{k}\equiv k u_*^2/g}$','Interpreter','LaTeX')
ax_struc(3).ax = gca;

nexttile(4)
hold on
for i = 1:length(wave_age_centers)
    n = length(wave_age_centers)-i+1;
    plot([delta_k_binned(ind_match_s(n),n) -0.5],k_s_hat_binned(n)*[1 1],'Color',[cmap_binned(n,:) alpha_vec(n)],'linewidth',3)
    plot([delta_k_binned(ind_match_s(n),n) -0.5],0.96*k_s_hat_binned(n)*[1 1],'k','linewidth',0.5)
    plot([delta_k_binned(ind_match_s(n),n) -0.5],1.04*k_s_hat_binned(n)*[1 1],'k','linewidth',0.5)
    plot(-0.5,k_s_hat_binned(n),'s','markersize',7,'markerfacecolor',cmap_binned(n,:),'markeredgecolor','k','linewidth',0.5)
end
plot([0 0],k_hat_lims,'--','Color',0.5*[1 1 1],'linewidth',2)
plot(delta_k_binned,k_hat_binned,'k-','linewidth',3)
for i = 1:length(wave_age_centers)
    n = length(wave_age_centers)-i+1;
    plot(delta_k_binned(:,n),k_hat_binned,'Color',cmap_binned(n,:),'linewidth',2)
end
hold off
colororder(cmap_binned)
box on
ylim(k_hat_lims)
xlim(Delta_lims)
xlabel('\Delta [rad]')
ylabel('$\mathrm{\hat{k}\equiv k u_*^2/g}$','Interpreter','LaTeX')
ax_struc(4).ax = gca;
clim(clims)
colormap(cmap_binned)

tile_cleaner(ax_struc,tlayout)
tlayout.TileSpacing = 'tight';

cbar = colorbar(ax_struc(3).ax);

cbar.Layout.Tile = 'east';
cbar.Layout.TileSpan = [2 2];
cbar.Ticks = clims(1):d_wave_age:clims(end);
set(get(cbar,'Title'),'String','c_E/u_*')

boxdims = [0.42 0.593 0.185 0.10];
str = 'upper limit of equilibrium range';
a = annotation('textbox',boxdims,'String',str);
a.BackgroundColor = 'w';
a.FontSize = fsize*1;
a.HorizontalAlignment = 'center';
a.VerticalAlignment = 'middle';

boxdims = [0.42 0.215 0.185 0.10];
str = {'upper limit of','saturation range'};
a = annotation('textbox',boxdims,'String',str);
a.BackgroundColor = 'w';
a.FontSize = fsize*1;
a.HorizontalAlignment = 'center';
a.VerticalAlignment = 'middle';

for n = 1:4
    if mod(n,2)
        ax_struc(n).ax.XTick = dir_ticks;
    else
        ax_struc(n).ax.XTick = Delta_ticks;
    end
    ax_struc(n).ax.YScale = 'log';
end

for n = 1:length(labels)
    nexttile(n)
    text(text_x,text_y,labels{n},'FontSize',fsize,'HorizontalAlignment','center','Units','normalized')
end

figure(fignum+1);clf
tlayout = tiledlayout(2,2);
ax_struc = struct();

nexttile(1)
hold on
plot([-180 180],[1 1],'k:','linewidth',2)
plot(D_f_50th_binned*180/pi,k_hat_disp_binned./k_n_hat_binned,'k-','linewidth',3)
for i = 1:length(wave_age_centers)
    n = length(wave_age_centers)-i+1;
    plot(D_f_50th_binned(:,n)*180/pi,k_hat_binned./k_n_hat_binned(n),'Color',cmap_binned(n,:),'linewidth',2)
    S_ds_fill = fill(S_ds_id_range(n)+[-1 -1 1 1]*2.5,[k_ds_hat(1,n) k_ds_hat(5,n) k_ds_hat(5,n) k_ds_hat(1,n)]/k_n_hat_binned(n),cmap_binned(n,:));
    S_ds_fill.FaceAlpha = 0.5;
    sg_norm_fill =  fill([100 100 125 125],[0.95 1.05 1.05 0.95]*k_sg_norm_n_s_binned(1,n),cmap_binned(n,:));
    sg_norm_fill.FaceAlpha = 0.5;
    sg_norm_fill.LineStyle = '-';
end
text(152,4.5e1,{'short gravity'},'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',fsize,'Color',0.35*[1 1 1])
hold off
colororder(cmap_binned)
box on
ylim([1e-1 1e2])
xlim(halfwidth_lims)
xlabel('\theta_{halfwidth} [\circ]')
ylabel('$\mathrm{\hat{k}_{disp}/\hat{k}_n}$','Interpreter','LaTeX')
ax_struc(1).ax = gca;

nexttile(2)
hold on
plot([-10 10],[1 1],'k:','linewidth',2)
plot([0 0],[1e-1 1e3],'--','Color',0.5*[1 1 1],'linewidth',2)
plot(delta_f_binned,k_hat_disp_binned./k_n_hat_binned,'k-','linewidth',3)
for i = 1:length(wave_age_centers)
    n = length(wave_age_centers)-i+1;
    plot(delta_f_binned(:,n),k_hat_binned./k_n_hat_binned(n),'Color',cmap_binned(n,:),'linewidth',2)
end
hold off
box on
ylim([1e-1 1e2])
xlim(Delta_lims)
xlabel('\Delta')
ylabel('$\mathrm{\hat{k}_{disp}/\hat{k}_n}$','Interpreter','LaTeX')
ax_struc(2).ax = gca;
clim(clims)
colormap(cmap_binned)

nexttile(3)
hold on
plot([-180 180],[1 1],'k:','linewidth',2)
plot(D_k_50th_binned*180/pi,k_hat_binned./k_s_hat_binned,'k-','linewidth',3)
for i = 1:length(wave_age_centers*180/pi)
    n = length(wave_age_centers)-i+1;
    plot(D_k_50th_binned(:,n)*180/pi,k_hat_binned./k_s_hat_binned(n),'Color',cmap_binned(n,:),'linewidth',2)
    S_ds_fill = fill(S_ds_id_range(n)+[-1 -1 1 1]*2.5,[k_ds_hat(1,n) k_ds_hat(5,n) k_ds_hat(5,n) k_ds_hat(1,n)]/k_s_hat_binned(n),cmap_binned(n,:));
    S_ds_fill.FaceAlpha = 0.5;
    sg_norm_fill =  fill([100 100 125 125],[0.95 1.05 1.05 0.95]*k_sg_norm_n_s_binned(2,n),cmap_binned(n,:));
    sg_norm_fill.FaceAlpha = 0.5;
    sg_norm_fill.LineStyle = '-';
end
text(152,1e1,{'short gravity'},'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',fsize,'Color',0.35*[1 1 1])
hold off
colororder(cmap_binned)
box on
ylim(k_hat_norm_lims)
xlim(halfwidth_lims)
xlabel('\theta_{halfwidth} [\circ]')
ylabel('$\mathrm{\hat{k}/\hat{k}_s}$','Interpreter','LaTeX')
clim(clims)
colormap(cmap_binned)
ax_struc(3).ax = gca;

nexttile(4)
hold on
plot([-10 10],[1 1],'k:','linewidth',2)
plot([0 0],[1e-2 1e3],'--','Color',0.5*[1 1 1],'linewidth',2)
plot(delta_k_binned,k_hat_binned./k_s_hat_binned,'k-','linewidth',3)
for i = 1:length(wave_age_centers)
    n = length(wave_age_centers)-i+1;
    plot(delta_k_binned(:,n),k_hat_binned./k_s_hat_binned(n),'Color',cmap_binned(n,:),'linewidth',2)
    gc_norm_fill =  fill([0.05 0.05 0.4 0.4],[0.95 1.05 1.05 0.95]*k_gc_norm_n_s_binned(2,n),cmap_binned(n,:));
    gc_norm_fill.FaceAlpha = 0.5;
    gc_norm_fill.LineStyle = '-';
end
text(-0.4,3.5e1,{'gravity-capillary'},'HorizontalAlignment','center','FontSize',fsize,'Color',0.35*[1 1 1])
hold off
colororder(cmap_binned)
box on
ylim(k_hat_norm_lims)
xlim(Delta_lims)
xlabel('\Delta')
ylabel('$\mathrm{\hat{k}/\hat{k}_s}$','Interpreter','LaTeX')
ax_struc(4).ax = gca;

cbar = colorbar(ax_struc(3).ax);

cbar.Layout.Tile = 'east';
cbar.Layout.TileSpan = [2 2];
cbar.Ticks = clims(1):d_wave_age:clims(end);
set(get(cbar,'Title'),'String','c_E/u_*')

tile_cleaner(ax_struc,tlayout)
tlayout.TileSpacing = 'tight';

boxdims = [0.43 0.68 0.16 0.06];
str = 'upper limit of equilibrium range';
a = annotation('textbox',boxdims,'String',str);
a.BackgroundColor = 'w';
a.FontSize = fsize*1;
a.HorizontalAlignment = 'center';
a.VerticalAlignment = 'middle';

boxdims = [0.44 0.28 0.15 0.06];
str = 'upper limit of saturation range';
a = annotation('textbox',boxdims,'String',str);
a.BackgroundColor = 'w';
a.FontSize = fsize*1;
a.HorizontalAlignment = 'center';
a.VerticalAlignment = 'middle';

for n = 1:4
    if mod(n,2)
        ax_struc(n).ax.XTick = dir_ticks;
    else
        ax_struc(n).ax.XTick = Delta_ticks;
    end
    ax_struc(n).ax.YScale = 'log';
end

for n = 1:length(labels)
    nexttile(n)
    text(text_x,text_y,labels{n},'FontSize',fsize,'HorizontalAlignment','center','Units','normalized')
end