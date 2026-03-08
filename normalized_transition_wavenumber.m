%
function normalized_transition_wavenumber(fignum,fsize)

in_nc_name = 'data/ASIT2019_supporting_environmental_observations.nc';

g = 9.81;

water_depth_m = 15;

EC_ustar_m_s = ncread(in_nc_name,'EC_ustar_m_s');

U_sfc_mag_m_s = ncread(in_nc_name,'U_sfc_mag_m_s');

load('data/ASIT2019_combined_wavenumber_elevation_spectra.mat')

load('data/frequency_spect_range_limits.mat')
load('data/wavenumber_spect_range_limits.mat')

f_p = f_E;

f_p = f_p(:);
f_eq_start = f_eq_start(:);
f_eq_end = f_eq_end(:);

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

B_k_block = k_rad_m_combined.^3.*F_k_block;
E_k_block = k_rad_m_combined.^2.5.*F_k_block;

B = NaN*EC_ustar_m_s;
E = B;

for n = 1:length(B)

    inds_consider = k_rad_m_combined >= k_eq_end_disp(n) & k_rad_m_combined <= k_sat_end(n);

    B(n) = mean(B_k_block(inds_consider,n),'all','omitnan');

    inds_consider = k_rad_m_combined >= k_eq_start_disp(n) & k_rad_m_combined <= k_eq_end_disp(n);

    E(n) = mean(E_k_block(inds_consider,n),'all','omitnan');

end

beta = 2*E*g^0.5./EC_ustar_m_s;

%

cmap = viridis(7);
cerulean = cmap(3,:);

our_color = cerulean;
RM_color = [1 1 1]*0.6;
% our_color = [104 71 141]/255;
% RM_color = [0 124 124]/255;

msize = 7;
lw_thin = 0.5;
lw_thick = 1.5;

fA = 0.25;

load('data/RM1010_ko.mat')

load('data/RM2010_ko_kp_wave_age.mat')

wave_age_full = 2*pi*f_E(:)./k_p_disp(:)./EC_ustar_m_s(:);

ko_HW_full = (2*B./beta).^2*g./EC_ustar_m_s(:).^2;

%%

wave_age_bins = 10:10:100;
cmap = flipud(magma(length(wave_age_bins)-1));

figure(fignum);clf
tlayout = tiledlayout(1,2);
ax_struc = struct();

nexttile(1)
hold on
h_RF05 = plot(RF05(:,1),RF05(:,2),'o','markerfacecolor',RM_color,'markeredgecolor',RM_color,'markersize',msize,'linewidth',lw_thin);
h_RF07 = plot(RF07(:,1),RF07(:,2),'d','markerfacecolor',RM_color,'markeredgecolor',RM_color,'markersize',msize,'linewidth',lw_thin);
h_RF09_k0 = plot(RF09_k0(:,1),RF09_k0(:,2),'p','markerfacecolor',RM_color,'markeredgecolor',RM_color,'markersize',msize,'linewidth',lw_thin);
h_RF09_k2 = plot(RF09_k2(:,1),RF09_k2(:,2),'p','markerfacecolor','w','markeredgecolor',RM_color,'markersize',msize,'linewidth',lw_thick);
h_RF10_k0 = plot(RF10_k0(:,1),RF10_k0(:,2),'s','markerfacecolor',RM_color,'markeredgecolor',RM_color,'markersize',msize,'linewidth',lw_thin);
h_RF10_k2 = plot(RF10_k2(:,1),RF10_k2(:,2),'s','markerfacecolor','w','markeredgecolor',RM_color,'markersize',msize,'linewidth',lw_thick);
h_RM2010_CI = fill([RM2010_wave_age; flipud(RM2010_wave_age)],[RM2010_ko_kp_CI(:,1); flipud(RM2010_ko_kp_CI(:,2))],RM_color);
h_RM2010_fit = plot(RM2010_wave_age,RM2010_ko_kp,'Color',RM_color,'linewidth',3);
h_ours = plot(wave_age_full,k_eq_end_disp./k_p_disp(:),'o','markerfacecolor',our_color,'markeredgecolor',our_color,'markersize',msize,'linewidth',lw_thin);
hold off
box on
ax_struc(1).ax = gca;
ax_struc(1).ax.XScale = 'log';
ax_struc(1).ax.YScale = 'log';
xlim([10 100])
ylim([1 100])
xlabel('c_E/u_*')
ylabel('k_n/k_E, obs.')

H = [h_ours h_RM2010_fit];
L = {'present study','Romero & Melville [2010]'};
legend(H,L,'Location','southwest')

h_RM2010_CI.FaceAlpha = fA;

nexttile(2)
hold on
plot(ko_HW_full./k_p_disp,k_eq_end_disp./k_p_disp,'o','markerfacecolor','k','markeredgecolor','k','markersize',msize)
plot([1 100],[1 100],'--','Color',0.5*[1 1 1],'linewidth',2)
scatter(ko_HW_full./k_p_disp,k_eq_end_disp./k_p_disp,0.8*msize^2,(wave_age_full),'filled')
hold off
box on
cbar = colorbar;
clim([10 100])
colormap(cmap)
set(get(cbar,'Title'),'String','c_E/u_*')
ax_struc(2).ax = gca;
ax_struc(2).ax.XScale = 'log';
ax_struc(2).ax.YScale = 'log';
xlim([1 100])
ylim([1 100])
xlabel('k_n/k_E, Hwang & Wang [2001]')
ylabel('k_n/k_E, obs.')

ax_struc(1).ax.XTick = 10.^[1 2];
ax_struc(1).ax.XTickLabel = {'10','100'};

ax_struc(1).ax.YTick = 10.^[0 1 2];
ax_struc(1).ax.YTickLabel = {'1','10','100'};

ax_struc(2).ax.XTick = 10.^[0 1 2];
ax_struc(2).ax.XTickLabel = {'1','10','100'};

tile_cleaner(ax_struc,tlayout)

labels = {'(a)','(b)'};
for n = 1:2
    nexttile(n)
    pbaspect([1 1 1])
    text(0.05,0.95,labels{n},'Units','Normalized','HorizontalAlignment','Center','FontSize',fsize)
end

tlayout.TileSpacing = 'tight';

ko_HW_full = ko_HW_full(:);
k_eq_end_disp = k_eq_end_disp(:);

inds_keep = ~isnan(ko_HW_full) & ~isnan(k_eq_end_disp);
[r2,rmse] = rsquare(ko_HW_full(inds_keep),k_eq_end_disp(inds_keep))
