%
function normalized_transition_wavenumber(fignum,fsize)

% Conventions: k_n, k_p, and the subrange limits all come from the frequency
% spectrum through linear dispersion (no current), so the ratio is internally
% consistent; wave age c_p/u* uses c_p from the Young (1995) peak frequency f_p

supporting_nc_name = 'data/ASIT2019_supporting_environmental_observations.nc';

g = 9.81;

water_depth_m = 15;

% Wind forcing per the source selected in the aa_step01 preamble
wind = wind_forcing(supporting_nc_name);
EC_ustar_m_s = wind.ustar;

U_sfc_mag_m_s = ncread(supporting_nc_name,'U_sfc_mag_m_s');

load('data/ASIT2019_combined_wavenumber_elevation_spectra.mat')

load('data/frequency_spect_range_limits.mat')
load('data/wavenumber_spect_range_limits.mat')

load('data/ASIT2019_combined_frequency_slope_spectra.mat')

f_p = sum(f_Hz_combined.*F_f_block.^4,1,'omitnan')./sum(F_f_block.^4,1,'omitnan');
f_p = f_p(:);

f_eq_start = f_eq_start(:);
f_eq_end = f_eq_end(:);

% k_n and the equilibrium-range start both come from the frequency spectrum
% through linear dispersion, stored in wavenumber_spect_range_limits.mat
k_eq_end = k_eq_end(:);
k_eq_start_disp = k_eq_start(:);

% Peak wavenumber from f_p through linear dispersion, matching k_n/k_eq/k_sat
k_p = dispersion_wavenumber(f_p,water_depth_m);

% Wave age phase speed: linear dispersion of f_p (no current)
[c_p,~] = lindisp_with_current(2*pi*f_p,water_depth_m,0);
c_p = c_p(:);

B_k_block = k_rad_m_combined.^3.*F_k_block;
E_k_block = k_rad_m_combined.^2.5.*F_k_block;

B = NaN*EC_ustar_m_s;
E = B;

for n = 1:length(B)

    inds_consider = k_rad_m_combined >= k_eq_end(n) & k_rad_m_combined <= k_sat_end(n);

    B(n) = mean(B_k_block(inds_consider,n),'all','omitnan');

    inds_consider = k_rad_m_combined >= k_eq_start_disp(n) & k_rad_m_combined <= k_eq_end(n);

    E(n) = mean(E_k_block(inds_consider,n),'all','omitnan');

end

beta = 2*E*g^0.5./EC_ustar_m_s;

%

% Our runs are colored by their wave-wind misalignment, theta_m-theta_wind; runs
% with no misalignment (no gated wind record) hold their place in the scatter as
% open gray circles rather than dropping out
mc = wave_wind_misalignment_color();
dtheta_m = mc.value;

no_dir_color = [1 1 1]*0.35;
% Green for the published comparison, so it separates from the present study
% rather than blending into the gray guide lines
RM_color = [0.13 0.55 0.24];

msize = 9;
lw_thin = 0.5;
lw_thick = 1.5;

fA = 0.25;

load('data/RM1010_ko.mat')

load('data/RM2010_ko_kp_wave_age.mat')

wave_age_full = c_p./EC_ustar_m_s(:);

ko_HW_full = (2*B./beta).^2*g./EC_ustar_m_s(:).^2;

%%

scat_size = 0.8*msize^2;

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
scatter(wave_age_full,k_eq_end./k_p(:),scat_size,'MarkerEdgeColor',no_dir_color,'LineWidth',lw_thin);
scatter(wave_age_full,k_eq_end./k_p(:),scat_size,dtheta_m(:),'filled','MarkerEdgeColor','k','LineWidth',lw_thin);
h_ours = plot(NaN,NaN,'o','markerfacecolor','none','markeredgecolor','k','markersize',msize,'linewidth',lw_thin);
hold off
box on
ax_struc(1).ax = gca;
ax_struc(1).ax.XScale = 'log';
ax_struc(1).ax.YScale = 'log';
colormap(ax_struc(1).ax,mc.cmap)
clim(mc.clims)
xlim([10 100])
ylim([1 200])
xlabel('$\mathrm{c_p/u_*}$','Interpreter','LaTeX')
ylabel('$\mathrm{k_n/k_p,\ obs.}$','Interpreter','LaTeX')

H = [h_ours h_RM2010_fit];
L = {'present study','Romero & Melville [2010]'};
legend(H,L,'Location','southwest')

h_RM2010_CI.FaceAlpha = fA;

nexttile(2)
hold on
plot([1 200],[1 200],'--','Color',0.5*[1 1 1],'linewidth',2)
scatter(ko_HW_full./k_p,k_eq_end./k_p,scat_size,'MarkerEdgeColor',no_dir_color,'LineWidth',lw_thin)
scatter(ko_HW_full./k_p,k_eq_end./k_p,scat_size,dtheta_m(:),'filled','MarkerEdgeColor','k','LineWidth',lw_thin)
hold off
box on
ax_struc(2).ax = gca;
colormap(ax_struc(2).ax,mc.cmap)
clim(mc.clims)
cbar = colorbar;
set(get(cbar,'Title'),'String',mc.label,'Interpreter','LaTeX')
cbar.Ticks = mc.ticks;
cbar.TickLabels = mc.ticklabels;
ax_struc(2).ax.XScale = 'log';
ax_struc(2).ax.YScale = 'log';
xlim([1 200])
ylim([1 200])
xlabel('$\mathrm{k_n/k_p,\ Hwang\ \&\ Wang\ [2001]}$','Interpreter','LaTeX')
ylabel('$\mathrm{k_n/k_p,\ obs.}$','Interpreter','LaTeX')

% Wave-age axis: integer, explicitly log-spaced ticks, as in fig 6
ax_struc(1).ax.XTick = nice_log_ticks(ax_struc(1).ax.XLim);
ax_struc(1).ax.XTickLabel = arrayfun(@(v) sprintf('%d',v),ax_struc(1).ax.XTick,'UniformOutput',false);
ax_struc(1).ax.XMinorTick = 'off';

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
k_eq_end = k_eq_end(:);

inds_keep = ~isnan(ko_HW_full) & ~isnan(k_eq_end);
[r2,rmse] = rsquare(ko_HW_full(inds_keep),k_eq_end(inds_keep))

% Manuscript fragment: Hwang & Wang model vs observed k_n agreement
write_tex_macros('figs/tex/values_fig_knkp.tex', ...
    {'FigKnKpRsq','FigKnKpRmse'}, ...
    {sprintf('%.2f',r2),sprintf('%.2f',rmse)})
