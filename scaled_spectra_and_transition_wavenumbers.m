% 
function scaled_spectra_and_transition_wavenumbers(fignum,fsize)


supporting_nc_name = 'data/ASIT2019_supporting_environmental_observations.nc';

g = 9.81;

% Wind forcing per the source selected in the aa_step01 preamble
wind = wind_forcing(supporting_nc_name);
EC_ustar_m_s = wind.ustar;

load('data/ASIT2019_combined_wavenumber_elevation_spectra.mat')

load('data/ASIT2019_combined_frequency_slope_spectra.mat')

f_p = sum(f_Hz_combined.*F_f_block.^4,1,'omitnan')./sum(F_f_block.^4,1,'omitnan');
f_p = f_p(:);

load('data/frequency_spect_range_limits.mat')
load('data/wavenumber_spect_range_limits.mat')

load('data/LenainMelville2017_kn.mat')

% k_n taken directly from the EWDM+Pyxis combined wavenumber elevation spectrum
k_eq_end = k_eq_end(:);

supporting_nc_name = 'data/ASIT2019_supporting_environmental_observations.nc';

wse_m_Riegl = ncread(supporting_nc_name,'wse_m_Riegl');
Hs = 4*squeeze(std(wse_m_Riegl,[],2,'omitnan'));
Hs(Hs<0.2) = NaN;
Hs_median = median(Hs,2);
ustar_norm = EC_ustar_m_s./sqrt(g*Hs_median);
% ustar_norm(ustar_norm>0.11) = NaN;

d_ustar_norm_centers = 0.02;
ustar_norm_centers = 0.05:d_ustar_norm_centers:0.13;
% inds_consider = ustar_norm >= ustar_norm_centers(1)-d_ustar_norm_centers/2 & ustar_norm < ustar_norm_centers(end)+d_ustar_norm_centers/2;
% frac_retain = sum(inds_consider)/sum(~isnan(ustar_norm));

cmap = (matter(length(ustar_norm_centers)+1));
cmap(1,:) = [];

% Individual runs are colored by their wave-wind misalignment, theta_m-theta_wind
mc = wave_wind_misalignment_color();
dtheta_m = mc.value;

B_k_binned = NaN*ones(length(k_rad_m_combined),length(ustar_norm_centers));
k_binned = B_k_binned;
k_n_binned = NaN*ustar_norm_centers;
k_n_std_e = k_n_binned;
dtheta_m_binned = k_n_binned;

for n = 1:length(ustar_norm_centers)

    ustar_norm_low = ustar_norm_centers(n) - d_ustar_norm_centers/2;
    ustar_norm_high = ustar_norm_centers(n) + d_ustar_norm_centers/2;

    inds_consider = ustar_norm >= ustar_norm_low & ustar_norm < ustar_norm_high;

    k_mean = mean(k_rad_m_combined.*EC_ustar_m_s(inds_consider)'.^2/g,2,'omitnan');
    Bk_mean = mean(k_rad_m_combined.^3.*F_k_block(:,inds_consider),2,'omitnan');

    [k_mean,order] = sort(k_mean);
    Bk_mean = Bk_mean(order);

    B_k_binned(:,n) = Bk_mean;
    k_binned(:,n) = k_mean;

    k_n_binned(n) = median(k_eq_end(inds_consider),'omitnan');
    k_n_std_e(n) = std(k_eq_end(inds_consider),'omitnan')/sqrt(sum(inds_consider));

    % Misalignment wraps, so the bin's representative direction is a circular
    % median: the arithmetic one would place a bin split across +/-180 near zero
    dtheta_m_binned(n) = circular_median(dtheta_m(inds_consider));

end

% Runs with no wave-wind misalignment (no gated wind record) hold their place in
% the scatter as open gray circles rather than dropping out
no_dir_color = [1 1 1]*0.55;
CI_color = [1 1 1]*0.3;
% Green for the published comparison, matching the treatment in figure 8
LM_color = [0.13 0.55 0.24];
msize = 9;
scat_size = 0.8*msize^2;

fA_marker = 0.85;

% low cutoff: k ~ 0.03 rad/m (f~0.05 Hz)
start_index = find(k_rad_m_combined >= 0.03,1,'first');
% high cutoff: k = 371 rad/m (k of gravity-capillary minimum phase speed)
end_index = find(k_rad_m_combined <= 371,1,'last');

figure(fignum);clf
tlayout = tiledlayout(2,1);

nexttile(1)
hold on
loglog(k_binned(start_index:end_index,:),B_k_binned(start_index:end_index,:),'k','linewidth',4)
loglog(k_binned(start_index:end_index,:),B_k_binned(start_index:end_index,:),'linewidth',3);colororder(cmap)
hold off
box on
colormap(gca,cmap)
cbar = colorbar;
cbar.Location = 'eastoutside';
set(get(cbar,'Title'),'String','$\mathrm{u_*/\sqrt{gH_s}}$','Interpreter','LaTeX')
cbar.Ticks = [ustar_norm_centers(1)-d_ustar_norm_centers/2 ustar_norm_centers+d_ustar_norm_centers/2];
clim([ustar_norm_centers(1) ustar_norm_centers(end)] + [-1 1]*d_ustar_norm_centers/2)
xlim([0.999e-5 1e1])
ylim([1e-6 1e-1])
xlabel('$\mathrm{\hat{k}\equiv ku_*^2/g}$','Interpreter','LaTeX')
ylabel('$\mathrm{B(k)\ [rad]}$','Interpreter','LaTeX')
ax = gca;
ax.XTick = 10.^(-5:1:2);
ax.YTick = 10.^(-6:1:-1);
ax.XScale = 'log';
ax.YScale = 'log';

nexttile(2)
hold on
% h_LM2017 = scatter(LM2017_ustar_norm,LM2017_kn,scat_size,LM_color,'filled');
h_LM2017 = plot(LM2017_ustar_norm,LM2017_kn,'o','markersize',msize*1.2,'markerfacecolor','w','markeredgecolor',LM_color,'linewidth',2);
h_ours_all = scatter(ustar_norm(:),k_eq_end(:),scat_size,'MarkerEdgeColor',no_dir_color,'linewidth',2);
h_ours_dir = scatter(ustar_norm(:),k_eq_end(:),scat_size,dtheta_m(:),'filled','MarkerEdgeColor','none');
% for i = 1:length(k_n_binned)
%     plot(ustar_norm_centers(i)*[1 1],k_n_binned(i)+[-1 1]*k_n_std_e(i)*1.96,'-','Color',CI_color,'linewidth',2)
% end
% h_ours_binned_outline = scatter(ustar_norm_centers,k_n_binned,scat_size*3,dtheta_m_binned,'filled','MarkerEdgeColor','k','linewidth',4);
% h_ours_binned = scatter(ustar_norm_centers,k_n_binned,scat_size*3,dtheta_m_binned,'filled','MarkerEdgeColor','w','linewidth',1.5);
% Legend proxy: the binned markers carry data in their fill, so the legend shows
% the outline alone rather than whichever color the first bin happens to take
% h_binned_proxy = plot(NaN,NaN,'s','markersize',msize*1.6,'markerfacecolor','w','markeredgecolor','k','linewidth',1.5);
hold off
box on
xlim([0 0.16])
ylim([2e-1 20])
ax = gca;
ax.YScale = 'log';
ax.YTick = [0.2 2 20];
ax.YTickLabel = {'0.2','2','20'};
colormap(ax,mc.cmap)
clim(mc.clims)
cbar_dir = colorbar;
cbar_dir.Location = 'eastoutside';
cbar_dir.Ticks = mc.ticks;
cbar_dir.TickLabels = mc.ticklabels;
set(get(cbar_dir,'Title'),'String',mc.label,'Interpreter','LaTeX')
xlabel('$\mathrm{u_*/\sqrt{gH_s}}$','Interpreter','LaTeX')
ylabel('$\mathrm{k_n\ [rad\ m^{-1}]}$','Interpreter','LaTeX')

h_ours_dir.MarkerFaceAlpha = fA_marker;

% h_ours_binned_outline.Marker = 's';
% h_ours_binned.Marker = 's';

H = [h_LM2017];
L = {'Lenain & Melville [2017]'};
legend(H,L,'Location','southwest')

labels = {'(a)','(b)'};
for n = 1:2
    nexttile(n)
    text(0.05,0.95,labels{n},'Units','Normalized','HorizontalAlignment','Center','FontSize',fsize)
end

tlayout.TileSpacing = 'compact';


% Circular median of a set of angles [deg]: the sample angle minimizing the
% summed angular distance to the rest, so a group straddling +/-180 is not
% dragged to zero the way an arithmetic median would be
function mu = circular_median(theta_deg)

theta_deg = theta_deg(isfinite(theta_deg));
theta_deg = theta_deg(:);

if isempty(theta_deg)
    mu = NaN;
    return
end

d = abs(mod(theta_deg - theta_deg' + 180,360) - 180);
[~,i] = min(sum(d,1));
mu = theta_deg(i);
