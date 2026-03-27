% 
function scaled_spectra_and_transition_wavenumbers(fignum,fsize)


in_nc_name = 'data/ASIT2019_supporting_environmental_observations.nc';

g = 9.81;

water_depth_m = 15;

EC_ustar_m_s = ncread(in_nc_name,'EC_ustar_m_s');

U_sfc_mag_m_s = ncread(in_nc_name,'U_sfc_mag_m_s');

load('data/ASIT2019_combined_wavenumber_elevation_spectra.mat')

load('data/ASIT2019_combined_frequency_slope_spectra.mat')

f_p = sum(f_Hz_combined.*F_f_block.^4,1,'omitnan')./sum(F_f_block.^4,1,'omitnan');
f_p = f_p(:);

load('data/frequency_spect_range_limits.mat')
load('data/wavenumber_spect_range_limits.mat')

load('data/LenainMelville2017_kn.mat')

cmap = viridis(7);
cerulean = cmap(3,:);

f_eq_end = f_eq_end(:);

k_eq_end_disp = NaN*f_p;

for n = 1:length(f_eq_end)

    try

        [c,~] = lindisp_with_current(2*pi*f_eq_end(n),water_depth_m,U_sfc_mag_m_s(n));
        k_eq_end_disp(n) = 2*pi*f_eq_end(n)./c;

    end

end

in_nc_name = 'data/ASIT2019_supporting_environmental_observations.nc';

wse_m_Riegl = ncread(in_nc_name,'wse_m_Riegl');
Hs = 4*squeeze(std(wse_m_Riegl,[],2,'omitnan'));
Hs(Hs<0.2) = NaN;
Hs_median = median(Hs,2);
ustar_norm = EC_ustar_m_s./sqrt(g*Hs_median);

d_ustar_norm_centers = 0.03;
ustar_norm_centers = 0.03:d_ustar_norm_centers:0.15; % retain 90% of points
% inds_consider = ustar_norm >= ustar_norm_centers(1)-d_ustar_norm_centers/2 & ustar_norm < ustar_norm_centers(end)+d_ustar_norm_centers/2;
% frac_retain = sum(inds_consider)/sum(~isnan(ustar_norm));

cmap = (matter(length(ustar_norm_centers)+1));
cmap(1,:) = [];

B_k_binned = NaN*ones(length(k_rad_m_combined),length(ustar_norm_centers));
k_binned = B_k_binned;
k_n_binned = NaN*ustar_norm_centers;
k_n_std_e = k_n_binned;

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

    k_n_binned(n) = mean(k_eq_end_disp(inds_consider),'omitnan');
    k_n_std_e(n) = std(k_eq_end_disp(inds_consider),'omitnan')/sqrt(sum(inds_consider));

end

our_color = cerulean;
LM_color = [1 1 1]*0.4;
% our_color = [0.6 0 0];
% our_color = [104 71 141]/255;
% LM_color = [0 124 124]/255;
msize = 7;
scat_size = 0.8*msize^2;

fA = 0.25;

% low cutoff: k ~ 0.03 rad/m (f~0.05 Hz)
start_index = 67;
% high cutoff: k = 371 rad/m (k of gravity-capillary minimum phase speed)
end_index = 1196;

figure(fignum);clf
tlayout = tiledlayout(2,1);

nexttile(1)
hold on
loglog(k_binned(start_index:end_index,:),B_k_binned(start_index:end_index,:),'k','linewidth',4)
loglog(k_binned(start_index:end_index,:),B_k_binned(start_index:end_index,:),'linewidth',3);colororder(cmap)
hold off
box on
colormap(cmap)
cbar = colorbar;
cbar.Location = 'northoutside';
set(get(cbar,'Label'),'String','$u_*/\sqrt{gH_s}$','Interpreter','LaTeX')
cbar.Ticks = [ustar_norm_centers(1)-d_ustar_norm_centers/2 ustar_norm_centers+d_ustar_norm_centers/2];
clim([ustar_norm_centers(1) ustar_norm_centers(end)] + [-1 1]*d_ustar_norm_centers/2)
xlim([0.999e-5 1e1])
ylim([1e-6 1e-1])
xlabel('$ku_*^2/g$','Interpreter','LaTeX')
ylabel('$B(k)\ \mathrm{[rad]}$','Interpreter','LaTeX')
ax = gca;
ax.XTick = 10.^(-5:1:2);
ax.XScale = 'log';
ax.YScale = 'log';

nexttile(2)
hold on
% h_LM2017 = scatter(LM2017_ustar_norm,LM2017_kn,scat_size,LM_color,'filled');
h_LM2017 = plot(LM2017_ustar_norm,LM2017_kn,'o','markersize',msize,'markerfacecolor','w','markeredgecolor',LM_color,'linewidth',2);
h_ours_all = scatter(ustar_norm,k_eq_end_disp,scat_size,our_color,'filled');
for i = 1:length(k_n_binned)
    plot(ustar_norm_centers(i)*[1 1],k_n_binned(i)+[-1 1]*k_n_std_e(i)*1.96,'-','Color',our_color,'linewidth',2)
end
h_ours_binned = scatter(ustar_norm_centers,k_n_binned,scat_size*3,our_color,'filled');
hold off
box on
xlim([0 0.2])
ylim([0 10])
xlabel('$u_*/\sqrt{gH_s}$','Interpreter','LaTeX')
ylabel('$k_n\ \mathrm{[rad}\ \mathrm{m}^{-1}]$','Interpreter','LaTeX')

h_ours_all.MarkerFaceAlpha = fA;

h_ours_binned.Marker = 's';

H = [h_ours_binned h_LM2017];
L = {'present study (binned)','Lenain & Melville [2017]'};
legend(H,L)

labels = {'(a)','(b)'};
for n = 1:2
    nexttile(n)
    text(0.05,0.95,labels{n},'Units','Normalized','HorizontalAlignment','Center','FontSize',fsize)
end

tlayout.TileSpacing = 'compact';