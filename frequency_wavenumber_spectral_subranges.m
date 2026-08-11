%
function frequency_wavenumber_spectral_subranges(fignum,fsize)

cmap = flipud(spectral(7));
violet = cmap(1,:);
teal = cmap(2,:);
crimson = cmap(7,:);

supporting_nc_name = 'data/ASIT2019_supporting_environmental_observations.nc';

% Wind forcing per the source selected in the aa_step01 preamble
wind = wind_forcing(supporting_nc_name);
EC_U10_m_s = wind.U10;
EC_ustar_m_s = wind.ustar;

% Oncoming-wind gate
% wind_dir_deg_coming_from = ncread(supporting_nc_name,'COARE_Wdir');
% mask = oncoming_wind_mask(wind_dir_deg_coming_from);
%
% EC_ustar_m_s = mask.*EC_ustar_m_s;
% EC_U10_m_s = mask.*EC_U10_m_s;

U_sfc_mag_m_s = ncread(supporting_nc_name,'U_sfc_mag_m_s');

load('data/frequency_spect_range_limits.mat')
load('data/wavenumber_spect_range_limits.mat')

load('data/ASIT2019_combined_frequency_slope_spectra.mat')

water_depth_m = 15;

f_p = sum(f_Hz_combined.*F_f_block.^4,1,'omitnan')./sum(F_f_block.^4,1,'omitnan');
f_p = f_p(:);

f_eq_start = f_eq_start(:);
f_eq_end = f_eq_end(:);
f_sat_end = f_sat_end(:);

% Equilibrium-range limits are stored dispersion-derived (no current). k_sat_end
% is stored as the direct k-space extraction, so its dispersion counterpart is
% computed here -- panel (c) contrasts the two.
k_eq_start_disp = k_eq_start(:);
k_eq_end_disp = k_eq_end(:);
k_sat_end_disp = dispersion_wavenumber(f_sat_end,water_depth_m);

% Wave age phase speed: linear dispersion of f_p (no current)
[c_p,~] = lindisp_with_current(2*pi*f_p,water_depth_m,0);
c_p = c_p(:);

g = 9.81;

omega_eq_start_norm = 2*pi*f_eq_start(:).*EC_U10_m_s/g;
omega_eq_end_norm = 2*pi*f_eq_end(:).*EC_U10_m_s/g;
omega_sat_end_norm = 2*pi*f_sat_end(:).*EC_U10_m_s/g;
omega_sat_end_disp_norm = sqrt(k_sat_end(:)*g).*EC_U10_m_s/g;

k_eq_start_disp_norm = k_eq_start_disp.*EC_ustar_m_s.^2/g;
k_eq_end_disp_norm = k_eq_end_disp.*EC_ustar_m_s.^2/g;
k_sat_end_disp_norm = k_sat_end_disp.*EC_ustar_m_s.^2/g;
k_sat_end_norm = k_sat_end(:).*EC_ustar_m_s.^2/g;

text_x = 0.94;
text_y = 0.95;
labels = {'(a)','(b)','(c)','(d)'};

scat_size = 64;
lw_thick = 3;
lw_thin = 2;

s = load('data/global_figure_settings.mat');
wave_age_edges = wave_age_bin_edges();

wave_age = c_p./EC_ustar_m_s;

wave_age_lims = [wave_age_edges(1:end-1)' wave_age_edges(2:end)'];
wave_age_mean = NaN*wave_age_lims(:,1);
f_block = NaN*ones(size(wave_age_lims,1),3);
omega_block = f_block;
k_block = wave_age_mean;
k_block_disp = NaN*ones(size(wave_age_lims,1),3);
mask = 0*EC_U10_m_s+1;
mask(isnan(f_eq_start)|isnan(f_eq_end) | isnan(f_sat_end)) = NaN;
for n = 1:size(wave_age_lims,1)
    inds_within = wave_age >= wave_age_lims(n,1) & wave_age < wave_age_lims(n,2);
    wave_age_mean(n) = median(mask(inds_within).*wave_age(inds_within),'all','omitnan');
    f_block(n,:) = [median(f_eq_start(inds_within),'all','omitnan') median(f_eq_end(inds_within),'all','omitnan') median(f_sat_end(inds_within),'all','omitnan')];
    omega_block(n,:) = [median(omega_eq_start_norm(inds_within),'all','omitnan') median(omega_eq_end_norm(inds_within),'all','omitnan') median(omega_sat_end_norm(inds_within),'all','omitnan')];
    k_block(n) = median(k_sat_end_norm(inds_within),'all','omitnan');
    k_block_disp(n,:) = [median(k_eq_start_disp_norm(inds_within),'all','omitnan') median(k_eq_end_disp_norm(inds_within),'all','omitnan') median(k_sat_end_disp_norm(inds_within),'all','omitnan')];
end

fA = 0.2;

wave_age_lims = s.wave_age_lims;
wc = wave_age_color();

figure(fignum);clf
tlayout = tiledlayout(3,1);
ax_struc = struct();

nexttile(1)
hold on
h_eq_start = scatter(wave_age,f_eq_start,scat_size,violet,'filled');
h_eq_end = scatter(wave_age,f_eq_end,scat_size,teal,'filled');
h_sat_end = scatter(wave_age,f_sat_end,scat_size,crimson,'filled');
h_eq_start_binned = scatter(wave_age_mean,f_block(:,1),2*scat_size,violet);
h_eq_end_binned = scatter(wave_age_mean,f_block(:,2),2*scat_size,teal);
h_sat_end_binned = scatter(wave_age_mean,f_block(:,3),2*scat_size,crimson);
hold off
box on
xlim(wave_age_lims)
set(gca,'XScale',wc.axis_scale,'XTick',wc.axis_ticks,'XMinorTick','off')
ylim([0 3])
xlabel('$\mathrm{c_p/u_*}$','Interpreter','LaTeX')
ylabel('$\mathrm{f\ [Hz]}$','Interpreter','LaTeX')
H = [h_eq_start_binned h_eq_end_binned h_sat_end_binned];
L = {'start of eq. range','eq. \rightarrow sat.','end of sat. range'};
legend(H,L,'location','northwest')

h_eq_start_binned.Marker = 's';
h_eq_end_binned.Marker = 's';
h_sat_end_binned.Marker = 's';

h_eq_start_binned.LineWidth = lw_thin;
h_eq_end_binned.LineWidth = lw_thin;
h_sat_end_binned.LineWidth = lw_thin;

h_eq_start.MarkerFaceAlpha = fA;
h_eq_end.MarkerFaceAlpha = fA;
h_sat_end.MarkerFaceAlpha = fA/2;

ax_struc(1).ax = gca;

nexttile(2)
hold on
h_eq_start = scatter(wave_age,omega_eq_start_norm,scat_size,violet,'filled');
h_eq_end = scatter(wave_age,omega_eq_end_norm,scat_size,teal,'filled');
h_sat_end = scatter(wave_age,omega_sat_end_norm,scat_size,crimson,'filled');
h_eq_start_binned = scatter(wave_age_mean,omega_block(:,1),2*scat_size,violet);
h_eq_end_binned = scatter(wave_age_mean,omega_block(:,2),2*scat_size,teal);
h_sat_end_binned = scatter(wave_age_mean,omega_block(:,3),2*scat_size,crimson);
hold off
box on
xlim(wave_age_lims)
set(gca,'XScale',wc.axis_scale,'XTick',wc.axis_ticks,'XMinorTick','off')
ylim([0 10])
xlabel('$\mathrm{c_p/u_*}$','Interpreter','LaTeX')
ylabel('$\mathrm{\hat{\omega}\equiv\omega U_{10}/g}$','Interpreter','LaTeX')

h_eq_start_binned.Marker = 's';
h_eq_end_binned.Marker = 's';
h_sat_end_binned.Marker = 's';

h_eq_start_binned.LineWidth = lw_thin;
h_eq_end_binned.LineWidth = lw_thin;
h_sat_end_binned.LineWidth = lw_thin;

h_eq_start.MarkerFaceAlpha = fA;
h_eq_end.MarkerFaceAlpha = fA;
h_sat_end.MarkerFaceAlpha = fA/2;

ax_struc(2).ax = gca;

nexttile(3)
hold on
h_eq_start = scatter(wave_age,k_eq_start_disp_norm,scat_size,violet,'filled');
h_eq_end = scatter(wave_age,k_eq_end_disp_norm,scat_size,teal,'filled');
h_sat_end = scatter(wave_age,k_sat_end_disp_norm,scat_size,crimson,'filled');
h_eq_start_binned = scatter(wave_age_mean,k_block_disp(:,1),2*scat_size,violet);
h_eq_end_binned = scatter(wave_age_mean,k_block_disp(:,2),2*scat_size,teal);
h_sat_end_binned = scatter(wave_age_mean,k_block_disp(:,3),2*scat_size,crimson);
h_sat_end_true_binned = scatter(wave_age_mean,k_block,2*scat_size,crimson,'filled');
plot([15 70],5e-2*[15 70].^-1,'k--','linewidth',2)
hold off
box on
xlim(wave_age_lims)
set(gca,'XScale',wc.axis_scale,'XTick',wc.axis_ticks,'XMinorTick','off')
ylim([1e-4 1e-1]*3)
xlabel('$\mathrm{c_p/u_*}$','Interpreter','LaTeX')
ylabel('$\mathrm{\hat{k}\equiv k u_*^2/g}$','Interpreter','LaTeX')
text(50,0.15,'direct','Color',crimson,'FontWeight','bold','FontSize',fsize,'HorizontalAlignment','center')
text(16,1.5e-3,'$\propto u_*^{-1}$','Interpreter','LaTeX','FontSize',fsize,'HorizontalAlignment','Center')

h_eq_start_binned.Marker = 's';
h_eq_end_binned.Marker = 's';
h_sat_end_binned.Marker = 's';
h_sat_end_true_binned.Marker = 's';

h_eq_start_binned.LineWidth = lw_thin;
h_eq_end_binned.LineWidth = lw_thin;
h_sat_end_binned.LineWidth = lw_thin;

h_eq_start.MarkerFaceAlpha = fA;
h_eq_end.MarkerFaceAlpha = fA;
h_sat_end.MarkerFaceAlpha = fA/2;

ax_struc(3).ax = gca;
ax_struc(3).ax.YScale = 'log';
ax_struc(3).ax.YTick = 10.^(-4:1:1);

for n = 1:3
    nexttile(n)
    text(text_x,text_y,labels{n},'Units','normalized','FontSize',fsize,'HorizontalAlignment','center')
end

tile_cleaner(ax_struc,tlayout)

tlayout.TileSpacing = 'tight';
