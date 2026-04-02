%
function frequency_wavenumber_spectral_subranges(fignum,fsize)

cmap = flipud(spectral(7));
violet = cmap(1,:);
teal = cmap(2,:);
crimson = cmap(7,:);

in_nc_name = 'data/ASIT2019_supporting_environmental_observations.nc';

kappa = 0.4;

EC_U_m_s = ncread(in_nc_name,'EC_U_m_s');
EC_ustar_m_s = ncread(in_nc_name,'EC_ustar_m_s');
EC_z_m_above_water = ncread(in_nc_name,'EC_z_m_above_water');
EC_z0_m = EC_z_m_above_water.*exp(-kappa*EC_U_m_s./EC_ustar_m_s);
EC_U10_m_s = EC_ustar_m_s/kappa.*log(10./EC_z0_m);

U_sfc_mag_m_s = ncread(in_nc_name,'U_sfc_mag_m_s');

load('data/frequency_spect_range_limits.mat')
load('data/wavenumber_spect_range_limits.mat')

load('data/ASIT2019_combined_frequency_slope_spectra.mat')

water_depth_m = 15;

f_p = sum(f_Hz_combined.*F_f_block.^4,1,'omitnan')./sum(F_f_block.^4,1,'omitnan');
f_p = f_p(:);

f_eq_start = f_eq_start(:);
f_eq_end = f_eq_end(:);
f_sat_end = f_sat_end(:);

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

c_p = 2*pi*f_p./k_p_disp;

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

d_wave_age = 10;

wave_age = c_p./EC_ustar_m_s;

wave_age_left = 10:10:70;
wave_age_right = wave_age_left + 10;
wave_age_lims = [wave_age_left(:) wave_age_right(:)];
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

wave_age_lims = [0 80];

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
hold off
box on
xlim(wave_age_lims)
ylim([1e-4 1e-1]*3)
xlabel('$\mathrm{c_p/u_*}$','Interpreter','LaTeX')
ylabel('$\mathrm{\hat{k}\equiv k u_*^2/g}$','Interpreter','LaTeX')
text(50,0.15,'direct','Color',crimson,'FontWeight','bold','FontSize',fsize,'HorizontalAlignment','center')

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
