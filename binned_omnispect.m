function binned_omnispect(fignum,fsize)

supporting_nc_name = 'data/ASIT2019_supporting_environmental_observations.nc';
short_wave_spectra_nc_name = 'data/ASIT2019_wave_spectra_stats_timeseries_empirical_gain.nc';
long_wave_spectra_nc_name = 'data/ASIT2019_EPSS_directional_spectra.nc';

g = 9.81;

f_Hz_EPSS = ncread(long_wave_spectra_nc_name,'frequency');

% Auto-detect direction units (legacy: deg/per-deg, current: rad/per-rad)
dir_EPSS = double(ncread(long_wave_spectra_nc_name,'direction'));
if max(abs(dir_EPSS)) > 2*pi
    theta_rad_EPSS = dir_EPSS*pi/180;
    epss_per_rad_scale = 180/pi;
else
    theta_rad_EPSS = dir_EPSS;
    epss_per_rad_scale = 1;
end

f_Hz_Pyxis = ncread(short_wave_spectra_nc_name,'f_Hz');
k_rad_m_Pyxis = ncread(short_wave_spectra_nc_name,'k_rad_m');
theta_rad_Pyxis = ncread(short_wave_spectra_nc_name,'theta_rad');

kappa = 0.4;

% Wind forcing per the source selected in the aa_step01 preamble
wind = wind_forcing(supporting_nc_name);
EC_ustar_m_s = wind.ustar;
EC_U10_m_s = wind.U10;

SKTHETA = double(ncread(short_wave_spectra_nc_name,'S_k_theta'));
SFTHETA = double(ncread(short_wave_spectra_nc_name,'S_f_theta'));
FFTHETA = double(ncread(long_wave_spectra_nc_name,'F_f_d'))*epss_per_rad_scale;

SKTHETA(SKTHETA==0) = NaN;
SFTHETA(SFTHETA==0) = NaN;
FFTHETA(FFTHETA==0) = NaN;

U_E_m_s = double(ncread(supporting_nc_name,'U_E_m_s'));
U_N_m_s = double(ncread(supporting_nc_name,'U_N_m_s'));
z_m = double(ncread(supporting_nc_name,'z_m'));
tail_flag = false;
k_max = 10;

U10_high_cutoff = 15;

f_cut = 0.35;

ind_trim = find(f_Hz_Pyxis>f_cut,1,'first');

FFTHETA = FFTHETA(:,:,f_Hz_EPSS<f_cut);
f_Hz_EPSS = f_Hz_EPSS(f_Hz_EPSS<f_cut);

ind_cut = length(f_Hz_EPSS);

water_depth_m = 15;

N = length(EC_U10_m_s);

% EWDM long-wave wavenumber spectra (handoff to Pyxis direct at k = 2 rad/m)
k_EWDM = ncread(long_wave_spectra_nc_name,'wavenumber');
F_k_omni_all = ncread(long_wave_spectra_nc_name,'F_k');   % MATLAB (run, wavenumber)
keepE = k_EWDM < 2;
keepP = k_rad_m_Pyxis >= 2;

F_f_block = NaN*ones(length(f_Hz_EPSS)+length(f_Hz_Pyxis(ind_trim:end-1)),N);
S_f_block = F_f_block;
F_k_block = NaN*ones(nnz(keepE)+nnz(keepP),N);

for example_run_ind = 1:N

    try

        wdir_deg = double(ncread(supporting_nc_name,'COARE_Wdir'));

        wdir_deg = mod(wdir_deg(example_run_ind)+180,360);
        wdir_rad = pi/180*wdir_deg;

        S_k_theta = squeeze(SKTHETA(:,:,example_run_ind));

        S_f_theta = squeeze(SFTHETA(:,:,example_run_ind));
        F_f_theta_m2_Hz_rad = squeeze(FFTHETA(example_run_ind,:,:))';

        [S_f_theta,theta_rad] = convert_dirspect_to_downwind(S_f_theta,theta_rad_Pyxis,wdir_rad);
        [F_f_theta_m2_Hz_rad,theta_rad_EPSS_dw] = convert_dirspect_to_downwind(F_f_theta_m2_Hz_rad,theta_rad_EPSS,wdir_rad);
        % Resample E-PSS onto the Pyxis downwind direction grid used downstream
        F_f_theta_m2_Hz_rad = regrid_directional_spectrum(F_f_theta_m2_Hz_rad,theta_rad_EPSS_dw,theta_rad);

        [c_disp_EPSS,~] = lindisp_with_current(2*pi*f_Hz_EPSS,water_depth_m,0);
        [c_disp_Pyxis,~] = lindisp_with_current(2*pi*f_Hz_Pyxis,water_depth_m,0);

        k_disp_EPSS = 2*pi*f_Hz_EPSS./c_disp_EPSS;
        k_disp_Pyxis = 2*pi*f_Hz_Pyxis./c_disp_Pyxis;

        dtheta = median(diff(theta_rad_Pyxis));

        S_f_Pyxis = sum(S_f_theta,2,'omitnan')*dtheta;
        F_f_Pyxis = k_disp_Pyxis.^-2.*S_f_Pyxis;
        F_f_EPSS = sum(F_f_theta_m2_Hz_rad,2,'omitnan')*dtheta;

        f_Hz_combined = [f_Hz_EPSS; f_Hz_Pyxis(ind_trim:end-1)];
        F_f_combined = [F_f_EPSS; F_f_Pyxis(ind_trim:end-1)];
        S_f_combined = [k_disp_EPSS.^2.*F_f_EPSS; S_f_Pyxis(ind_trim:end-1)];

        [S_k_theta,~] = convert_dirspect_to_downwind(S_k_theta,theta_rad_Pyxis,wdir_rad);

        % Long waves: EWDM omni wavenumber elevation spectrum (k<2); short
        % waves: Pyxis direct (k>=2), F_k = sum(k.*(k^-2 S_k),theta)dtheta
        F_k_EWDM_omni = F_k_omni_all(example_run_ind,:).';
        F_k_Pyxis_omni = sum(k_rad_m_Pyxis.^-1.*S_k_theta,2,'omitnan')*dtheta;

        k_rad_m_combined = [k_EWDM(keepE); k_rad_m_Pyxis(keepP)];
        F_k_combined = [F_k_EWDM_omni(keepE); F_k_Pyxis_omni(keepP)];

        F_f_block(:,example_run_ind) = F_f_combined(:);
        S_f_block(:,example_run_ind) = S_f_combined(:);
        F_k_block(:,example_run_ind) = F_k_combined(:);

    end

end

F_f_block(:,EC_U10_m_s>U10_high_cutoff) = NaN;
F_k_block(:,EC_U10_m_s>U10_high_cutoff) = NaN;

F_f_block(F_f_block==0) = NaN;
S_f_block(S_f_block==0) = NaN;
F_k_block(F_k_block==0) = NaN;

[c_m_s_disp,~] = lindisp_with_current(2*pi*f_Hz_combined,water_depth_m,0);

k_rad_m_disp = 2*pi*f_Hz_combined./c_m_s_disp;

save('data/ASIT2019_combined_frequency_slope_spectra.mat','F_f_block','S_f_block','f_Hz_combined','-v7.3')
save('data/ASIT2019_combined_wavenumber_elevation_spectra.mat','F_k_block','k_rad_m_combined','-v7.3')

%% FIGURES

text_x = 0.05;
text_y = 0.95;
labels = {'(a)','(b)','(c)','(d)','(e)','(f)'};

f_p = sum(f_Hz_combined.*F_f_block.^4,1,'omitnan')./sum(F_f_block.^4,1,'omitnan');
[c_p,~] = lindisp_with_current(2*pi*f_p,water_depth_m,0);

waveage = c_p./EC_ustar_m_s(:);

s = load('data/global_figure_settings.mat');
[waveage_edges,waveage_centers] = wave_age_bin_edges();
nU = length(waveage_centers);

wc = wave_age_color();
Ulims = wc.clims;

cmap_binned = (magma(nU));

F_f_binned = NaN*ones(size(F_f_combined,1),nU);
S_f_binned = F_f_binned;
F_k_binned = NaN*ones(size(F_k_combined,1),nU);

% Runs actually contributing to each bin mean, i.e. in the bin AND carrying a
% spectrum, which is the number worth reporting rather than the bin occupancy
has_spectrum = any(isfinite(F_f_block),1)';

n_per_bin = NaN*ones(1,nU);

for n = 1:nU

    inds_consider = waveage >= waveage_edges(n) & waveage < waveage_edges(n+1);
    F_f_binned(:,n) = mean(F_f_block(:,inds_consider),2,'omitnan');
    S_f_binned(:,n) = mean(S_f_block(:,inds_consider),2,'omitnan');
    F_k_binned(:,n) = mean(F_k_block(:,inds_consider),2,'omitnan');

    n_per_bin(n) = sum(inds_consider & has_spectrum);

end

F_f_binned = fliplr(F_f_binned);
F_k_binned = fliplr(F_k_binned);

f_ind_cut = ind_cut;
k_ind_cut = nnz(keepE);

F_f_binned(f_ind_cut,:) = NaN;
F_k_binned(k_ind_cut,:) = NaN;

f_eq = [2e-1 5e-1];
f_sat = [5e-1 1.5e0];

k_eq = [1e-1 1e0];
k_sat = [1e0 2e1];

figure(fignum);clf
tlayout = tiledlayout(2,2);
ax_struc = struct();

nexttile(1)
hold on
loglog(f_Hz_combined,F_f_binned,'k-','linewidth',3)
loglog(f_Hz_combined,F_f_binned,'linewidth',2)
plot(f_eq,1e-2*f_eq.^-4,'k--','linewidth',2.5)
plot(f_sat,5e-3*f_sat.^-5,'k:','linewidth',3)
hold off
colororder(cmap_binned)
colormap(cmap_binned)
clim(Ulims)
xlim([1e-2 2e1])
ylim([1e-10 1e2])
ylabel('$\mathrm{F(f)\ [m^2Hz^{-1}]}$','Interpreter','latex')
text(mean(f_eq),10.^mean(log10(1e-2*f_eq.^-4))*7,'$\mathrm{f^{-4}}$','FontSize',fsize,'HorizontalAlignment','center','Interpreter','latex')
text(mean(f_sat),10.^mean(log10(5e-3*f_sat.^-5))*7,'$\mathrm{f^{-5}}$','FontSize',fsize,'HorizontalAlignment','center','Interpreter','latex')

% Sample size per wave-age bin
ax_panel = gca;
ax_panel.XScale = 'log';
ax_panel.YScale = 'log';

bin_label_size = fsize*1.25;
x_left = 0.035;
y_base = 0.045;
lateral_gap = 0.04;

% One token per bin, in the color of its curve; the curves are drawn in
% reverse bin order, so the colormap is indexed the same way
label_string = cell(1,nU+1);
label_color = cell(1,nU+1);
label_string{1} = 'N  =';
label_color{1} = [0 0 0];
for n = 1:nU
    label_string{n+1} = num2str(n_per_bin(n));
    label_color{n+1} = cmap_binned(nU+1-n,:);
end

% Lay the tokens out left to right using throwaway text objects, whose extents
% are then reused to place the bordered text and its backing box
h_measure = gobjects(1,length(label_string));
label_extent = NaN*ones(length(label_string),4);

x_cursor = x_left;
for n = 1:length(label_string)

    h_measure(n) = text(x_cursor,y_base,label_string{n},'Units','normalized', ...
        'FontSize',bin_label_size,'FontWeight','bold', ...
        'HorizontalAlignment','left','VerticalAlignment','bottom');

    label_extent(n,:) = h_measure(n).Extent;

    x_cursor = x_cursor + label_extent(n,3) + lateral_gap;

end

token_x = label_extent(:,1)';

delete(h_measure)

% Backing box, in normalized units
box_pad = 0.014;
box_x = [min(label_extent(:,1))-box_pad x_cursor-lateral_gap+box_pad];
box_y = [min(label_extent(:,2))-box_pad max(label_extent(:,2)+label_extent(:,4))+box_pad];

% Both axes are logarithmic, so normalized positions convert through log10
x_limits = xlim;
y_limits = ylim;
box_x_data = 10.^(log10(x_limits(1)) + box_x*diff(log10(x_limits)));
box_y_data = 10.^(log10(y_limits(1)) + box_y*diff(log10(y_limits)));
token_x_data = 10.^(log10(x_limits(1)) + token_x*diff(log10(x_limits)));
token_y_data = 10.^(log10(y_limits(1)) + y_base*diff(log10(y_limits)));

hold on
patch(box_x_data([1 2 2 1]),box_y_data([1 1 2 2]),0.5*[1 1 1], ...
    'EdgeColor','k','LineWidth',0.75,'FaceAlpha',0.15);

for n = 1:length(label_string)
    textborder(token_x_data(n),token_y_data,label_string{n}, ...
        label_color{n},'k','FontSize',bin_label_size,'FontWeight','bold', ...
        'HorizontalAlignment','left','VerticalAlignment','bottom')
end
hold off

nexttile(2)
hold on
loglog(k_rad_m_combined,F_k_binned,'k-','linewidth',3)
loglog(k_rad_m_combined,F_k_binned,'linewidth',2)
plot(k_eq,4e-2*k_eq.^-2.5,'k--','linewidth',2.5)
plot(k_sat,4e-2*k_sat.^-3,'k:','linewidth',3)
hold off
colororder(cmap_binned)
xlim([1e-2 1e3])
ylim([1e-12 1e2])
ylabel('$\mathrm{F(k)\ [m^3]}$','Interpreter','latex')
text(1.5*10^mean(log10(k_eq)),10.^mean(log10(4e-2*k_eq.^-2.5))*20,'$\mathrm{k^{-2.5}}$','FontSize',fsize,'HorizontalAlignment','center','Interpreter','latex')
text(1.5*10^mean(log10(k_sat)),10.^mean(log10(4e-2*k_sat.^-3))*20,'$\mathrm{k^{-3}}$','FontSize',fsize,'HorizontalAlignment','center','Interpreter','latex')

nexttile(3)
hold on
loglog(f_Hz_combined,k_rad_m_disp.^2.*f_Hz_combined.*F_f_binned,'k-','linewidth',3)
loglog(f_Hz_combined,k_rad_m_disp.^2.*f_Hz_combined.*F_f_binned,'linewidth',2)
plot(f_eq,1.2*3e-2*f_eq.^1,'k--','linewidth',2.5)
plot(f_sat,1.2*1.5e-2*f_sat.^0,'k:','linewidth',3)
hold off
colororder(cmap_binned)
xlim([1e-2 2e1])
ylim([1e-4 5e-2])
xlabel('$\mathrm{f\ [Hz]}$','Interpreter','latex')
ylabel('$\mathrm{(2\pi f)^5g^{-2}F(f)\ [rad]}$','Interpreter','latex')
text(mean(f_eq),1.3*10.^mean(log10(3e-2*f_eq.^1))*1.5,'$\mathrm{f^{1}}$','FontSize',fsize,'HorizontalAlignment','center','Interpreter','latex')
text(mean(f_sat),1.3*10.^mean(log10(1.5e-2*f_sat.^0))*1.5,'$\mathrm{f^{0}}$','FontSize',fsize,'HorizontalAlignment','center','Interpreter','latex')

nexttile(4)
hold on
loglog(k_rad_m_combined,k_rad_m_combined.^3.*F_k_binned,'k-','linewidth',3)
loglog(k_rad_m_combined,k_rad_m_combined.^3.*F_k_binned,'linewidth',2)
plot(k_eq,1.8e-2*k_eq.^0.5,'k--','linewidth',2.5)
plot(k_sat,1.7e-2*k_sat.^0,'k:','linewidth',3)
hold off
colororder(cmap_binned)
xlim([1e-2 1e3])
ylim([1e-4 5e-2])
xlabel('$\mathrm{k\ [rad\ m^{-1}]}$','Interpreter','latex')
ylabel('$\mathrm{k^3F(k)\ [rad]}$','Interpreter','latex')
text(10^mean(log10(k_eq)),1.1*10.^mean(log10(1.8e-2*k_eq.^0.5))*1.5,'$\mathrm{k^{0.5}}$','FontSize',fsize,'HorizontalAlignment','center','Interpreter','latex')
text(10^mean(log10(k_sat)),1.1*10.^mean(log10(1.7e-2*k_sat.^0))*1.5,'$\mathrm{k^{0}}$','FontSize',fsize,'HorizontalAlignment','center','Interpreter','latex')

for n = 1:4
    ax_struc(n).ax = nexttile(n);
    ax_struc(n).ax.XScale = 'log';
    ax_struc(n).ax.YScale = 'log';
    ax_struc(n).ax.XTick = 10.^(-3:1:3);
    box on
end

ax_struc(2).ax.YAxisLocation = 'right';
ax_struc(4).ax.YAxisLocation = 'right';

for n = 1:4
    nexttile(n)
    text(text_x,text_y,labels{n},'Color',[0 0 0],'Units','normalized','FontSize',fsize,'HorizontalAlignment','center')
    clim(Ulims)
end

tlayout.TileSpacing = 'tight';

for n = 1:2
    ax_struc(n).ax.XTickLabel = '';
    ax_struc(n).ax.YTick = 10.^(-20:5:5);
end

cbar = colorbar(ax_struc(2).ax);

cbar.Layout.Tile = 'north';
cbar.Layout.TileSpan = [2 2];
cbar.Ticks = wc.ticks;
cbar.TickLabels = wc.ticklabels;
cbar.TickLabels = flipud(cbar.TickLabels);
cbar.Direction = 'reverse';
set(get(cbar,'Label'),'String','$\mathrm{c_p/u_*}$','Interpreter','LaTeX')
