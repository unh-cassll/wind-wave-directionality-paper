%
function directional_spectra_spreading_revisit(fignum,fsize)

supporting_data_nc_name = 'data/ASIT2019_supporting_environmental_observations.nc';
short_wave_nc_name = 'ASIT2019_wave_spectra_stats_timeseries_empirical_gain.nc';

s = load('data/global_figure_settings.mat');
example_run_ind = s.example_run_ind;

g = 9.81;

wdir_deg = ncread(supporting_data_nc_name,'COARE_Wdir');
wdir_deg = wdir_deg(example_run_ind);
wdir_deg = mod(wdir_deg+180,360);
wdir_rad = pi/180*wdir_deg;

k_rad_m_Pyxis = ncread(short_wave_nc_name,'k_rad_m');
theta_rad_Pyxis = ncread(short_wave_nc_name,'theta_rad');

S_k_theta = double(ncread(short_wave_nc_name,'S_k_theta'));

S_k_theta = squeeze(S_k_theta(:,:,example_run_ind));

[S_k_theta,theta_rad] = convert_dirspect_to_downwind(S_k_theta,theta_rad_Pyxis,wdir_rad);

dir_ticks = 180*(-1:0.25:1);

S_k_theta(S_k_theta==0) = NaN;
BKTHETA = k_rad_m_Pyxis.^2.*S_k_theta;

BKTHETA = (BKTHETA-min(BKTHETA,[],'all','omitnan'))/(max(BKTHETA,[],'all','omitnan')-min(BKTHETA,[],'all','omitnan'));

D_k_theta = smoothdata2(k_rad_m_Pyxis.*S_k_theta,'movmedian',{3,3},'omitnan')./max(smoothdata2(k_rad_m_Pyxis.*S_k_theta,'movmedian',{3,3},'omitnan'),[],2,'omitnan');

s = load('data/global_figure_settings.mat');
khigh = s.k_high;
klow = 6.48;

diff_k = khigh/5;
k_ticks = 0:diff_k:khigh;

fA = 0.6;

figure(fignum);clf

ax_struc = struct();

tlayout = tiledlayout(3,1);

ax_struc(1).ax = nexttile(1);
hold on
pcolor(theta_rad*180/pi,k_rad_m_Pyxis,log10(BKTHETA))
f = fill([-180 180 180 -180],[k_rad_m_Pyxis(1) k_rad_m_Pyxis(1) klow klow],'w');
f.FaceAlpha = fA;
f.LineStyle = 'none';
hold off
box on
ax = gca;
ax.YScale = 'log';
cbar_spect = colorbar;
cbar_spect.Location = 'eastoutside';
shading('flat')
ylabel('$\mathrm{\theta\ [^\circ]}$','Interpreter','latex')
ylabel('$\mathrm{k\ [rad\ m^{-1}]}$','Interpreter','latex')

ax_struc(2).ax = nexttile(2);
hold on
pcolor(theta_rad*180/pi,k_rad_m_Pyxis,D_k_theta)
f = fill([-180 180 180 -180],[k_rad_m_Pyxis(1) k_rad_m_Pyxis(1) klow klow],'w');
f.FaceAlpha = fA;
f.LineStyle = 'none';
hold off
box on
ax = gca;
ax.YScale = 'log';
shading('flat')
ylabel('$\mathrm{\theta\ [^\circ]}$','Interpreter','latex')
ylabel('$\mathrm{k\ [rad\ m^{-1}]}$','Interpreter','latex')

ax_struc(3).ax = nexttile(3);
hold on
pcolor(theta_rad*180/pi,k_rad_m_Pyxis,D_k_theta)
f = fill([-180 180 180 -180],[k_rad_m_Pyxis(1) k_rad_m_Pyxis(1) klow klow],'w');
f.FaceAlpha = fA;
f.LineStyle = 'none';
hold off
box on
cbar_spread = colorbar;
shading('flat')
ylabel('$\mathrm{\theta\ [^\circ]}$','Interpreter','latex')
ylabel('$\mathrm{k\ [rad\ m^{-1}]}$','Interpreter','latex')

khigh = s.k_high;

for n = 1:2
    nexttile(n)
    ylim([1e0 1e3])
    ax_struc(n).YScale = 'log';
    ax_struc(n).ax.YTick = 10.^(0:1:3);
end

for n = 3
    nexttile(n)
    ylim([0 khigh])
    ax_struc(n).ax.YTick = k_ticks;
end

for n = [2 3]
    colormap(ax_struc(n).ax,flipud(spectral))
    nexttile(n)
    clim([0.01 1])
end

for n = 1
    nexttile(n)
    clim([-2 0])
end

cbar_spect.Ticks = [-2 -1 0];
cbar_spect.TickLabels = {'1%','10%','100%'};
set(get(cbar_spect,'label'),'String','$B^*(X,\theta)$','Interpreter','latex')

cbar_spread.Ticks = [0.01 0.25 0.5 0.75 1];
cbar_spread.TickLabels = {'1%','25%','50%','75%','100%'};
set(get(cbar_spread,'label'),'String','$D(X,\theta)$','Interpreter','LaTeX')

for n = 1:3
    nexttile(n)
    hold on
    plot([0 0],[1e-1 1e3],'--','Color',0.5*[1 1 1],'linewidth',2)
    plot([-1 1]*180,khigh*[1 1],'k-','linewidth',4)
    plot([-1 1]*180,khigh*[1 1],'w:','linewidth',3)
    hold off
end

for n = 1:3
    nexttile(n)
    xlim([-1 1]*180)
    ax_struc(n).ax.XTick = dir_ticks;
end

text_x = 0.04;
text_y = 1.05;
labels = {'(a)','(b)','(c)','(d)','(e)','(f)'};

for n = 1:3
    nexttile(n)
    text(text_x,text_y,labels{n},'Units','normalized','FontSize',fsize,'HorizontalAlignment','center')
end

cbar_spread.Layout.Tile = 2;
cbar_spread.Layout.TileSpan = [2 2];

tile_cleaner(ax_struc,tlayout)
