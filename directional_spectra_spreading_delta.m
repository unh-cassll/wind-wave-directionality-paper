%
function directional_spectra_spreading_delta(fignum,fsize)

supporting_data_nc_name = 'data/ASIT2019_supporting_environmental_observations.nc';
short_wave_nc_name = 'ASIT2019_wave_spectra_stats_timeseries_empirical_gain.nc';

load('data/directional_spreading_quantification_data.mat')

s = load('data/global_figure_settings.mat');
example_run_ind = s.example_run_ind;

cmap = flipud(spectral(7));
crimson = cmap(7,:);

g = 9.81;

wdir_deg = ncread(supporting_data_nc_name,'COARE_Wdir');
wdir_deg = wdir_deg(example_run_ind);
wdir_deg = mod(wdir_deg+180,360);
wdir_rad = pi/180*wdir_deg;

f_Hz_Pyxis = ncread(short_wave_nc_name,'f_Hz');
k_rad_m_Pyxis = ncread(short_wave_nc_name,'k_rad_m');
theta_rad_Pyxis = ncread(short_wave_nc_name,'theta_rad');

S_f_theta = double(ncread(short_wave_nc_name,'S_f_theta'));
S_k_theta = double(ncread(short_wave_nc_name,'S_k_theta'));

S_f_theta = squeeze(S_f_theta(:,:,example_run_ind));
S_k_theta = squeeze(S_k_theta(:,:,example_run_ind));

[S_f_theta,theta_rad] = convert_dirspect_to_downwind(S_f_theta,theta_rad_Pyxis,wdir_rad);
[S_k_theta,~] = convert_dirspect_to_downwind(S_k_theta,theta_rad_Pyxis,wdir_rad);

f_cut = 0.6;

S_f_theta(f_Hz_Pyxis<f_cut,:) = NaN;

dtheta = median(diff(theta_rad_Pyxis));
theta_rad = theta_rad + dtheta/2;

s = load('data/global_figure_settings.mat');
klow = s.k_low;
khigh = s.k_high;
fhigh = s.f_high;
flow = sqrt(9.81*klow)/(2*pi);

S_k_theta(k_rad_m_Pyxis<klow,:) = NaN;
S_f_theta(f_Hz_Pyxis<flow,:) = NaN;

% Compute D
D_f_theta = smoothdata2(S_f_theta,'movmedian',{5,3},'omitnan')./max(smoothdata2(S_f_theta,'movmedian',{5,3},'omitnan'),[],2,'omitnan');
D_k_theta = smoothdata2(k_rad_m_Pyxis.*S_k_theta,'movmedian',{3,3},'omitnan')./max(smoothdata2(k_rad_m_Pyxis.*S_k_theta,'movmedian',{3,3},'omitnan'),[],2,'omitnan');

dir_ticks = 180*(-1:0.25:1);

S_k_theta(S_k_theta==0) = NaN;
BKTHETA = k_rad_m_Pyxis.^2.*S_k_theta;
BFTHETA = (2*pi*f_Hz_Pyxis).^2/g.*S_f_theta;

BKTHETA_smooth = smoothdata2(BKTHETA,'movmean',{3,3});
BFTHETA_smooth = smoothdata2(BFTHETA,'movmean',{5,3});

BKTHETA = (BKTHETA-min(BKTHETA,[],'all','omitnan'))/(max(BKTHETA,[],'all','omitnan')-min(BKTHETA,[],'all','omitnan'));
BFTHETA = (BFTHETA-min(BFTHETA,[],'all','omitnan'))/(max(BFTHETA,[],'all','omitnan')-min(BFTHETA,[],'all','omitnan'));

Dk50_left = smoothdata(squeeze(D_spread_holder_struc(example_run_ind).D_k_limits(:,1,1)),'movmean',5)-10*pi/180;
Dk50_right = smoothdata(squeeze(D_spread_holder_struc(example_run_ind).D_k_limits(:,2,1)),'movmean',5)-10*pi/180;

Dk50_left(k_rad_m_Pyxis<klow) = NaN;
Dk50_right(k_rad_m_Pyxis<klow) = NaN;

Df50_left = smoothdata(squeeze(D_spread_holder_struc(example_run_ind).D_f_limits(:,1,1)),'movmean',11)-10*pi/180;
Df50_right = smoothdata(squeeze(D_spread_holder_struc(example_run_ind).D_f_limits(:,2,1)),'movmean',11)-10*pi/180;

Df50_left(f_Hz_Pyxis<flow) = NaN;
Df50_right(f_Hz_Pyxis<flow) = NaN;

inds_0_deg = 35:37;
inds_90_deg = [17:19 53:55];

Delta_k = (mean(BKTHETA_smooth(:,inds_0_deg),2,'omitnan')-mean(BKTHETA_smooth(:,inds_90_deg),2,'omitnan'))./(mean(BKTHETA_smooth(:,inds_0_deg),2,'omitnan')+mean(BKTHETA_smooth(:,inds_90_deg),2,'omitnan'));
Delta_f = (mean(BFTHETA_smooth(:,inds_0_deg),2,'omitnan')-mean(BFTHETA_smooth(:,inds_90_deg),2,'omitnan'))./(mean(BFTHETA_smooth(:,inds_0_deg),2,'omitnan')+mean(BFTHETA_smooth(:,inds_90_deg),2,'omitnan'));

diff_k = khigh/5;
diff_f = fhigh/5;
k_ticks = 0:diff_k:khigh;
f_ticks = 0:diff_f:fhigh;

figure(fignum);clf

ax_struc = struct();

tlayout = tiledlayout(2,3);

ax_struc(1).ax = nexttile(1);
pcolor((theta_rad-dtheta/2)*180/pi,k_rad_m_Pyxis,log10(BKTHETA))
cbar_spect = colorbar;
cbar_spect.Location = 'northoutside';
shading('flat')
xlabel('\theta [\circ]')
ylabel('k [rad m^{-1}]')

ax_struc(2).ax = nexttile(2);
hold on
pcolor((theta_rad-dtheta/2)*180/pi,k_rad_m_Pyxis,D_k_theta)
plot(180/pi*Dk50_left,k_rad_m_Pyxis,':','Color','k','linewidth',2)
plot(180/pi*Dk50_right,k_rad_m_Pyxis,':','Color','k','linewidth',2)
hold off
box on
cbar_spread = colorbar;
cbar_spread.Location = 'northoutside';
shading('flat')
xlabel('\theta [\circ]')
ylabel('k [rad m^{-1}]')

ax_struc(3).ax = nexttile(3);
plot(Delta_k,k_rad_m_Pyxis,'linewidth',3,'Color',crimson)
xlabel('\Delta(k) [rad]')
ylabel('k [rad m^{-1}]')

ax_struc(4).ax = nexttile(4);
pcolor((theta_rad-dtheta/2)*180/pi,f_Hz_Pyxis,log10(BFTHETA))
shading('flat')
xlabel('\theta [\circ]')
ylabel('f [Hz]')

ax_struc(5).ax = nexttile(5);
hold on
pcolor((theta_rad-dtheta/2)*180/pi,f_Hz_Pyxis,D_f_theta)
plot(180/pi*Df50_left,f_Hz_Pyxis,':','Color','k','linewidth',2)
plot(180/pi*Df50_right,f_Hz_Pyxis,':','Color','k','linewidth',2)
hold off
box on
shading('flat')
xlabel('\theta [\circ]')
ylabel('f [Hz]')

ax_struc(6).ax = nexttile(6);
plot(Delta_f,f_Hz_Pyxis,'linewidth',3,'Color',crimson)
xlabel('\Delta [rad]')
ylabel('f [Hz]')

for n = 1:3
    nexttile(n)
    ylim([0 khigh])
    ax_struc(n).ax.YTick = k_ticks;
end

for n = 4:6
    nexttile(n)
    ylim([0 fhigh])
    ax_struc(n).ax.YTick = f_ticks;
end

for n = [2 5]
    colormap(ax_struc(n).ax,flipud(spectral))
    nexttile(n)
    clim([0 1])
end

for n = [1 4]
    nexttile(n)
    clim([-2 0])
end

cbar_spect.Ticks = [-2 -1 0];
cbar_spect.TickLabels = {'1%','10%','100%'};
% cbar_spect.Position = [0.12 0.86 0.235 0.0258];
set(get(cbar_spect,'label'),'String','$B^*(X,\theta)$','Interpreter','latex')

cbar_spread.Ticks = [0.25 0.5 0.75];
cbar_spread.TickLabels = {'25%','50%','75%'};
set(get(cbar_spread,'label'),'String','$D(X,\theta)$','Interpreter','LaTeX')

for n = [1 2 4 5]
    nexttile(n)
    xlim([-1 1]*180)
    hold on
    plot([1 1]*0,[-1 1]*1e10,'--','Color',0.5*[1 1 1],'linewidth',2)
    hold off
    ax_struc(n).ax.XTick = dir_ticks;
end

for n = [3 6]
    nexttile(n)
    xlim([-1 1]*0.8)
    hold on
    plot([1 1]*0,[-1 1]*1e10,'--','Color',0.5*[1 1 1],'linewidth',2)
    hold off
end

text_x = 0.07;
text_y = 0.95;
labels = {'(a)','(b)','(c)','(d)','(e)','(f)'};

fill_vals = [1 1 0 1 1 0];
for n = 1:6
    nexttile(n)
    text(text_x,text_y,labels{n},'Color',fill_vals(n)*[1 1 1],'Units','normalized','FontSize',fsize,'HorizontalAlignment','center')
end

tile_cleaner(ax_struc,tlayout)

tlayout.TileSpacing = 'tight';