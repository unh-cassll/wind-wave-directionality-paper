%
function wind_wave_subrange_directions(fignum,fsize)

clc

load('frequency_spect_range_limits.mat');
k_limits = load('wavenumber_spect_range_limits.mat');
k_sat_end = k_limits.k_sat_end;

short_wave_nc_name = 'data/ASIT2019_wave_spectra_stats_timeseries_empirical_gain.nc';
long_wave_nc_name = 'data/ASIT2019_EPSS_directional_spectra.nc';
supporting_nc_name = 'data/ASIT2019_supporting_environmental_observations.nc';

integrated_wave_breaking_quantities = load('data/integrated_wave_breaking_quantities.mat');

wind_driven_current_data = load('ASIT2019_skin_Doppler_subsurface_velocities.mat');

current_Winddir = wind_driven_current_data.depth_01cm_Winddir;
current_U_dir = wind_driven_current_data.depth_01cm_Udir;
current_Wavedir = wind_driven_current_data.depth_01cm_Wavedir;

f_Hz_Pyxis = double(ncread(short_wave_nc_name,'f_Hz'));
k_rad_m = double(ncread(short_wave_nc_name,'k_rad_m'));
theta_rad_Pyxis = double(ncread(short_wave_nc_name,'theta_rad'));

wind_dir_deg_coming_from = ncread(supporting_nc_name,'COARE_Wdir');
wind_dir_deg_going_to = mod(wind_dir_deg_coming_from+180,360);
EC_ustar_m_s = ncread(supporting_nc_name,'EC_ustar_m_s');

SFTHETA = double(ncread(short_wave_nc_name,'S_f_theta'));
SKTHETA = double(ncread(short_wave_nc_name,'S_k_theta'));

f_Hz_EPSS = double(ncread(long_wave_nc_name,'frequency'));
theta_rad_EPSS = double(ncread(long_wave_nc_name,'direction'))*pi/180;
FFTHETA = double(ncread(long_wave_nc_name,'F_f_d'))*180/pi;

f_Hz_EPSS = f_Hz_EPSS(2:end);
FFTHETA = FFTHETA(:,:,2:end);

FFTHETA = permute(FFTHETA,[3 2 1]);

SFTHETA(SFTHETA==0) = NaN;
FFTHETA(FFTHETA==0) = NaN;

f_cut_low = 0.75;
f_cut_high = 1.25;

FFTHETA = FFTHETA(f_Hz_EPSS<f_cut_low,:,:);
f_Hz_EPSS = f_Hz_EPSS(f_Hz_EPSS<f_cut_low);

SFTHETA = SFTHETA(f_Hz_Pyxis>f_cut_high,:,:);
f_Hz_Pyxis = f_Hz_Pyxis(f_Hz_Pyxis>f_cut_high);

water_depth_m = 15;
[c_disp,~] = lindisp_with_current(2*pi*f_Hz_Pyxis,water_depth_m,0);
k_disp_Pyxis = 2*pi*f_Hz_Pyxis./c_disp;

FFTHETA_EPSS_big = [FFTHETA FFTHETA FFTHETA];
FFTHETA_direct_big = k_disp_Pyxis.^-2.*[SFTHETA SFTHETA SFTHETA];

theta_rad_EPSS_big = [theta_rad_EPSS(:)'-2*pi theta_rad_EPSS(:)' theta_rad_EPSS(:)'+2*pi];
theta_rad_direct_big = [theta_rad_Pyxis(:)'-2*pi theta_rad_Pyxis(:)' theta_rad_Pyxis(:)'+2*pi];

inds_keep_EPSS = theta_rad_EPSS_big >= -pi & theta_rad_EPSS_big < pi;
inds_keep_direct = theta_rad_direct_big >= -pi & theta_rad_direct_big < pi;

FFTHETA_combined = [FFTHETA_EPSS_big(:,inds_keep_EPSS,:); FFTHETA_direct_big(:,inds_keep_direct,:)];

f_Hz_combined = [f_Hz_EPSS; f_Hz_Pyxis];

theta_rad = theta_rad_EPSS_big(inds_keep_EPSS);

T_E = squeeze(sum(f_Hz_combined.^-1.*FFTHETA_combined,[1 2],'omitnan'))./squeeze(sum(FFTHETA_combined,[1 2],'omitnan'));
f_E = T_E.^-1;

[c_E,~] = lindisp_with_current(2*pi*f_E,water_depth_m,0);

wave_age = c_E./EC_ustar_m_s;

Sm = squeeze(mean(sin(theta_rad(:)').*f_Hz_combined.^-1.*FFTHETA_combined,[1 2],'omitnan'));
Cm = squeeze(mean(cos(theta_rad(:)').*f_Hz_combined.^-1.*FFTHETA_combined,[1 2],'omitnan'));

Vm = atan2(Sm,Cm);

theta_E = 180/pi*Vm;

theta_eq = NaN*theta_E;
theta_sat = theta_eq;
theta_sg = theta_eq;
theta_gc = theta_eq;

inds_gc = k_rad_m > 117;

for n = 1:length(f_eq_start)

    Fftheta_slice = squeeze(FFTHETA_combined(:,:,n));

    inds_eq = f_Hz_combined >= f_eq_start(n) & f_Hz_combined <= f_eq_end(n);

    Sm = mean(sin(theta_rad(:)').*f_Hz_combined(inds_eq).^-1.*Fftheta_slice(inds_eq,:),'all','omitnan');
    Cm = mean(cos(theta_rad(:)').*f_Hz_combined(inds_eq).^-1.*Fftheta_slice(inds_eq,:),'all','omitnan');
    Vm = atan2(Sm,Cm);
    theta_eq(n) = 180/pi*Vm;

    Sktheta_slice = squeeze(SKTHETA(:,:,n));

    inds_sat = k_rad_m <= k_sat_end(n);

    Sm = mean(sin(theta_rad_Pyxis(:)').*k_rad_m(inds_sat).^-2.5.*Sktheta_slice(inds_sat,:),'all','omitnan');
    Cm = mean(cos(theta_rad_Pyxis(:)').*k_rad_m(inds_sat).^-2.5.*Sktheta_slice(inds_sat,:),'all','omitnan');
    Vm = atan2(Sm,Cm);
    theta_sat(n) = 180/pi*Vm;

    inds_sg = k_rad_m > k_sat_end(n) & k_rad_m < 117;

    Sm = mean(sin(theta_rad_Pyxis(:)').*k_rad_m(inds_sg).^-2.5.*Sktheta_slice(inds_sg,:),'all','omitnan');
    Cm = mean(cos(theta_rad_Pyxis(:)').*k_rad_m(inds_sg).^-2.5.*Sktheta_slice(inds_sg,:),'all','omitnan');
    Vm = atan2(Sm,Cm);
    theta_sg(n) = 180/pi*Vm;

    Sktheta_big = [Sktheta_slice Sktheta_slice Sktheta_slice];
    theta_rad_big = [theta_rad_Pyxis(:)'-2*pi theta_rad_Pyxis(:)' theta_rad_Pyxis(:)'+2*pi] - Vm;
    mask = NaN*Sktheta_big;
    mask(:,theta_rad_big>=-pi/2&theta_rad_big<pi/2) = 1;

    Sm = mean(mask(inds_gc,:).*sin(theta_rad_big(:)').*k_rad_m(inds_gc).^-2.5.*Sktheta_big(inds_gc,:),'all','omitnan');
    Cm = mean(mask(inds_gc,:).*cos(theta_rad_big(:)').*k_rad_m(inds_gc).^-2.5.*Sktheta_big(inds_gc,:),'all','omitnan');
    Vm = atan2(Sm,Cm);
    theta_gc(n) = 180/pi*Vm + theta_sg(n);

end

% Remove duplicate points
theta_E(136:141) = NaN;

theta_E = mod(theta_E+360,360);
theta_eq = mod(theta_eq+360,360);
theta_sat = mod(theta_sat+360,360);
theta_sg = mod(theta_sg+360,360);
theta_gc = mod(theta_gc+360,360);

% Remove runs for which wind is not coming towards the anemometer (+/-90 degrees)
diving_board_heading = 255;
winddir_rel_diving_board = compute_relative_angle(wind_dir_deg_coming_from*0+diving_board_heading,wind_dir_deg_coming_from)';
bad_winddir_inds = abs(winddir_rel_diving_board) > 90;
wind_dir_deg_going_to(bad_winddir_inds) = NaN;

winddir_rel_diving_board = compute_relative_angle(integrated_wave_breaking_quantities.wdir_deg*0+diving_board_heading,integrated_wave_breaking_quantities.wdir_deg+180)';
bad_winddir_inds = abs(winddir_rel_diving_board) > 90;
integrated_wave_breaking_quantities.wdir_deg(bad_winddir_inds) = NaN;

% Remove runs for which breaking crests are moving offshore (+/-90 degrees)
integrated_wave_breaking_quantities.theta_br(integrated_wave_breaking_quantities.theta_br<270 & integrated_wave_breaking_quantities.theta_br > 90) = NaN;

dir_diff_wind_De = compute_relative_angle(theta_E,wind_dir_deg_going_to);
dir_diff_wind_Dslope_eq = compute_relative_angle(theta_eq,wind_dir_deg_going_to);
dir_diff_wind_Dslope_sat = compute_relative_angle(theta_sat,wind_dir_deg_going_to);
dir_diff_wind_Dslope_sg = compute_relative_angle(theta_sg,wind_dir_deg_going_to);
dir_diff_wind_Dslope_gc = compute_relative_angle(theta_gc,wind_dir_deg_going_to);

dir_diff_wind_D_ds = compute_relative_angle(integrated_wave_breaking_quantities.theta_br,integrated_wave_breaking_quantities.wdir_deg);

t = table(dir_diff_wind_De(:),dir_diff_wind_Dslope_eq(:),dir_diff_wind_Dslope_sat(:),dir_diff_wind_Dslope_sg(:),dir_diff_wind_Dslope_gc(:));
t.Properties.VariableNames = {'\theta_{E}-\theta_{wind} [\circ]','\theta_{eq.}-\theta_{wind} [\circ]','\theta_{sat.}-\theta_{wind} [\circ]','\theta_{s.g.}-\theta_{wind} [\circ]','\theta_{g-c}-\theta_{wind} [\circ]'};

inds_remove = isnan(t.(1));
for n = 1:size(t,2)
    holder = t.(n);
    holder(inds_remove) = NaN;
    t.(n) = holder;
end

% Plot

labels = {'(a)','(b)','(c)','(d)','(e)','(f)'};
label_x = 0.015;
label_y = 0.93;

names = {'equilibrium','saturation','short gravity','gravity-capillary','breaking','current'};

grayish = 0.5*[1 1 1];

dir_max = 180;

dD = 15;

D_ref = t.(1);
D_ref_masked = D_ref;
D_ref_masked(D_ref<-dir_max) = NaN;
D_ref_masked(D_ref>dir_max) = NaN;

holder_struc = struct();
for n = 1:size(t,2)

    D_particular = t.(n);
    holder_struc(n).D_ref = D_ref_masked;
    holder_struc(n).D_particular = D_particular;
    holder_struc(n).wave_age = wave_age;

end

dir_diff_wind_De_breakers = compute_relative_angle(integrated_wave_breaking_quantities.D_E_deg,integrated_wave_breaking_quantities.wdir_deg);
D_ref_breakers = dir_diff_wind_De_breakers;
D_particular_breakers = dir_diff_wind_D_ds;

holder_struc(n+1).D_ref = D_ref_breakers;
holder_struc(n+1).D_particular = D_particular_breakers;
holder_struc(n+1).wave_age = integrated_wave_breaking_quantities.c_E_m_s./integrated_wave_breaking_quantities.ustar_m_s;

D_ref_current = compute_relative_angle(current_Wavedir,current_Winddir);
D_particular_current = compute_relative_angle(current_U_dir,current_Winddir);

holder_struc(n+2).D_ref = D_ref_current;
holder_struc(n+2).D_particular = D_particular_current;
holder_struc(n+2).wave_age = wave_age;

% order:
% mean
% equilibrium
% saturation
% short gravity
% gravity-capillary
% breakers
% current

msize = 6;
lw = 0.5;

fA = 0.2;

wave_age_lims = [10 70];

figure(fignum);clf
tlayout = tiledlayout(3,2);
ax_struc = struct();

stats_holder = struct();

for n = 1:6

    nexttile(n)
    hold on
    plot([-180 180],[0 0],'--','Color',0.5*[1 1 1],'linewidth',2)
    plot([0 0],[-180 180],'--','Color',0.5*[1 1 1],'linewidth',2)

    x = holder_struc(n+1).D_ref(:);
    y = holder_struc(n+1).D_particular(:);
    c = holder_struc(n+1).wave_age(:);
    inds_keep = ~isnan(x) & ~isnan(y) & ~isnan(c);
    x = x(inds_keep);
    y = y(inds_keep);
    c = c(inds_keep);
    xplot = linspace(-90+dD/2,90-dD/2,100);

    inds_fit = x>-90+dD/2 & x<90-dD/2;
    D_fit_object = fitlm(x(inds_fit),y(inds_fit));
    [D_fit,D_fit_95CI] = predict(D_fit_object,xplot(:));
    D_lower = D_fit_95CI(:,1)';
    D_upper = D_fit_95CI(:,2)';

    stats_holder(n).line_slope = D_fit_object.Coefficients.Estimate(2);
    stats_holder(n).line_intercept = D_fit_object.Coefficients.Estimate(1);
    stats_holder(n).line_slope_pValue = D_fit_object.Coefficients.pValue(2);
    stats_holder(n).line_intercept_pValue = D_fit_object.Coefficients.pValue(1);
    stats_holder(n).R2 = D_fit_object.Rsquared.Adjusted;

    plot(x,y,'o','markerfacecolor','w','markeredgecolor','k','markersize',msize,'linewidth',lw)
    scatter(x,y,0.8*msize^2,c,'filled')
    colormap(flipud(magma))
    clim(wave_age_lims)

    if n == 1
        cbar = colorbar;
    end

    plot(xplot,D_fit,'-','Color','k','linewidth',2)
    f = fill([xplot fliplr(xplot)],[D_upper fliplr(D_lower)],grayish);
    f.FaceAlpha = fA;

    hold off
    box on
    pbaspect([1 0.5 1])

    xlim([-1 1]*dir_max)
    ylim([-1 1]*dir_max/2)
    ax_struc(n).ax = gca;
    ax_struc(n).ax.YTick = -360:45:360;
    ax_struc(n).ax.XTick = -360:45:360;
    ax_struc(n).ax.YDir = 'reverse';

    xlabel('\theta_E-\theta_{wind} [\circ]')
    ylabel('\theta_{particular}-\theta_{wind} [\circ]')

end

tile_cleaner(ax_struc,tlayout)

for n = 1:length(ax_struc)
    nexttile(n)
    text(label_x,label_y,labels{n},'Units','normalized','FontSize',fsize)
    text(1-label_x, label_y, names{n},'color', [38 38 38]/255,'FontSize',fsize*1.25,'Units','Normalized','HorizontalAlignment','right')
end

tlayout.TileSpacing = 'tight';

cbar.Layout.Tile = 'north';
set(get(cbar,'Label'),'String','c_E/u_*')

stats_table = struct2table(stats_holder);
stats_table.Properties.VariableNames = {'Line_Slope', 'Line_Intercept', 'Line_Slope_pValue', 'Line_Intercept_pValue', 'R2'};
row_names = {'Equilibrium', 'Saturation', 'Short Gravity', 'Gravity-Capillary', 'Breakers', 'Current'};
stats_table.Properties.RowNames = row_names;
disp(stats_table);

% writetable(stats_table,'data/off_wind_wave_angle_table.csv','WriteRowNames',true)
