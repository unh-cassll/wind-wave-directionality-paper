
function spectral_subrange_slope_extraction_example(fignum,fsize)

load('codes/coolwarm.mat')

short_wave_nc_name = 'data/ASIT2019_wave_spectra_stats_timeseries_empirical_gain.nc';
long_wave_nc_name = 'data/ASIT2019_EPSS_directional_spectra.nc';

g = 9.81;

f_Hz_Pyxis = double(ncread(short_wave_nc_name,'f_Hz'));

theta_rad_Pyxis = double(ncread(short_wave_nc_name,'theta_rad'));

SFTHETA = double(ncread(short_wave_nc_name,'S_f_theta'));

dtheta = median(diff(theta_rad_Pyxis));
Sf = squeeze(sum(SFTHETA,2,'omitnan'))*dtheta;

f_Hz_EPSS = double(ncread(long_wave_nc_name,'frequency'));
theta_rad_EPSS = double(ncread(long_wave_nc_name,'direction'))*pi/180;
F_f_m2_Hz_rad_EPSS = double(ncread(long_wave_nc_name,'F_f_d'))*180/pi;

dtheta = median(diff(theta_rad_EPSS));

Ff_EPSS = squeeze(sum(F_f_m2_Hz_rad_EPSS,2,'omitnan'))'*dtheta;

f_cut = 0.2;
inds_keep = f_Hz_EPSS < f_cut;
f_Hz_EPSS = f_Hz_EPSS(inds_keep);
Ff_EPSS = Ff_EPSS(inds_keep,:);

inds_keep_Pyxis = f_Hz_Pyxis > f_cut;
f_Hz_Pyxis_keep = f_Hz_Pyxis(inds_keep_Pyxis);
Sf_Pyxis_keep = Sf(inds_keep_Pyxis,:);

f_Hz_combined = [f_Hz_EPSS(:); f_Hz_Pyxis_keep(:)];
Ff_combined = [Ff_EPSS; ((2*pi*f_Hz_Pyxis_keep(:)).^2/g).^-2.*Sf_Pyxis_keep;];

inds_keep = f_Hz_combined <= 10;
f_Hz_combined = f_Hz_combined(inds_keep);
Ff_combined = Ff_combined(inds_keep,:);

particular_ind = 66;

f_new = logspace(log10(5e-2),log10(5e0),64)';
Ff_new = interp1(f_Hz_combined,Ff_combined,f_new);

log_f = log10(f_new);
log_Ff = log10(Ff_new);

indices_per_window = 16;

colors = flipud(spectral(7));
violet = colors(1,:);
teal = colors(2,:);

cmap = viridis(7);
cerulean = cmap(3,:);

text_x = 0.05;
text_y = 0.95;
labels = {'(a)','(b)','(c)','(d)','(e)','(f)'};


log_Ff_bit = log_Ff(:,particular_ind);

ind_p = find(log_Ff_bit==max(log_Ff_bit),1,'first');

f_slope = NaN*log_f;
R2 = f_slope;
SE = f_slope;


for i = indices_per_window/2:length(log_f)-indices_per_window/2+1

    istart = i - indices_per_window/2 + 1;
    iend = i + indices_per_window/2 - 1;

    x = log_f(istart:iend);
    y = log_Ff_bit(istart:iend);

    lm = fitlm(x,y);

    f_slope(i) = lm.Coefficients.Estimate(2);
    R2(i) = lm.Rsquared.Ordinary;
    SE(i) = lm.Coefficients.SE(2);

end

f_slope = f_slope + 5;

f_slope(f_Hz_combined<1.5*f_Hz_combined(ind_p)) = NaN;

x = log_f;
y = f_slope;
z = SE;

inds_keep = ~isnan(x) & ~isnan(y) & x < 1 & ~isinf(x) & ~isinf(y);
x = x(inds_keep);
y = y(inds_keep);
z = z(inds_keep);
ind_start = find(y == max(y));
ind_end = find(y == min(y));
x = x(ind_start:ind_end);
y = y(ind_start:ind_end);
z = z(ind_start:ind_end);

ind_eq_start = find(y-2*z>1,1,'last');
ind_eq_end = find(y-2*z>0,1,'last');
ind_sat_end = find(y+2*z>-0.25,1,'last');

f_eq_start = 10.^x(ind_eq_start);
f_eq_end = 10.^x(ind_eq_end);
f_sat_end = 10.^x(ind_sat_end);

figure(fignum);clf
tlayout = tiledlayout(2,1);
ax_struc = struct();

nexttile(1)
hold on
plot([1e-2 1e2],[1 1]*1,'r--','linewidth',2)
plot([1e-2 1e2],[1 1]*0,'r:','linewidth',2)
plot([1 1]*f_eq_start,[-1 1]*3,'k--','linewidth',2)
plot([1 1]*f_eq_end,[-1 1]*3,'k-','linewidth',2)
plot([1 1]*f_sat_end,[-1 1]*3,'k:','linewidth',3)
plot(f_new,f_slope,'.-','Color',cerulean,'markersize',15)
for i = 1:length(f_slope)
    plot(f_new(i)*[1 1],f_slope(i)*[1 1]+[-1 1]*SE(i)*2,'-','linewidth',1.5,'Color',cerulean)
end
hold off
box on
ax_struc(1).ax = gca;
ax_struc(1).ax.XScale = 'log';
ax_struc(1).ax.YTick = -3:1:3;
xlim([1e-2 1e1]*2)
ylim([-1 1]*3)
xlabel('$\mathrm{f\ [Hz]}$','Interpreter','LaTeX')
ylabel('$\mathrm{\lq\lq n"\ such\ that\ B(f)\propto f^{n}}$','Interpreter','LaTeX')

nexttile(2)
plot(f_Hz_EPSS,(2*pi*f_Hz_EPSS).^5/g^2.*Ff_EPSS(:,particular_ind),'-',f_Hz_Pyxis,((2*pi*(f_Hz_Pyxis)).^1).*Sf(:,particular_ind),'-','linewidth',3)
hold on
plot([1 1]*f_eq_start,[1e-5 1e0],'k--','linewidth',2)
plot([1 1]*f_eq_end,[1e-5 1e0],'k-','linewidth',2)
plot([1 1]*f_sat_end,[1e-5 1e0],'k:','linewidth',3)
s_combined = scatter(f_new,(2*pi*f_new).^5/g^2.*Ff_new(:,particular_ind),60,f_slope,'filled');
hold off
box on
ax_struc(2).ax = gca;
ax_struc(2).ax.XScale = 'log';
ax_struc(2).ax.YScale = 'log';
xlim([1e-2 1e1]*2)
ylim([1e-4 1e0])
clim([-1 1]*3)
colormap(coolwarm)
xlabel('$\mathrm{f\ [Hz]}$','Interpreter','LaTeX')
ylabel('$\mathrm{B(f)\ [rad]}$','Interpreter','LaTeX')
cbar = colorbar;
set(get(cbar,'Title'),'String','$\mathrm{\lq\lq n"}$','Interpreter','LaTeX')

text(0.08,0.0035,'direct','FontSize',fsize,'fontweight','bold','HorizontalAlignment','center','color',teal)
text(0.08,0.00035,'E-PSS','FontSize',fsize,'fontweight','bold','HorizontalAlignment','center','color',violet)

s_combined.Marker = 's';
s_combined.LineWidth = 1.5;

tile_cleaner(ax_struc,tlayout)

for n = 1:length(ax_struc)
    nexttile(n)
    text(text_x,text_y,labels{n},'HorizontalAlignment','center','units','normalized','fontsize',fsize)
    ax_struc(n).ax.XTick = 10.^(-3:1:3);
end
