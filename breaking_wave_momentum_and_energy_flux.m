%
%
function breaking_wave_momentum_and_energy_flux(fignum,fsize)

load('codes/coolwarm.mat')
load('data/JGRO_LIU25_Data_Compiled.mat')
integrated_wave_breaking_quantities = load('data/integrated_wave_breaking_quantities.mat');

text_x = 0.03;
text_y = 0.95;
labels = {'(a)','(b)','(c)','(d)'};

U10_br = integrated_wave_breaking_quantities.U10_m_s;
ustar_br = integrated_wave_breaking_quantities.ustar_m_s;
theta_m_br = integrated_wave_breaking_quantities.D_m_deg;
c_p_br = integrated_wave_breaking_quantities.c_p_m_s;
Momentum_flux_br = integrated_wave_breaking_quantities.Momentum_flux_ds_int;
theta_br = integrated_wave_breaking_quantities.theta_br;
Wdir_br = integrated_wave_breaking_quantities.wdir_deg;

Wdir_br(Wdir_br<270&Wdir_br>90) = NaN;

U10_br(U10_br<0) = NaN;

Momentum_flux_br(Momentum_flux_br>0.5) = NaN;

off_wind_theta_br = compute_relative_angle(theta_br,Wdir_br);
off_wind_theta_m = compute_relative_angle(theta_m_br,Wdir_br);

% stress

viscosity_dynamic = 1e-3;    % Pa*s
rho_a = 1.25;

shear_fit_toskin_5_8_12 = struc_VisStress.shear_fit_toskin_5_8_12  ;
U10i = struc_VisStress.U10i  ;

wave_age_lims = [10 60];
d_wave_age = 10;
wave_age_centers = wave_age_lims(1)+d_wave_age/2:d_wave_age:wave_age_lims(2)-d_wave_age/2;

off_wind_theta_br(isnan(off_wind_theta_m)) = NaN;
Momentum_flux_br(isnan(off_wind_theta_m)) = NaN;

binned_theta_m_minus_theta_wind = NaN*ones(length(wave_age_centers),1);
binned_theta_wave_minus_theta_wind = binned_theta_m_minus_theta_wind;
binned_momentum_flux = binned_theta_m_minus_theta_wind;
binned_wave_age = binned_theta_m_minus_theta_wind;
binned_binsize = binned_theta_m_minus_theta_wind;
binned_ustar = binned_theta_m_minus_theta_wind;
binned_u10 = binned_theta_m_minus_theta_wind;

wave_age = c_p_br./ustar_br;

x = Momentum_flux_br(:);
y = off_wind_theta_br(:);
z = wave_age(:);
q = ustar_br(:);
u10 = U10_br(:);

inds_keep = ~isnan(x) & ~isnan(y) & ~isnan(z);
x = x(inds_keep);
y = y(inds_keep);
z = z(inds_keep);
q = q(inds_keep);
u10 = u10(inds_keep);

for m = 1:length(wave_age_centers)

    wave_age_low = wave_age_centers(m) - d_wave_age/2;
    wave_age_high = wave_age_centers(m) + d_wave_age/2;
    inds_consider = z >= wave_age_low & z < wave_age_high;

    binned_momentum_flux(m) = median(x(inds_consider),'all','omitnan');
    binned_theta_wave_minus_theta_wind(m) = median(y(inds_consider),'all','omitnan');
    binned_wave_age(m) = median(z(inds_consider),'all','omitnan');
    binned_binsize(m) = length(z(inds_consider));
    binned_ustar(m) = median(q(inds_consider),'all','omitnan');
    binned_u10(m) = median(u10(inds_consider),'all','omitnan');

end

binned_viscous_stress = interp1(U10i,viscosity_dynamic*shear_fit_toskin_5_8_12(:,2),binned_u10,'pchip');

msize = 10;

figure(fignum);clf
tlayout = tiledlayout(2,1);
ax_struc = struct();

nexttile()
hold on
plot([0 100],[0 0],'--','linewidth',2,'Color',0.5*[1 1 1])
plot(wave_age,off_wind_theta_br,'o','markerfacecolor','k','markeredgecolor','k','markersize',msize,'linewidth',0.5)
scatter(wave_age,off_wind_theta_br,0.65*msize^2,off_wind_theta_m,'filled')
plot(binned_wave_age,binned_theta_wave_minus_theta_wind,'o-','color','k','markerfacecolor','k','markeredgecolor','k','markersize',msize*2.25,'linewidth',2)
for i = 1:length(binned_binsize)
    text(binned_wave_age(i),binned_theta_wave_minus_theta_wind(i),num2str(binned_binsize(i)),'Color','w','FontSize',fsize,'FontWeight','bold','HorizontalAlignment','center')
end
hold off
box on
ylim([-90 90]/2)

cbar = colorbar;
cbar.Layout.Tile = 'north';

ax_struc(1).ax = gca;
ax_struc(1).ax.YTick = -45:15:45;
ax_struc(1).ax.YDir = 'reverse';

cbar.Ticks = -90:30:90;
set(get(cbar,'Label'),'String','$\mathrm{\theta_m-\theta_{wind}\ [^\circ]}$','Interpreter','LaTeX')

xlabel('$\mathrm{c_p/u_*}$','Interpreter','latex')
ylabel('$\mathrm{\theta_{br}-\theta_{wind}\ [^\circ]}$','Interpreter','latex')

nexttile()
hold on
plot([0 100],[0 0],'--','linewidth',2,'Color',0.5*[1 1 1])
plot(wave_age,Momentum_flux_br,'o','markerfacecolor','k','markeredgecolor','k','markersize',msize,'linewidth',0.5)
scatter(wave_age,Momentum_flux_br,0.65*msize^2,off_wind_theta_m,'filled')
plot(binned_wave_age,binned_momentum_flux,'o-','color','k','markerfacecolor','k','markeredgecolor','k','markersize',msize*1.25,'linewidth',2)
plot(binned_wave_age,rho_a*binned_ustar.^2,'s:','color','k','markerfacecolor','none','markeredgecolor','k','markersize',msize*1.75,'linewidth',2.5)
plot(binned_wave_age,binned_viscous_stress,'d--','color','k','markerfacecolor','none','markeredgecolor','k','markersize',msize*1.75,'linewidth',2.5)
hold off
box on
ylim([1e-3 1e1])

% text(binned_wave_age(end)+3,binned_momentum_flux(end),'$\tau_{\mathrm{br}}$','FontSize',fsize*1.25,'Interpreter','latex')
% text(binned_wave_age(end)+3,binned_viscous_stress(end),'$\tau_{\nu}$','FontSize',fsize*1.25,'Interpreter','latex')
% text(binned_wave_age(end)+3,1.25*rho_a*binned_ustar(end).^2,'$\tau_{\mathrm{total}}$','FontSize',fsize*1.25,'Interpreter','latex')

ax_struc(2).ax = gca;
ax_struc(2).ax.YScale = 'log';

xlabel('$\mathrm{c_p/u_*}$','Interpreter','latex')
ylabel('$\mathrm{\tau\ [N\ m^{-2}]}$','Interpreter','latex')

for n = 1:2
    nexttile(n)
    colormap(coolwarm)
    clim([-90 90])
    xlim([0 80])
    text(text_x,text_y,labels{n},'Units','normalized','FontSize',fsize,'HorizontalAlignment','left')
end

tile_cleaner(ax_struc,tlayout)
tlayout.TileSpacing = 'tight';