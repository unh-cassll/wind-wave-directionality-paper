%
%
function breaking_wave_momentum_and_energy_flux(fignum,fsize)

load('codes/coolwarm.mat')
load('data/JGRO_LIU25_Data_Compiled.mat')
integrated_wave_breaking_quantities = load('data/integrated_wave_breaking_quantities.mat');

s = load('data/global_figure_settings.mat');
wave_age_lims = s.wave_age_lims;
wc = wave_age_color();

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

% Oncoming-wind gate; Wdir_br is the direction the wind blows towards
% Breaking records carry bulk COARE wind, so they take the COARE half-angle
gset = load('data/global_figure_settings.mat');
Wdir_br = oncoming_wind_mask(mod(Wdir_br+180,360),gset.ASIT_COARE_wind_gate_deg).*Wdir_br;

U10_br(U10_br<0) = NaN;

Momentum_flux_br(Momentum_flux_br>0.5) = NaN;

off_wind_theta_br = compute_relative_angle(theta_br,Wdir_br);
off_wind_theta_m = compute_relative_angle(theta_m_br,Wdir_br);

% stress

viscosity_dynamic = 1e-3;    % Pa*s
rho_a = 1.25;

shear_fit_toskin_5_8_12 = struc_VisStress.shear_fit_toskin_5_8_12  ;
U10i = struc_VisStress.U10i  ;

% Wave-age binning from the shared scheme
[wave_age_edges,wave_age_centers] = wave_age_bin_edges();
n_wave_age_bins = numel(wave_age_centers);

off_wind_theta_br(isnan(off_wind_theta_m)) = NaN;
Momentum_flux_br(isnan(off_wind_theta_m)) = NaN;

binned_theta_m_minus_theta_wind = NaN*ones(n_wave_age_bins,1);
binned_theta_wave_minus_theta_wind = binned_theta_m_minus_theta_wind;
binned_momentum_flux = binned_theta_m_minus_theta_wind;
binned_wave_age = binned_theta_m_minus_theta_wind;
binned_binsize = binned_theta_m_minus_theta_wind;
binned_ustar = binned_theta_m_minus_theta_wind;
binned_u10 = binned_theta_m_minus_theta_wind;

% Interquartile range within each bin, drawn as a whisker on its own marker.
% The wave-age range is the horizontal one, and applies to both panels
binned_theta_iqr = NaN*ones(n_wave_age_bins,2);
binned_momentum_flux_iqr = binned_theta_iqr;
binned_wave_age_iqr = binned_theta_iqr;

wave_age = c_p_br./ustar_br;

x = Momentum_flux_br(:);
y = off_wind_theta_br(:);
z = wave_age(:);
q = ustar_br(:);
u10 = U10_br(:);

% Cap applied to the binning sample, not just the axis, so the binned momentum,
% total and viscous stresses all terminate inside the plot box. The excluded old
% seas are also the runs where the Lambda(c) integral goes negative through the
% breaking strength b, which is unplottable on the log axis.
inds_keep = ~isnan(x) & ~isnan(y) & ~isnan(z) & ...
    z >= wave_age_lims(1) & z <= wave_age_lims(2);
x = x(inds_keep);
y = y(inds_keep);
z = z(inds_keep);
q = q(inds_keep);
u10 = u10(inds_keep);

for m = 1:n_wave_age_bins

    wave_age_low = wave_age_edges(m);
    wave_age_high = wave_age_edges(m+1);
    if m < n_wave_age_bins
        inds_consider = z >= wave_age_low & z < wave_age_high;
    else
        inds_consider = z >= wave_age_low & z <= wave_age_high;
    end

    binned_momentum_flux(m) = median(x(inds_consider),'all','omitnan');
    binned_theta_wave_minus_theta_wind(m) = median(y(inds_consider),'all','omitnan');
    binned_wave_age(m) = median(z(inds_consider),'all','omitnan');
    binned_binsize(m) = length(z(inds_consider));
    binned_ustar(m) = median(q(inds_consider),'all','omitnan');
    binned_u10(m) = median(u10(inds_consider),'all','omitnan');

    if sum(inds_consider) > 1
        binned_theta_iqr(m,:) = prctile(y(inds_consider),[25 75]);
        binned_momentum_flux_iqr(m,:) = prctile(x(inds_consider),[25 75]);
        binned_wave_age_iqr(m,:) = prctile(z(inds_consider),[25 75]);
    end

end

binned_viscous_stress = interp1(U10i,viscosity_dynamic*shear_fit_toskin_5_8_12(:,2),binned_u10,'pchip');

msize = 10;

figure(fignum);clf
tlayout = tiledlayout(2,1);
ax_struc = struct();

nexttile()
hold on
plot(wave_age_lims,[0 0],'--','linewidth',2,'Color',0.5*[1 1 1])
plot(wave_age,off_wind_theta_br,'o','markerfacecolor','k','markeredgecolor','k','markersize',msize,'linewidth',0.5)
scatter(wave_age,off_wind_theta_br,0.65*msize^2,off_wind_theta_m,'filled')
for i = 1:n_wave_age_bins
    plot(binned_wave_age(i)*[1 1],binned_theta_iqr(i,:),'-','Color','k','linewidth',2)
    plot(binned_wave_age_iqr(i,:),binned_theta_wave_minus_theta_wind(i)*[1 1],'-','Color','k','linewidth',2)
end
plot(binned_wave_age,binned_theta_wave_minus_theta_wind,'o','color','k','markerfacecolor','k','markeredgecolor','k','markersize',msize*1.25,'linewidth',2)
% Bin size floats clear of the marker so the whiskers stay visible underneath
text(12,30,'N = ','Color','k','FontSize',fsize,'FontWeight','bold','HorizontalAlignment','center')
for i = 1:length(binned_binsize)
    text(binned_wave_age(i),30,num2str(binned_binsize(i)),'Color','k','FontSize',fsize,'FontWeight','bold','HorizontalAlignment','center')
end
hold off
box on
ylim([-90 90])

cbar = colorbar;
cbar.Layout.Tile = 'north';

ax_struc(1).ax = gca;
ax_struc(1).ax.YTick = -90:30:90;
% ax_struc(1).ax.YDir = 'reverse';

cbar.Ticks = -90:30:90;
set(get(cbar,'Label'),'String','$\mathrm{\theta_m-\theta_{wind}\ [^\circ]}$','Interpreter','LaTeX')

xlabel('$\mathrm{c_p/u_*}$','Interpreter','latex')
ylabel('$\mathrm{\theta_{br}-\theta_{wind}\ [^\circ]}$','Interpreter','latex')

nexttile()
hold on
plot(wave_age_lims,[0 0],'--','linewidth',2,'Color',0.5*[1 1 1])
plot(wave_age,Momentum_flux_br,'o','markerfacecolor','k','markeredgecolor','k','markersize',msize,'linewidth',0.5)
scatter(wave_age,Momentum_flux_br,0.65*msize^2,off_wind_theta_m,'filled')
for i = 1:n_wave_age_bins
    plot(binned_wave_age(i)*[1 1],binned_momentum_flux_iqr(i,:),'-','Color','k', ...
        'linewidth',2,'HandleVisibility','off')
    plot(binned_wave_age_iqr(i,:),binned_momentum_flux(i)*[1 1],'-','Color','k', ...
        'linewidth',2,'HandleVisibility','off')
end
h_tau_br = plot(binned_wave_age,binned_momentum_flux,'o','color','k','markerfacecolor','k','markeredgecolor','k','markersize',msize*1.25,'linewidth',2);
h_tau_tot = plot(binned_wave_age,rho_a*binned_ustar.^2,'s','color','k','markerfacecolor','none','markeredgecolor','k','markersize',msize*1.75,'linewidth',2.5);
h_tau_visc = plot(binned_wave_age,binned_viscous_stress,'d','color','k','markerfacecolor','none','markeredgecolor','k','markersize',msize*1.75,'linewidth',2.5);
hold off
box on
ylim([5e-5 1e0])

H = [h_tau_tot h_tau_br h_tau_visc];
L = {'$\tau_{\mathrm{tot}}$','$\tau_{\mathrm{br}}$','$\tau_{\mathrm{visc}}$'};

legend(H,L,'Interpreter','latex','Location','southwest')

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
    xlim(wave_age_lims)
    ax_struc(n).ax.XScale = wc.axis_scale;
    ax_struc(n).ax.XTick = wc.axis_ticks;
    text(text_x,text_y,labels{n},'Units','normalized','FontSize',fsize,'HorizontalAlignment','left')
end

tile_cleaner(ax_struc,tlayout)
tlayout.TileSpacing = 'tight';