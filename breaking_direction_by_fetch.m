%
% Breaking direction relative to the wind, and breaking stress, against wave
% age, with runs split into discrete wind-fetch categories: the ray-cast
% distance from ASIT to the regional coastline along the wind coming-from
% azimuth, via asit_fetch_vs_azimuth. Medians on the shared log-spaced
% wave-age bins
%
% Nathan Laxague 2026
%
function breaking_direction_by_fetch(fignum,fsize)

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

% Offshore-wind gate (fig 12's crest-gate counterpart), currently disabled
% Wdir_br(Wdir_br<270&Wdir_br>90) = NaN;

% Oncoming-wind gate; Wdir_br is the direction the wind blows towards.
% Breaking records carry bulk COARE wind, so they take the COARE half-angle
Wdir_br = oncoming_wind_mask(mod(Wdir_br+180,360),s.ASIT_COARE_wind_gate_deg).*Wdir_br;

U10_br(U10_br<0) = NaN;

Momentum_flux_br(Momentum_flux_br>0.5) = NaN;

off_wind_theta_br = compute_relative_angle(theta_br,Wdir_br);
off_wind_theta_m = compute_relative_angle(theta_m_br,Wdir_br);

% Fetch along the wind coming-from azimuth
fetch_km = asit_fetch_vs_azimuth(mod(Wdir_br+180,360));

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
% seas are also the runs whose Lambda(c) integral goes negative through the
% breaking strength b, and so cannot be drawn on the log axis
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

msize = 15;

% Discrete fetch categories, as in wind_fetch_category, but on the breaking
% records' own wind direction. Colors from the figure_style sweep (crimson,
% teal, violet); black marks the tau_tot/tau_visc references
fetch_offshore_km = 10;
fetch_onshore_km = 100;

color_order = get(groot,'DefaultAxesColorOrder');
cat_colors = color_order([4 2 1],:);
cat_names = {'offshore','cross-shore','onshore'};

fetch_cat = 0*fetch_km;
fetch_cat(fetch_km < fetch_offshore_km) = 1;
fetch_cat(fetch_km >= fetch_offshore_km & fetch_km < fetch_onshore_km) = 2;
fetch_cat(fetch_km >= fetch_onshore_km) = 3;
fetch_cat(~isfinite(fetch_km)) = 0;

ref_color = [0 0 0];

figure(fignum);clf
tlayout = tiledlayout(2,1);
ax_struc = struct();

nexttile()
hold on
plot(wave_age_lims,[0 0],'--','linewidth',2,'Color',0.5*[1 1 1])
h_cat = category_bin_medians(wave_age,off_wind_theta_br,fetch_cat,cat_colors,msize,wave_age_edges);
hold off
box on
ylim([-90 90])
ax_struc(1).ax = gca;
ax_struc(1).ax.YTick = -90:30:90;

xlabel('$\mathrm{c_p/u_*}$','Interpreter','latex')
ylabel('$\mathrm{\theta_{br}-\theta_{wind}\ [^\circ]}$','Interpreter','latex')

legend(h_cat,cat_names,'Location','northeast','FontSize',0.85*fsize)

nexttile()
hold on
% Reference stresses drawn first so the category medians sit on top
h_tau_tot = plot(binned_wave_age,rho_a*binned_ustar.^2,'s:','color',ref_color, ...
    'markerfacecolor','none','markeredgecolor',ref_color,'markersize',msize*1.25,'linewidth',2.5);
h_tau_visc = plot(binned_wave_age,binned_viscous_stress,'d:','color',ref_color, ...
    'markerfacecolor','none','markeredgecolor',ref_color,'markersize',msize*1.25,'linewidth',2.5);
category_bin_medians(wave_age,Momentum_flux_br,fetch_cat,cat_colors,msize,wave_age_edges)
hold off
box on
ylim([5e-4 1e0])

legend([h_tau_tot h_tau_visc],{'$\tau_{\mathrm{tot}}$','$\tau_{\mathrm{visc}}$'}, ...
    'Interpreter','latex','Location','southwest')

ax_struc(2).ax = gca;
ax_struc(2).ax.YScale = 'log';

xlabel('$\mathrm{c_p/u_*}$','Interpreter','latex')
ylabel('$\mathrm{\tau\ [N\ m^{-2}]}$','Interpreter','latex')

for n = 1:2
    nexttile(n)
    xlim(wave_age_lims)
    ax_struc(n).ax.XScale = wc.axis_scale;
    ax_struc(n).ax.XTick = wc.axis_ticks;
    text(text_x,text_y,labels{n},'Units','normalized','FontSize',fsize,'HorizontalAlignment','left')
end

tile_cleaner(ax_struc,tlayout)
tlayout.TileSpacing = 'tight';

n_per_cat = NaN*ones(1,length(cat_names));
for cat_num = 1:length(cat_names)
    n_per_cat(cat_num) = sum(fetch_cat(:)==cat_num & isfinite(wave_age(:)) & isfinite(off_wind_theta_br(:)));
end
disp(array2table(n_per_cat,'VariableNames',cat_names))


% Medians per category on the shared log-spaced wave-age bins, requiring at
% least three runs per bin. Square at the median, IQR whiskers in both variables
function h = category_bin_medians(x,y,cat_vec,colors,msize,edges)

min_per_bin = 3;

n_cat = size(colors,1);

h = gobjects(1,n_cat);

for cat_num = 1:n_cat

    for bin_num = 1:length(edges)-1

        inds_bin = cat_vec(:)==cat_num & isfinite(x(:)) & isfinite(y(:)) & ...
            x(:) >= edges(bin_num) & x(:) < edges(bin_num+1);

        if sum(inds_bin) < min_per_bin
            continue
        end

        x_med = median(x(inds_bin));
        y_med = median(y(inds_bin));

        x_iqr = prctile(x(inds_bin),[25 75]);
        y_iqr = prctile(y(inds_bin),[25 75]);

        plot(x_iqr,y_med*[1 1],'-','Color',colors(cat_num,:),'linewidth',2)
        plot(x_med*[1 1],y_iqr,'-','Color',colors(cat_num,:),'linewidth',2)

        h(cat_num) = plot(x_med,y_med,'s','markerfacecolor',colors(cat_num,:),'markeredgecolor','k', ...
            'markersize',1.5*msize,'linewidth',0.5);

    end

end
