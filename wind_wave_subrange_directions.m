%
function wind_wave_subrange_directions(fignum,fsize)

clc

s = load('data/global_figure_settings.mat');
wave_age_lims = s.wave_age_lims;

% Subrange directions and their misalignment with the wind
m = subrange_wind_misalignment();

wave_age = m.wave_age;

t = table(m.dtheta_m,m.dtheta_eq,m.dtheta_sat,m.dtheta_sg,m.dtheta_gc);
t.Properties.VariableNames = {'\theta_{m}-\theta_{wind} [\circ]','\theta_{eq.}-\theta_{wind} [\circ]','\theta_{sat.}-\theta_{wind} [\circ]','\theta_{s.g.}-\theta_{wind} [\circ]','\theta_{g-c}-\theta_{wind} [\circ]'};

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
    holder_struc(n).blocks = m.block_id;

end

holder_struc(n+1).D_ref = m.breakers.D_ref;
holder_struc(n+1).D_particular = m.breakers.D_particular;
holder_struc(n+1).wave_age = m.breakers.wave_age;
holder_struc(n+1).blocks = m.breakers.block_id;

holder_struc(n+2).D_ref = m.current.D_ref;
holder_struc(n+2).D_particular = m.current.D_particular;
holder_struc(n+2).wave_age = m.current.wave_age;
holder_struc(n+2).blocks = m.block_id;

% order:
% mean
% equilibrium
% saturation
% short gravity
% gravity-capillary
% breakers
% current

msize = 8;
lw = 0.5;

fA = 0.2;

figure(fignum);clf
tlayout = tiledlayout(3,2);
ax_struc = struct();

stats_holder = struct();

rng(0)   % reproducible wild bootstrap for the tabulated slope CIs

for n = 1:6

    nexttile(n)
    hold on
    plot([-180 180],[0 0],'--','Color',0.5*[1 1 1],'linewidth',2)
    plot([0 0],[-180 180],'--','Color',0.5*[1 1 1],'linewidth',2)

    x = holder_struc(n+1).D_ref(:);
    y = holder_struc(n+1).D_particular(:);
    c = holder_struc(n+1).wave_age(:);
    blk = holder_struc(n+1).blocks(:);
    inds_keep = ~isnan(x) & ~isnan(y) & ~isnan(c);
    x = x(inds_keep);
    y = y(inds_keep);
    c = c(inds_keep);
    blk = blk(inds_keep);
    xplot = linspace(-90+dD/2,90-dD/2,100);

    inds_fit = x>-90+dD/2 & x<90-dD/2;
    D_fit_object = fitlm(x(inds_fit),y(inds_fit));
    [D_fit,D_fit_95CI] = predict(D_fit_object,xplot(:));
    D_lower = D_fit_95CI(:,1)';
    D_upper = D_fit_95CI(:,2)';

    % Slope CI from a wild bootstrap over acquisition blocks, so that the
    % interval tolerates the hour-to-hour correlation between runs
    cr = cluster_robust_slope(x(inds_fit),y(inds_fit),blk(inds_fit));

    stats_holder(n).line_slope = D_fit_object.Coefficients.Estimate(2);
    stats_holder(n).line_intercept = D_fit_object.Coefficients.Estimate(1);
    stats_holder(n).slope_CI_lower = cr.CI_wild(1);
    stats_holder(n).slope_CI_upper = cr.CI_wild(2);
    stats_holder(n).n_blocks = cr.G;
    stats_holder(n).R2 = D_fit_object.Rsquared.Adjusted;

    plot(x,y,'o','markerfacecolor','k','markeredgecolor','k','markersize',msize,'linewidth',lw)
    scatter(x,y,0.6*msize^2,c,'filled')
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
    pbaspect([1 0.75 1])

    xlim([-1 1]*dir_max)
    ylim([-1 1]*dir_max*3/4)
    ax_struc(n).ax = gca;
    ax_struc(n).ax.YTick = -360:45:360;
    ax_struc(n).ax.XTick = -360:45:360;
    ax_struc(n).ax.YDir = 'reverse';

    xlabel('$\mathrm{\theta_m-\theta_{wind}\ [^\circ]}$','Interpreter','latex')
    ylabel('$\mathrm{\theta_{particular}-\theta_{wind}\ [^\circ]}$','Interpreter','latex')

end

tile_cleaner(ax_struc,tlayout)

for n = 1:length(ax_struc)
    nexttile(n)
    text(label_x,label_y,labels{n},'Units','normalized','FontSize',fsize)
    text(1-label_x, label_y, names{n},'color', [38 38 38]/255,'FontSize',fsize*1.25,'Units','Normalized','HorizontalAlignment','right')
end

tlayout.TileSpacing = 'tight';

cbar.Layout.Tile = 'north';
set(get(cbar,'Label'),'String','$\mathrm{c_p/u_*}$','Interpreter','LaTeX')

stats_table = struct2table(stats_holder);
stats_table.Properties.VariableNames = {'Line_Slope', 'Line_Intercept', 'Slope_CI_Lower', 'Slope_CI_Upper', 'N_Blocks', 'R2'};
row_names = {'Equilibrium', 'Saturation', 'Short Gravity', 'Gravity-Capillary', 'Breakers', 'Current'};
stats_table.Properties.RowNames = row_names;
disp(stats_table);

% writetable(stats_table,'data/off_wind_wave_angle_table.csv','WriteRowNames',true)
