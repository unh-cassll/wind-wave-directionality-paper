% Subrange direction relative to the wind, against wave age, with runs split
% into the wind-fetch categories of breaking_direction_by_fetch: offshore
% < 10 km, cross-shore 10-100 km, onshore >= 100 km, via asit_fetch_vs_azimuth
%
% Breaking records categorize on their own wind record; the spectral and
% current panels share the run-grid wind
%
% Nathan Laxague 2026
%
function wind_wave_subrange_directions_by_fetch(fignum,fsize)

s = load('data/global_figure_settings.mat');
wave_age_lims = s.wave_age_lims;

% Subrange directions and their misalignment with the wind
m = subrange_wind_misalignment();

wave_age = m.wave_age;

% Cap the wave age at the axis limit so nothing plots outside the box
wave_age(wave_age>wave_age_lims(2)) = NaN;
m.breakers.wave_age(m.breakers.wave_age>wave_age_lims(2)) = NaN;
m.current.wave_age(m.current.wave_age>wave_age_lims(2)) = NaN;

% Wind-fetch categories, as in breaking_direction_by_fetch; crimson, teal and
% violet from the figure_style color order. cat_names are struct-field safe,
% cat_labels are what the legend shows
cat_names = {'offshore','crossshore','onshore'};
cat_labels = {'offshore','cross-shore','onshore'};
n_cat = length(cat_names);

color_order = get(groot,'DefaultAxesColorOrder');
cat_colors = color_order([4 2 1],:);

% Run grid takes the gated run-grid wind; breaking records take their own
fetch_run_km = asit_fetch_vs_azimuth(mod(m.wind_dir_deg_going_to+180,360));
cat_run = fetch_category(fetch_run_km);

Q = load('data/integrated_wave_breaking_quantities.mat','wdir_deg');
fetch_br_km = asit_fetch_vs_azimuth(mod(Q.wdir_deg+180,360));
cat_br = fetch_category(fetch_br_km);

D_particulars = {m.dtheta_eq,m.dtheta_sat,m.dtheta_sg,m.dtheta_gc, ...
    m.breakers.D_particular,m.current.D_particular};
D_refs = {m.dtheta_m,m.dtheta_m,m.dtheta_m,m.dtheta_m, ...
    m.breakers.D_ref,m.current.D_ref};
blocks = {m.block_id,m.block_id,m.block_id,m.block_id, ...
    m.breakers.block_id,m.block_id};
cats = {cat_run,cat_run,cat_run,cat_run,cat_br,cat_run};
wave_ages = {wave_age,wave_age,wave_age,wave_age, ...
    m.breakers.wave_age,m.current.wave_age};

names = {'equilibrium','saturation','short gravity','gravity-capillary','breaking','current'};
labels = {'(a)','(b)','(c)','(d)','(e)','(f)'};
label_x = 0.015;
label_y = 0.93;

msize = 12;
lw = 0.5;

figure(fignum);clf
tlayout = tiledlayout(3,2);
ax_struc = struct();

count_holder = zeros(6,n_cat);
stats_holder = struct();

% Right-margin distribution glyphs: drawn just outside the plot box, one
% glyph column per summarized category
glyph_x = [106 115 125];
glyph_cats = [1 2 3];   % offshore, crossshore, onshore

for n = 1:6

    nexttile(n)
    hold on

    x = wave_ages{n}(:);
    plot(wave_age_lims,[0 0],'--','Color',0.5*[1 1 1],'linewidth',2)

    y = D_particulars{n}(:);
    refv = D_refs{n}(:);
    blk = blocks{n}(:);
    cat_vec = cats{n}(:);

    for cat_num = 1:n_cat

        inds_cat = cat_vec==cat_num & ~isnan(x) & ~isnan(y);
        count_holder(n,cat_num) = sum(inds_cat);

        % Stash the plotted populations for the summary-statistics table
        stats_holder(n).(['wa_' cat_names{cat_num}]) = x(inds_cat);
        stats_holder(n).(['dth_' cat_names{cat_num}]) = y(inds_cat);
        stats_holder(n).(['ref_' cat_names{cat_num}]) = refv(inds_cat);
        stats_holder(n).(['blk_' cat_names{cat_num}]) = blk(inds_cat);

        scatter(x(inds_cat),y(inds_cat),0.5*msize^2,cat_colors(cat_num,:), ...
            'filled','MarkerFaceAlpha',0.15,'MarkerEdgeColor','k', ...
            'MarkerEdgeAlpha',0.15,'LineWidth',lw);

    end

    % Per-category equal-count wave-age bin medians: at most five bins, holding
    % the per-bin occupancy at five or more by dropping bins instead
    max_bins = 5;
    min_per_bin = 5;

    for cat_num = 1:n_cat

        inds_med = cat_vec==cat_num & ~isnan(x) & ~isnan(y);

        n_bins = min(max_bins,floor(sum(inds_med)/min_per_bin));
        if n_bins < 1
            continue
        end

        [x_sort,order] = sort(x(inds_med));
        y_sort = y(inds_med);
        y_sort = y_sort(order);

        bin_bounds = round(linspace(0,length(x_sort),n_bins+1));

        for bin_num = 1:n_bins

            inds_bin = bin_bounds(bin_num)+1:bin_bounds(bin_num+1);

            x_med = median(x_sort(inds_bin));
            y_med = median(y_sort(inds_bin));

            x_iqr = prctile(x_sort(inds_bin),[25 75]);
            y_iqr = prctile(y_sort(inds_bin),[25 75]);

            plot(x_iqr,y_med*[1 1],'-','Color',cat_colors(cat_num,:),'linewidth',2)
            plot(x_med*[1 1],y_iqr,'-','Color',cat_colors(cat_num,:),'linewidth',2)
            plot(x_med,y_med,'s','markerfacecolor',cat_colors(cat_num,:), ...
                'markeredgecolor','k','markersize',1.5*msize,'linewidth',lw)

        end

    end

    % Margin glyphs: square at the circular median, thick line over the IQR,
    % thin line over the IDR, skipping empty categories
    for g = 1:length(glyph_cats)

        dth_g = stats_holder(n).(['dth_' cat_names{glyph_cats(g)}]);
        if isempty(dth_g)
            continue
        end
        [med_g,dev_g] = circular_median_deg(dth_g);

        % [10 25 75 90] percentiles: the IDR then the IQR
        q = prctile(dev_g,[10 25 75 90]);

        col = cat_colors(glyph_cats(g),:);

        plot(glyph_x(g)*[1 1],med_g+[q(1) q(4)],'-','Color',col, ...
            'linewidth',1.5,'Clipping','off')
        plot(glyph_x(g)*[1 1],med_g+[q(2) q(3)],'-','Color',col, ...
            'linewidth',5,'Clipping','off')
        plot(glyph_x(g),med_g,'s', ...
            'markerfacecolor',col,'markeredgecolor',col, ...
            'markersize',1.4*msize,'linewidth',lw,'Clipping','off')

    end

    hold off
    box on
    pbaspect([1 1*2/3 1])

    xlim(wave_age_lims)
    ylim([-90 90])
    ax_struc(n).ax = gca;
    ax_struc(n).ax.YTick = -90:45:90;
    ax_struc(n).ax.XScale = 'log';
    ax_struc(n).ax.XTick = [10 20 40 80];

    xlabel('$\mathrm{c_p/u_*}$','Interpreter','latex')
    ylabel('$\mathrm{\theta_{particular}-\theta_{wind}\ [^\circ]}$','Interpreter','latex')

    % Category legend in the current panel; NaN placeholders give
    % full-opacity swatches regardless of which medians are drawn
    if n == 6

        hold on
        h_key = gobjects(1,n_cat);
        for cat_num = 1:n_cat
            h_key(cat_num) = plot(NaN,NaN,'s','markerfacecolor',cat_colors(cat_num,:), ...
                'markeredgecolor','k','markersize',1.5*msize,'linewidth',lw);
        end
        hold off

        legend(h_key,cat_labels,'Location','southwest','FontSize',0.85*fsize)

    end

end

tile_cleaner(ax_struc,tlayout)

for n = 1:length(ax_struc)
    nexttile(n)
    text(1-label_x,label_y,labels{n},'Units','normalized','FontSize',fsize,'HorizontalAlignment','right')
    text(label_x,label_y,names{n},'color',[38 38 38]/255,'FontSize',fsize*1.25,'Units','Normalized','HorizontalAlignment','left')
end

tlayout.TileSpacing = 'compact';

count_table = array2table(count_holder,'VariableNames',cat_names,'RowNames',names);
disp(count_table)

% Per-regime, per-category summary statistics: circular median and dispersion
% of the misalignment, plus an LSQ regression of dtheta_particular on dtheta_m
% within each category, with a 95% CI from the wild bootstrap over acquisition
% blocks. Categories holding fewer than min_fit_n runs are tabulated with the
% case count alone: too few cases to support a slope, a median, or a percentile
% span. The wave-age t-tests print to screen only, documenting that the
% categories sample the same forcing regime
row_names = {'Equilibrium','Saturation','Short Gravity','Gravity-Capillary','Breaking','Current'};
min_fit_n = 5;

rng(0)   % reproducible wild bootstrap over acquisition blocks

% Forcing-regime equivalence, console only
p_wa = NaN*ones(6,2);
for n = 1:6
    [~,p_wa(n,1)] = ttest2(stats_holder(n).wa_onshore,stats_holder(n).wa_crossshore);
    [~,p_wa(n,2)] = ttest2(stats_holder(n).wa_onshore,stats_holder(n).wa_offshore);
end

fprintf('\nWave-age two-sample t-tests (forcing-regime equivalence):\n')
disp(array2table(p_wa,'VariableNames',{'p_WA_OnCross','p_WA_OnOff'},'RowNames',row_names))

if ~exist('figs/tex','dir'); mkdir('figs/tex'); end
fid = fopen('figs/tex/table_fetch_stats.tex','w');
fprintf(fid,'%% Auto-generated by wind_wave_subrange_directions_by_fetch.m; do not edit by hand\n');
fprintf(fid,'\\begin{tabular}{|l|l|l|l|l|l|l|l|}\n');
fprintf(fid,'\\hline\n');
fprintf(fid,['\\textbf{Subrange} & \\textbf{sector} & \\textbf{$N$} & ' ...
    '\\textbf{$\\tilde{\\Delta\\theta}$ [$^\\circ$]} & ' ...
    '\\textbf{IQR [$^\\circ$]} & \\textbf{IDR [$^\\circ$]} & ' ...
    '\\textbf{slope} & \\textbf{slope 95\\%% CI} \\\\ \\hline\n']);
stats_rows = cell(0,8);

for n = 1:6

    for cat_num = 1:n_cat

        dth_c = stats_holder(n).(['dth_' cat_names{cat_num}]);
        ref_c = stats_holder(n).(['ref_' cat_names{cat_num}]);
        blk_c = stats_holder(n).(['blk_' cat_names{cat_num}]);

        if isempty(dth_c)

            med_c = NaN;
            iqr_c = NaN;
            idr_c = NaN;

        else

            [med_c,dev] = circular_median_deg(dth_c);
            q = prctile(dev,[10 25 75 90]);
            iqr_c = q(3)-q(2);
            idr_c = q(4)-q(1);

        end

        % Within-category regression of the particular misalignment on the
        % mean-sea misalignment. The same case-count gate suppresses the
        % tabulated median and dispersion statistics
        inds_reg = isfinite(ref_c(:)) & isfinite(dth_c(:));

        if sum(inds_reg) >= min_fit_n

            cr = cluster_robust_slope(ref_c(inds_reg),dth_c(inds_reg),blk_c(inds_reg));
            slope_str = sprintf('%.3f',cr.slope);
            ci_str = sprintf('(%.3f, %.3f)',cr.CI_wild(1),cr.CI_wild(2));
            slope_val = cr.slope;
            ci_val = cr.CI_wild(:)';

            med_str = sprintf('%.1f',med_c);
            iqr_str = sprintf('%.1f',iqr_c);
            idr_str = sprintf('%.1f',idr_c);

        else

            % Too few cases to support a slope, a median or a percentile span,
            % so the row carries the count alone
            slope_str = '--';
            ci_str = '--';
            slope_val = NaN;
            ci_val = [NaN NaN];

            med_str = '--';
            iqr_str = '--';
            idr_str = '--';

        end

        % Subrange label on the middle row of each sector triplet
        if cat_num == 2
            row_head = row_names{n};
        else
            row_head = '';
        end

        % Horizontal rule after the last sector of each subrange
        line_end = ' \\';
        if cat_num == n_cat
            line_end = ' \\ \hline';
        end

        fprintf(fid,'%s & %s & %d & %s & %s & %s & %s & %s%s\n', ...
            row_head,cat_labels{cat_num},sum(inds_reg),med_str,iqr_str,idr_str,slope_str,ci_str,line_end);

        stats_rows(end+1,:) = {sprintf('%s, %s',row_names{n},cat_labels{cat_num}), ...
            med_c,iqr_c,idr_c,slope_val,ci_val(1),ci_val(2),sum(inds_reg)}; %#ok<AGROW>

    end

end
fprintf(fid,'\\end{tabular}\n');
fclose(fid);

disp(cell2table(stats_rows,'VariableNames', ...
    {'Regime_Sector','Circ_Median','IQR','IDR','Slope','CI_Lower','CI_Upper','N_fit'}))


% Discrete wind-fetch category from fetch [km]: 1 offshore, 2 cross-shore,
% 3 onshore, 0 where the fetch is undefined (no gated wind record)
function cat_vec = fetch_category(fetch_km)

fetch_offshore_km = 10;
fetch_onshore_km = 100;

fetch_km = fetch_km(:);

cat_vec = 0*fetch_km;
cat_vec(fetch_km < fetch_offshore_km) = 1;
cat_vec(fetch_km >= fetch_offshore_km & fetch_km < fetch_onshore_km) = 2;
cat_vec(fetch_km >= fetch_onshore_km) = 3;
cat_vec(~isfinite(fetch_km)) = 0;


% Circular median: the sample angle minimizing the mean absolute circular
% distance to all others. Also returns each angle's signed circular deviation
% from that median, for circular-aware percentile spans
function [med_deg,dev_deg] = circular_median_deg(deg)

th = pi/180*deg(:);

pairwise = abs(angle(exp(1i*(th - th.'))));

[~,ind_min] = min(mean(pairwise,1));

med_deg = deg(ind_min);
dev_deg = 180/pi*angle(exp(1i*(th - th(ind_min))));
