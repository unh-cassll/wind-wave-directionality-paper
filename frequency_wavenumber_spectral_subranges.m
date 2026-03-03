%
function frequency_wavenumber_spectral_subranges(fignum,fsize)

cmap = flipud(spectral(7));
violet = cmap(1,:);
teal = cmap(2,:);
crimson = cmap(7,:);

in_nc_name = 'data/ASIT2019_supporting_environmental_observations.nc';

kappa = 0.4;

EC_U_m_s = ncread(in_nc_name,'EC_U_m_s');
EC_ustar_m_s = ncread(in_nc_name,'EC_ustar_m_s');
EC_z_m_above_water = ncread(in_nc_name,'EC_z_m_above_water');
EC_z0_m = EC_z_m_above_water.*exp(-kappa*EC_U_m_s./EC_ustar_m_s);
EC_U10_m_s = EC_ustar_m_s/kappa.*log(10./EC_z0_m);

load('data/frequency_spect_range_limits.mat')

f_p = f_E;

f_eq_start = f_eq_start(:);
f_eq_end = f_eq_end(:);
f_sat_end = f_sat_end(:);


text_x = 0.05;
text_y = 0.95;
labels = {'(a)','(b)','(c)','(d)'};

scat_size = 49;

fit_U10 = (1:0.1:14)';

lm_eq_start = fitlm(EC_U10_m_s,log10(f_eq_start));
[f_eq_start_predict,f_eq_start_ci] = predict(lm_eq_start,fit_U10);
f_eq_start_predict = 10.^f_eq_start_predict;
f_eq_start_ci = 10.^f_eq_start_ci;

lm_eq_end = fitlm(EC_U10_m_s,log10(f_eq_end));
[f_eq_end_predict,f_eq_end_ci] = predict(lm_eq_end,fit_U10);
f_eq_end_predict = 10.^f_eq_end_predict;
f_eq_end_ci = 10.^f_eq_end_ci;

lm_sat_end = fitlm(EC_U10_m_s,log10(f_sat_end));
[f_sat_end_predict,f_sat_end_ci] = predict(lm_sat_end,fit_U10);
f_sat_end_predict = 10.^f_sat_end_predict;
f_sat_end_ci = 10.^f_sat_end_ci;

s = load('data/global_figure_settings.mat');
dU = s.dU;

U10_left = s.U_low:dU:s.U_high-dU;
U10_right = U10_left + dU;
U10_lims = [U10_left(:) U10_right(:)];
U10_mean = NaN*U10_lims(:,1);
f_block = NaN*ones(size(U10_lims,1),3);
f_norm_block = f_block;
mask = 0*EC_U10_m_s+1;
mask(isnan(f_eq_start)|isnan(f_eq_end)|isnan(f_sat_end)) = NaN;
for n = 1:size(U10_lims,1)
    inds_within = EC_U10_m_s >= U10_lims(n,1) & EC_U10_m_s < U10_lims(n,2);
    U10_mean(n) = median(mask(inds_within).*EC_U10_m_s(inds_within),'all','omitnan');
    f_block(n,:) = [median(f_eq_start(inds_within),'all','omitnan') median(f_eq_end(inds_within),'all','omitnan') median(f_sat_end(inds_within),'all','omitnan')];
    f_norm_block(n,:) = [median(f_eq_start(inds_within)./f_p(inds_within),'all','omitnan') median(f_eq_end(inds_within)./f_p(inds_within),'all','omitnan') median(f_sat_end(inds_within)./f_p(inds_within),'all','omitnan')];
end

fA = 0.2;

U_lims = [0 15];

yticks = 10.^(-1:1:2);
yticklabels = {'0.1','1','10','100'};

figure(fignum);clf
tlayout = tiledlayout(2,2);
ax_struc = struct();

nexttile(1)
hold on
h_eq_start = scatter(EC_U10_m_s,f_eq_start,scat_size,violet,'filled');
h_eq_end = scatter(EC_U10_m_s,f_eq_end,scat_size,teal,'filled');
h_sat_end = scatter(EC_U10_m_s,f_sat_end,scat_size,crimson,'filled');
h_eq_start_binned = scatter(U10_mean,f_block(:,1),2*scat_size,violet,'filled');
h_eq_end_binned = scatter(U10_mean,f_block(:,2),2*scat_size,teal,'filled');
h_sat_end_binned = scatter(U10_mean,f_block(:,3),2*scat_size,crimson,'filled');
plot(fit_U10,f_eq_start_predict,'-','Color',violet,'linewidth',1.5)
plot(fit_U10,f_eq_start_ci(:,1),'--','Color',violet,'linewidth',1.5)
plot(fit_U10,f_eq_start_ci(:,2),'--','Color',violet,'linewidth',1.5)
plot(fit_U10,f_eq_end_predict,'-','Color',teal,'linewidth',1.5)
plot(fit_U10,f_eq_end_ci(:,1),'--','Color',teal,'linewidth',1.5)
plot(fit_U10,f_eq_end_ci(:,2),'--','Color',teal,'linewidth',1.5)
plot(fit_U10,f_sat_end_predict,'-','Color',crimson,'linewidth',1.5)
plot(fit_U10,f_sat_end_ci(:,1),'--','Color',crimson,'linewidth',1.5)
plot(fit_U10,f_sat_end_ci(:,2),'--','Color',crimson,'linewidth',1.5)
hold off
box on
xlim(U_lims)
% ylim([0 3])
ylim([1e-1 5e0])
% xlabel('U_{10} [m s^{-1}]')
ylabel('f [Hz]')
H = [h_eq_start_binned h_eq_end_binned h_sat_end_binned];
L = {'start of eq. range','eq. \rightarrow sat.','end of sat. range'};
legend(H,L,'location','northeast')

h_eq_start.MarkerFaceAlpha = fA;
h_eq_end.MarkerFaceAlpha = fA;
h_sat_end.MarkerFaceAlpha = fA/2;

h_eq_start_binned.Marker = 's';
h_eq_end_binned.Marker = 's';
h_sat_end_binned.Marker = 's';

ax_struc(1).ax = gca;
ax_struc(1).ax.YScale = 'log';
ax_struc(1).ax.YTick = yticks;
ax_struc(1).ax.YTickLabel = yticklabels;

lm_eq_start = fitlm(EC_U10_m_s,log10(f_eq_start./f_p(:)));
[f_eq_start_predict,f_eq_start_ci] = predict(lm_eq_start,fit_U10);
f_eq_start_predict = 10.^f_eq_start_predict;
f_eq_start_ci = 10.^f_eq_start_ci;

lm_eq_end = fitlm(EC_U10_m_s,log10(f_eq_end./f_p(:)));
[f_eq_end_predict,f_eq_end_ci] = predict(lm_eq_end,fit_U10);
f_eq_end_predict = 10.^f_eq_end_predict;
f_eq_end_ci = 10.^f_eq_end_ci;

lm_sat_end = fitlm(EC_U10_m_s,log10(f_sat_end./f_p(:)));
[f_sat_end_predict,f_sat_end_ci] = predict(lm_sat_end,fit_U10);
f_sat_end_predict = 10.^f_sat_end_predict;
f_sat_end_ci = 10.^f_sat_end_ci;

nexttile(3)
hold on
h_eq_start = scatter(EC_U10_m_s,f_eq_start./f_p(:),scat_size,violet,'filled');
h_eq_end = scatter(EC_U10_m_s,f_eq_end./f_p(:),scat_size,teal,'filled');
h_sat_end = scatter(EC_U10_m_s,f_sat_end./f_p(:),scat_size,crimson,'filled');
h_eq_start_binned = scatter(U10_mean,f_norm_block(:,1),2*scat_size,violet,'filled');
h_eq_end_binned = scatter(U10_mean,f_norm_block(:,2),2*scat_size,teal,'filled');
h_sat_end_binned = scatter(U10_mean,f_norm_block(:,3),2*scat_size,crimson,'filled');
plot(fit_U10,f_eq_start_predict,'-','Color',violet,'linewidth',1.5)
plot(fit_U10,f_eq_start_ci(:,1),'--','Color',violet,'linewidth',1.5)
plot(fit_U10,f_eq_start_ci(:,2),'--','Color',violet,'linewidth',1.5)
plot(fit_U10,f_eq_end_predict,'-','Color',teal,'linewidth',1.5)
plot(fit_U10,f_eq_end_ci(:,1),'--','Color',teal,'linewidth',1.5)
plot(fit_U10,f_eq_end_ci(:,2),'--','Color',teal,'linewidth',1.5)
plot(fit_U10,f_sat_end_predict,'-','Color',crimson,'linewidth',1.5)
plot(fit_U10,f_sat_end_ci(:,1),'--','Color',crimson,'linewidth',1.5)
plot(fit_U10,f_sat_end_ci(:,2),'--','Color',crimson,'linewidth',1.5)
plot(U_lims,2*[1 1],'-.','Color',violet,'linewidth',2)
hold off
box on
xlim(U_lims)
% ylim([0.5 50])
ylim([1 20])
xlabel('U_{10} [m s^{-1}]')
ylabel('f/f_E')

text(0.5,1.75,'Banner [1990]','Color',violet,'FontSize',fsize)

h_eq_start.MarkerFaceAlpha = fA;
h_eq_end.MarkerFaceAlpha = fA;
h_sat_end.MarkerFaceAlpha = fA/2;

h_eq_start_binned.Marker = 's';
h_eq_end_binned.Marker = 's';
h_sat_end_binned.Marker = 's';

ax_struc(3).ax = gca;
ax_struc(3).ax.YScale = 'log';
ax_struc(3).ax.YTick = yticks;
ax_struc(3).ax.YTickLabel = yticklabels;

%

cmap = flipud(spectral(7));
violet = cmap(1,:);
teal = cmap(2,:);
goldenrod = cmap(5,:);
crimson = cmap(7,:);

in_nc_name = 'data/ASIT2019_supporting_environmental_observations.nc';

g = 9.81;

kappa = 0.4;

water_depth_m = 15;

EC_U_m_s = ncread(in_nc_name,'EC_U_m_s');
EC_ustar_m_s = ncread(in_nc_name,'EC_ustar_m_s');
EC_z_m_above_water = ncread(in_nc_name,'EC_z_m_above_water');
EC_z0_m = EC_z_m_above_water.*exp(-kappa*EC_U_m_s./EC_ustar_m_s);
EC_U10_m_s = EC_ustar_m_s/kappa.*log(10./EC_z0_m);

U_sfc_mag_m_s = ncread(in_nc_name,'U_sfc_mag_m_s');


load('data/frequency_spect_range_limits.mat')
load('data/wavenumber_spect_range_limits.mat')

f_p = f_E;

f_p = f_p(:);
f_eq_start = f_eq_start(:);
f_eq_end = f_eq_end(:);
f_sat_end = f_sat_end(:);
k_sat_end = k_sat_end(:);

k_p_disp = NaN*f_p;
k_eq_start_disp = k_p_disp;
k_eq_end_disp = k_p_disp;
k_sat_end_disp = k_p_disp;

for n = 1:length(U_sfc_mag_m_s)

    try

        [c,~] = lindisp_with_current(2*pi*f_p(n),water_depth_m,U_sfc_mag_m_s(n));
        k_p_disp(n) = 2*pi*f_p(n)./c;

        [c,~] = lindisp_with_current(2*pi*f_eq_start(n),water_depth_m,U_sfc_mag_m_s(n));
        k_eq_start_disp(n) = 2*pi*f_eq_start(n)./c;

        [c,~] = lindisp_with_current(2*pi*f_eq_end(n),water_depth_m,U_sfc_mag_m_s(n));
        k_eq_end_disp(n) = 2*pi*f_eq_end(n)./c;

        [c,~] = lindisp_with_current(2*pi*f_sat_end(n),water_depth_m,U_sfc_mag_m_s(n));
        k_sat_end_disp(n) = 2*pi*f_sat_end(n)./c;

    end

end


scat_size = 49;

fit_U10 = (1:0.1:14)';

s = load('data/global_figure_settings.mat');
dU = s.dU;

U10_left = s.U_low:dU:s.U_high-dU;
U10_right = U10_left + dU;
U10_lims = [U10_left(:) U10_right(:)];
U10_mean = NaN*U10_lims(:,1);
k_block = U10_mean;
k_block_disp = NaN*ones(size(U10_lims,1),3);
k_norm_block = k_block;
k_norm_block_disp = k_block_disp;
mask = 0*EC_U10_m_s+1;
mask(isnan(f_eq_start)|isnan(f_eq_end)|isnan(f_sat_end)) = NaN;
for n = 1:size(U10_lims,1)
    inds_within = EC_U10_m_s >= U10_lims(n,1) & EC_U10_m_s < U10_lims(n,2);
    U10_mean(n) = median(mask(inds_within).*EC_U10_m_s(inds_within),'all','omitnan');
    k_block(n) = median(k_sat_end(inds_within),'all','omitnan');
    k_block_disp(n,:) = [median(k_eq_start_disp(inds_within),'all','omitnan') median(k_eq_end_disp(inds_within),'all','omitnan') median(k_sat_end_disp(inds_within),'all','omitnan')];
    k_norm_block_disp(n,:) = [median(k_eq_start_disp(inds_within)./k_p_disp(inds_within),'all','omitnan') median(k_eq_end_disp(inds_within)./k_p_disp(inds_within),'all','omitnan') median(k_sat_end_disp(inds_within)./k_p_disp(inds_within),'all','omitnan')];
    k_norm_block(n,:) = median(k_sat_end(inds_within)./k_p_disp(inds_within),'all','omitnan');
end


lm_eq_start = fitlm(EC_U10_m_s,log10(k_eq_start_disp));
[k_eq_start_predict,~] = predict(lm_eq_start,fit_U10);
k_eq_start_predict = 10.^k_eq_start_predict;

lm_eq_end = fitlm(EC_U10_m_s,log10(k_eq_end_disp));
[k_eq_end_predict,~] = predict(lm_eq_end,fit_U10);
k_eq_end_predict = 10.^k_eq_end_predict;

lm_sat_end = fitlm(EC_U10_m_s,log10(k_sat_end_disp));
[k_sat_end_predict,~] = predict(lm_sat_end,fit_U10);
k_sat_end_predict = 10.^k_sat_end_predict;

lm_sat_end_true = fitlm(EC_U10_m_s,log10(k_sat_end));
[k_sat_end_true_predict,k_sat_end_true_ci] = predict(lm_sat_end_true,fit_U10);
k_sat_end_true_predict = 10.^k_sat_end_true_predict;
k_sat_end_true_ci = 10.^k_sat_end_true_ci;

fA = 0.2;

U_lims = [0 15];

yticks = 10.^(-1:1:3);
yticklabels = {'0.1','1','10','100','1000'};

nexttile(2)
hold on
h_sat_end_true = scatter(EC_U10_m_s,k_sat_end,scat_size,crimson,'filled');
h_eq_start_binned = scatter(U10_mean,k_block_disp(:,1),2*scat_size,violet);
h_eq_end_binned = scatter(U10_mean,k_block_disp(:,2),2*scat_size,teal);
h_sat_end_binned = scatter(U10_mean,k_block_disp(:,3),2*scat_size,crimson);
h_sat_end_true_binned = scatter(U10_mean,k_block,2*scat_size,crimson,'filled');
plot(fit_U10,k_eq_start_predict,':','Color',violet,'linewidth',2)
plot(fit_U10,k_eq_end_predict,':','Color',teal,'linewidth',2)
plot(fit_U10,k_sat_end_predict,':','Color',crimson,'linewidth',2)
plot(fit_U10,k_sat_end_true_predict,'-','Color',crimson,'linewidth',1.5)
plot(fit_U10,k_sat_end_true_ci(:,1),'--','Color',crimson,'linewidth',1.5)
plot(fit_U10,k_sat_end_true_ci(:,2),'--','Color',crimson,'linewidth',1.5)
hold off
box on
xlim(U_lims)
ylim([0.1 100])
ylabel('k [rad m^{-1}]')
text(10.75,20,'true','Color',crimson,'FontWeight','bold','FontSize',fsize,'HorizontalAlignment','center')

h_sat_end_true.MarkerFaceAlpha = fA/2;

h_eq_start_binned.Marker = 's';
h_eq_end_binned.Marker = 's';
h_sat_end_binned.Marker = 's';
h_sat_end_true_binned.Marker = 's';
h_eq_start_binned.LineWidth = 1.5;
h_eq_end_binned.LineWidth = 1.5;
h_sat_end_binned.LineWidth = 1.5;

ax_struc(2).ax = gca;
ax_struc(2).ax.YScale = 'log';
ax_struc(2).ax.YTick = yticks;
ax_struc(2).ax.YTickLabel = yticklabels;

lm_eq_start = fitlm(EC_U10_m_s,log10(k_eq_start_disp./k_p_disp(:)));
[k_eq_start_predict,~] = predict(lm_eq_start,fit_U10);
k_eq_start_predict = 10.^k_eq_start_predict;

lm_eq_end = fitlm(EC_U10_m_s,log10(k_eq_end_disp./k_p_disp(:)));
[k_eq_end_predict,~] = predict(lm_eq_end,fit_U10);
k_eq_end_predict = 10.^k_eq_end_predict;

lm_sat_end = fitlm(EC_U10_m_s,log10(k_sat_end_disp./k_p_disp(:)));
[k_sat_end_predict,~] = predict(lm_sat_end,fit_U10);
k_sat_end_predict = 10.^k_sat_end_predict;

lm_sat_end_true = fitlm(EC_U10_m_s,log10(k_sat_end./k_p_disp(:)));
[k_sat_end_true_predict,k_sat_end_true_ci] = predict(lm_sat_end_true,fit_U10);
k_sat_end_true_predict = 10.^k_sat_end_true_predict;
k_sat_end_true_ci = 10.^k_sat_end_true_ci;

nexttile(4)
hold on
h_sat_end_true = scatter(EC_U10_m_s,k_sat_end./k_p_disp(:),scat_size,crimson,'filled');
h_eq_start_binned = scatter(U10_mean,k_norm_block_disp(:,1),2*scat_size,violet);
h_eq_end_binned = scatter(U10_mean,k_norm_block_disp(:,2),2*scat_size,teal);
h_sat_end_binned = scatter(U10_mean,k_norm_block_disp(:,3),2*scat_size,crimson);
h_sat_end_true_binned = scatter(U10_mean,k_norm_block,2*scat_size,crimson,'filled');
plot(fit_U10,k_eq_start_predict,':','Color',violet,'linewidth',2)
plot(fit_U10,k_eq_end_predict,':','Color',teal,'linewidth',2)
plot(fit_U10,k_sat_end_predict,':','Color',crimson,'linewidth',2)
plot(fit_U10,k_sat_end_true_predict,'-','Color',crimson,'linewidth',1.5)
plot(fit_U10,k_sat_end_true_ci(:,1),'--','Color',crimson,'linewidth',1.5)
plot(fit_U10,k_sat_end_true_ci(:,2),'--','Color',crimson,'linewidth',1.5)
hold off
box on
xlim(U_lims)
ylim([1 1000])
xlabel('U_{10} [m s^{-1}]')
ylabel('k/k_E')

h_sat_end_true.MarkerFaceAlpha = fA/2;

h_eq_start_binned.Marker = 's';
h_eq_end_binned.Marker = 's';
h_sat_end_binned.Marker = 's';
h_sat_end_true_binned.Marker = 's';
h_eq_start_binned.LineWidth = 1.5;
h_eq_end_binned.LineWidth = 1.5;
h_sat_end_binned.LineWidth = 1.5;

ax_struc(4).ax = gca;
ax_struc(4).ax.YScale = 'log';
ax_struc(4).ax.YTick = yticks;
ax_struc(4).ax.YTickLabel = yticklabels;

for n = 2:2:4
    ax_struc(n).ax.YAxisLocation = 'right';
end

for n = 1:2
    ax_struc(n).ax.XTickLabel = '';
end

for n = 1:4
    nexttile(n)
    text(text_x,text_y,labels{n},'Units','normalized','FontSize',fsize,'HorizontalAlignment','center')
end

% tile_cleaner(ax_struc,tlayout)

tlayout.TileSpacing = 'tight';
