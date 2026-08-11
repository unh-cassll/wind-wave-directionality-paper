%
function all_wind_speed_stress_profiles(fignum,fsize)

load('data/ASIT2019_compiled_flux_timeseries.mat')

% Trim off tower measurements (keeping foremast)
T_degC(:,4:5) = [];
U_m_s(:,4:5) = [];
tw_covariance_m_s_degC(:,4:5) = [];
uw_covariance_m2_s2(:,4:5) = [];
vw_covariance_m2_s2(:,4:5) = [];
z_m(:,4:5) = [];

g = 9.81;
kappa = 0.4;

ustar_m_s = (uw_covariance_m2_s2.^2+vw_covariance_m2_s2.^2).^(0.25);

U10_extrap = median(ustar_m_s(:,1:3),2,'omitnan')/kappa.*log(10./z_m(:,3)) + U_m_s(:,3);

% Oncoming-wind gate, on the flux timeseries' own wind direction
mask = oncoming_wind_mask(D_deg_true(:));

U10_extrap = mask.*U10_extrap;

tau_downwind_N_m2 = -rho_a_kg_m3.*uw_covariance_m2_s2;
tau_crosswind_N_m2 = -rho_a_kg_m3.*vw_covariance_m2_s2;

ustar_m_s = (uw_covariance_m2_s2.^2+vw_covariance_m2_s2.^2).^0.25;
tau_N_m2 = sqrt(tau_downwind_N_m2.^2+tau_crosswind_N_m2.^2);

dU = 2;
N_winds = 6;
wind_speed_centers = dU:dU:N_winds*dU;

tau_profile = NaN*ones(length(wind_speed_centers),3);
z_profile = tau_profile;
U_profile = tau_profile;

N = NaN*ones(length(wind_speed_centers),1);

fit_struc = struct();

z_extend = (0:0.1:15)';

for n = 1:length(wind_speed_centers)

    U_low = wind_speed_centers(n) - dU/2;
    U_high = wind_speed_centers(n) + dU/2;

    inds_consider = U10_extrap >= U_low & U10_extrap < U_high;

    z_selection = z_m(inds_consider,:);
    U_selection = U_m_s(inds_consider,:);
    tau_selection = tau_N_m2(inds_consider,:);

    for m = 1:size(z_m,2)

        taubit = tau_selection(:,m);

        sorted_tau = sort(taubit);
        L = sum(~isnan(sorted_tau));
        tau99 = sorted_tau(floor(0.99*L));

        tau_selection(tau_selection > tau99) = NaN;

    end

    big_z = reshape(z_selection,[],1);
    big_U = reshape(U_selection,[],1);
    big_tau = reshape(tau_selection,[],1);

    U_profile(n,:) = mean(U_selection,1,'omitnan');
    tau_profile(n,:) = mean(tau_selection,1,'omitnan');
    z_profile(n,:) = mean(z_selection,1,'omitnan');

    inds_keep = ~isnan(big_tau);

    N(n) = sum(inds_keep);

    mdl = fitlm(big_z(inds_keep),big_tau(inds_keep));
    [tau_predict,tau_ci]=predict(mdl,z_extend(:));

    fit_struc(n).tau_predict = tau_predict;
    fit_struc(n).tau_ci = tau_ci;

end

cmap = magma(length(wind_speed_centers));

msize = 25;

figure(fignum);clf
set(gcf,'Position',[50 50 1100 550])
tlayout = tiledlayout(1,2);

nexttile(1)
hold on
plot(U_profile',z_profile','k.-','markersize',msize,'linewidth',4)
plot(U_profile',z_profile','.-','markersize',msize*0.92,'linewidth',2.75)
hold off
box on
xlim([0 15])
ylim([0 16])
colororder(cmap)
xlabel('$\mathrm{U\ [m\ s^{-1}]}$','Interpreter','LaTeX')
ylabel('$\mathrm{z\ [m]}$','Interpreter','LaTeX')
ax_struc(1).ax = gca;

profile_inds = z_extend >= min(z_profile(:,1)) & z_extend <= max(z_profile(:,3));

nexttile(2)
hold on
for n = 1:length(wind_speed_centers)
    tau_predict = fit_struc(n).tau_predict;
    tau_ci = fit_struc(n).tau_ci;
    plot(tau_predict,z_extend,'--','linewidth',1.5,'Color','k')
    X = [tau_ci(:,1); flipud(tau_ci(:,2))];
    Y = [z_extend(:); flipud(z_extend(:))];
    f = fill(X,Y,cmap(n,:));
    f.FaceAlpha = 0.2;
    points_profile = plot(tau_profile(n,:),z_profile(n,:),'o','markerfacecolor',cmap(n,:),'markeredgecolor','k','markersize',10,'linewidth',1.0);
    point_surface = plot(mean(tau_predict(profile_inds)),0,'p','markerfacecolor',cmap(n,:),'markeredgecolor','k','markersize',15,'linewidth',1.0);
end
hold off
box on
xlim([0 0.4])
ylim([0 16])
colororder(magma(length(wind_speed_centers)))
xlabel('$\mathrm{\tau\ [N\ m^{-2}]}$','Interpreter','LaTeX')
ylabel('$\mathrm{z\ [m]}$','Interpreter','LaTeX')

cbar = colorbar;
set(get(cbar,'Title'),'String','$\mathrm{U_{10}\ [m\ s^{-1}]}$','Interpreter','LaTeX')
colormap(magma(length(wind_speed_centers)))
clim([wind_speed_centers(1) wind_speed_centers(end)] + [-1 1]*dU/2)
cbar.Ticks = wind_speed_centers(1)-dU/2:dU:wind_speed_centers(end)+dU/2;

ax_struc(2).ax = gca;

panel_labels = {'(a)','(b)'};
for panel_ind = 1:2
    nexttile(panel_ind)
    text(0.05,1.05,panel_labels{panel_ind},'units','normalized','horizontalalignment','center','FontSize',fsize)
end

tile_cleaner(ax_struc,tlayout)

% nexttile(1)
% text(1.0,16,'N= ','Color','k','HorizontalAlignment','center','FontSize',16)
% for n = 1:length(N)
%     text(U_profile(n,end),16,num2str(N(n)),'Color','k','HorizontalAlignment','center','FontSize',16)
% end

