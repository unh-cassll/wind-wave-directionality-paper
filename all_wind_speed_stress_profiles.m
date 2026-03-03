%
function all_wind_speed_stress_profiles(fignum,fsize)

load('data/ASIT2019_compiled_flux_timeseries.mat')

g = 9.81;
kappa = 0.4;

ustar_m_s = (uw_covariance_m2_s2.^2+vw_covariance_m2_s2.^2).^(0.25);

U10_extrap = median(ustar_m_s(:,1:3),2,'omitnan')/kappa.*log(10./z_m(:,3)) + U_m_s(:,3);

tau_downwind_N_m2 = -rho_a_kg_m3.*uw_covariance_m2_s2;
tau_crosswind_N_m2 = -rho_a_kg_m3.*vw_covariance_m2_s2;

ustar_m_s = (uw_covariance_m2_s2.^2+vw_covariance_m2_s2.^2).^0.25;
tau_N_m2 = sqrt(tau_downwind_N_m2.^2+tau_crosswind_N_m2.^2);

dU = 2;
N_winds = 9;
wind_speed_centers = dU:dU:N_winds*dU;

tau_profile = NaN*ones(length(wind_speed_centers),5);
ustar_profile = tau_profile;
z_profile = tau_profile;
U_profile = tau_profile;

N = NaN*ones(length(wind_speed_centers),1);

for n = 1:length(wind_speed_centers)

    U_low = wind_speed_centers(n) - dU/2;
    U_high = wind_speed_centers(n) + dU/2;

    inds_consider = U10_extrap >= U_low & U10_extrap < U_high;

    U_profile(n,:) = median(U_m_s(inds_consider,:),1,'omitnan');
    ustar_profile(n,:) = median(ustar_m_s(inds_consider,:),1,'omitnan');
    tau_profile(n,:) = median(tau_N_m2(inds_consider,:),1,'omitnan');
    z_profile(n,:) = median(z_m(inds_consider,:),1,'omitnan');

    N(n) = sum(inds_consider);

end

cmap = magma(length(wind_speed_centers));

figure(fignum);clf
set(gcf,'Position',[50 50 1100 550])
tlayout = tiledlayout(1,2);

nexttile(1)
hold on
plot(U_profile',z_profile','k.-','linewidth',3.5)
plot(U_profile',z_profile','.-','linewidth',2)
hold off
box on
xlim([0 20])
ylim([0 25])
colororder(cmap)
xlabel('U [m s^{-1}]')
ylabel('z [m]')
ax_struc(1).ax = gca;

nexttile(2)
hold on
plot(tau_profile',z_profile','k.-','linewidth',3.5)
plot(tau_profile',z_profile','.-','linewidth',2)
hold off
box on
xlim([5e-3 1e0])
ylim([0 25])
colororder(magma(length(wind_speed_centers)))
xlabel('\tau [N m^{-2}]')
ylabel('z [m]')

cbar = colorbar;
set(get(cbar,'Title'),'String','U_{10} [m s^{-1}]')
colormap(magma(length(wind_speed_centers)))
clim([wind_speed_centers(1) wind_speed_centers(end)] + [-1 1]*dU/2)
cbar.Ticks = wind_speed_centers(1)-dU/2:dU:wind_speed_centers(end)+dU/2;

ax_struc(2).ax = gca;

ax_struc(2).ax.XScale = 'log';

tile_cleaner(ax_struc,tlayout)

nexttile(1)
for n = 1:length(N)
    textborder(U_profile(n,end),23.5,num2str(N(n)), cmap(n,:), [0 0 0],'HorizontalAlignment','center','FontSize',12)
end
for n = 1:length(N)
    textborder(U_profile(n,end),23.5,num2str(N(n)), cmap(n,:), [0 0 0],'HorizontalAlignment','center','FontSize',12)
end

dim = [0.532 0.29 0.16 0.08];
str = 'mast';
a = annotation('textbox',dim,'String',str);
a.Rotation = 90;
a.HorizontalAlignment = 'center';
a.VerticalAlignment = 'middle';
a.FontSize = fsize;

dim = [0.532 0.74 0.06 0.08];
str = 'tower';
a = annotation('textbox',dim,'String',str);
a.Rotation = 90;
a.HorizontalAlignment = 'center';
a.VerticalAlignment = 'middle';
a.FontSize = fsize;