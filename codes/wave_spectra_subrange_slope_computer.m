
close all;clear;clc

load ../../../../../colormaps/coolwarm.mat

short_wave_nc_name = '../data/ASIT2019_wave_spectra_stats_timeseries_empirical_gain.nc';
long_wave_nc_name = '../data/ASIT2019_EPSS_directional_spectra.nc';
supporting_nc_name = '../data/ASIT2019_supporting_environmental_observations.nc';

g = 9.81;

f_Hz_Pyxis = double(ncread(short_wave_nc_name,'f_Hz'));
k_rad_m_Pyxis = double(ncread(short_wave_nc_name,'k_rad_m'));

theta_rad_Pyxis = double(ncread(short_wave_nc_name,'theta_rad'));

SKTHETA = double(ncread(short_wave_nc_name,'S_k_theta'));
SFTHETA = double(ncread(short_wave_nc_name,'S_f_theta'));

dtheta = median(diff(theta_rad_Pyxis));

Sk = squeeze(sum(k_rad_m_Pyxis(:).*SKTHETA,2,'omitnan'))*dtheta;
Sf = squeeze(sum(SFTHETA,2,'omitnan'))*dtheta;

% k_rad_m_Pyxis(1) = [];
% Sk(1,:) = [];

avg_size = 2;

f_Hz_EPSS = double(ncread(long_wave_nc_name,'frequency'));
theta_rad_EPSS = double(ncread(long_wave_nc_name,'direction'))*pi/180;
F_f_m2_Hz_rad_EPSS = double(ncread(long_wave_nc_name,'F_f_d'))*180/pi;

dtheta = median(diff(theta_rad_EPSS));

Ff_EPSS = squeeze(sum(F_f_m2_Hz_rad_EPSS,2,'omitnan'))'*dtheta;

kappa = 0.4;

EC_U_m_s = ncread(supporting_nc_name,'EC_U_m_s');
EC_ustar_m_s = ncread(supporting_nc_name,'EC_ustar_m_s');
EC_z_m_above_water = ncread(supporting_nc_name,'EC_z_m_above_water');
EC_z0_m = EC_z_m_above_water.*exp(-kappa*EC_U_m_s./EC_ustar_m_s);
EC_U10_m_s = EC_ustar_m_s/kappa.*log(10./EC_z0_m);

U_sfc_mag_m_s = ncread(supporting_nc_name,'U_sfc_mag_m_s');
U_sfc_dir_deg = ncread(supporting_nc_name,'U_sfc_dir_deg');

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

T_s = f_Hz_combined.^-1;
T_s(1) = 0;
f_E = (trapz(f_Hz_combined,T_s.*Ff_combined)./trapz(f_Hz_combined,Ff_combined)).^-1;
f_p = NaN*f_E;
for n = 1:190
    try
        ind = find(Ff_EPSS(:,n)==max(Ff_EPSS(:,n),[],'all','omitnan'),1,'first');
        f_p(n) = f_Hz_EPSS(ind);
    end
end

water_depth_m = 15;
f_p = f_E;
[C_p,Cg_p] = lindisp_with_current(2*pi*f_p,water_depth_m,0);
k_p = 2*pi*f_p(:)./C_p(:);

% save('../data/ASIT2019_peak_wave_phase_speed.mat','f_p','k_p','C_p')

%%

f_new = logspace(log10(5e-2),log10(5e0),64)';
Ff_new = interp1(f_Hz_combined,Ff_combined,f_new);

Ff_new = smoothdata2(Ff_new,'movmedian',{11,1},'omitnan');

% log_f = log10(f_Hz_combined);
% log_Ff = log10(Ff_combined);

log_f = log10(f_new);
log_Ff = log10(Ff_new);

indices_per_window = 16;

f_holder_struc = struct();

colors = flipud(spectral(7));
violet = colors(1,:);
teal = colors(2,:);


for particular_ind = 1:size(Sf,2)

    e = [];

    try

        log_Ff_bit = log_Ff(:,particular_ind);  

        ind_p = find(log_Ff_bit==max(log_Ff_bit),1,'first');
        % log_Ff_bit(f_Hz_combined<f_Hz_combined(ind_p)) = NaN;

        f_slope = NaN*log_f;
        R2 = f_slope;
        SE = f_slope;

        % if particular_ind > 68

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

        % f_slope(f_new<1.25*f_new(ind_p)) = NaN;

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

        % buffer_inds = -1:1:1;
        % fit_eq_start = fit(y(buffer_inds+ind_eq_start)-2*z(buffer_inds+ind_eq_start),x(buffer_inds+ind_eq_start),'linear');
        % fit_eq_end = fit(y(buffer_inds+ind_eq_end),x(buffer_inds+ind_eq_end),'linear');
        % fit_sat_end = fit(y(buffer_inds+ind_sat_end)+2*z(buffer_inds+ind_sat_end),x(buffer_inds+ind_sat_end),'linear');

        f_eq_start = 10.^x(ind_eq_start);
        f_eq_end = 10.^x(ind_eq_end);
        f_sat_end = 10.^x(ind_sat_end);

        % f_eq_start = 10.^fit_eq_start(1);
        % f_eq_end = 10.^fit_eq_end(0);
        % f_sat_end = 10.^fit_sat_end(-0.5);

        % slope_fit = fit(y-2*z,x,'spline');
        % f_eq_start = 10.^slope_fit(1);
        % slope_fit = fit(y,x,'spline');
        % f_eq_end = 10.^slope_fit(0);
        % slope_fit = fit(y+2*z,x,'spline');
        % f_sat_end = 10.^slope_fit(-0.5);

        tail_slope = mean(f_new(f_new>f_sat_end),'omitnan');

        f_holder_struc(particular_ind).f_eq_start = f_eq_start;
        f_holder_struc(particular_ind).f_eq_end = f_eq_end;
        f_holder_struc(particular_ind).f_sat_end = f_sat_end;
        f_holder_struc(particular_ind).tail_slope = tail_slope;

        figure(1);clf
        tlayout = tiledlayout(2,1);
        ax_struc = struct();

        nexttile(1)
        hold on
        plot([1e-2 1e2],[1 1]*1,'r--','linewidth',2)
        plot([1e-2 1e2],[1 1]*0,'r:','linewidth',2)
        plot([1 1]*f_p(particular_ind),[-1 1]*3,'k-','linewidth',2)
        plot([1 1]*f_eq_start,[-1 1]*3,'k--','linewidth',2)
        plot([1 1]*f_eq_end,[-1 1]*3,'k-.','linewidth',2)
        plot([1 1]*f_sat_end,[-1 1]*3,'k:','linewidth',3)
        plot(f_new,f_slope,'.-','Color',violet,'markersize',15)
        for i = 1:length(f_slope)
            plot(f_new(i)*[1 1],f_slope(i)*[1 1]+[-1 1]*SE(i)*2,'-','linewidth',1.5,'Color',violet)
        end
        hold off
        box on
        ax_struc(1).ax = gca;
        ax_struc(1).ax.XScale = 'log';
        ax_struc(1).ax.YTick = -2:1:2;
        xlim([1e-2 1e1]*2)
        ylim([-1 1]*3)
        xlabel('f [Hz]')
        ylabel('"n" such that B(f)\proptof^{n}')

        nexttile(2)
        plot(f_Hz_EPSS,(2*pi*f_Hz_EPSS).^5/g^2.*Ff_EPSS(:,particular_ind),'-',f_Hz_Pyxis,((2*pi*(f_Hz_Pyxis)).^1).*Sf(:,particular_ind),'-','linewidth',2)
        hold on
        plot([1 1]*f_p(particular_ind),[1e-5 1e0],'k-','linewidth',2)
        plot([1 1]*f_eq_start,[1e-5 1e0],'k--','linewidth',2)
        plot([1 1]*f_eq_end,[1e-5 1e0],'k-.','linewidth',2)
        plot([1 1]*f_sat_end,[1e-5 1e0],'k:','linewidth',3)
        s_combined = scatter(f_new,(2*pi*f_new).^5/g^2.*Ff_new(:,particular_ind),60,f_slope,'filled');
        % s_Pyxis = scatter(f_Hz_Pyxis,((2*pi*(f_Hz_Pyxis)).^1).*Sf(:,particular_ind),60,f_slope_Pyxis);colorbar
        hold off
        box on
        ax_struc(2).ax = gca;
        ax_struc(2).ax.XScale = 'log';
        ax_struc(2).ax.YScale = 'log';
        xlim([1e-2 1e1]*2)
        ylim([1e-4 1e0])
        clim([-1 1]*3)
        colormap(coolwarm)
        xlabel('f [Hz]')
        ylabel('B(f) [rad]')
        cbar = colorbar;
        set(get(cbar,'Title'),'String','"n"')


        s_combined.Marker = 's';
        % s_Pyxis.Marker = 'o';
        s_combined.LineWidth = 1.5;

        tile_cleaner(ax_struc,tlayout)

    catch e

    end

    if ~isempty(e)

        f_holder_struc(particular_ind).f_eq_start = NaN;
    f_holder_struc(particular_ind).f_eq_end = NaN;
    f_holder_struc(particular_ind).f_sat_end = NaN;
    f_holder_struc(particular_ind).tail_slope = NaN;

    end

    % pause(0.1)
    pause()

end

f_eq_start = [f_holder_struc.f_eq_start];
f_eq_end = [f_holder_struc.f_eq_end];
f_sat_end = [f_holder_struc.f_sat_end];
tail_slope = [f_holder_struc.tail_slope]-5;

% save('../data/frequency_spect_range_limits.mat','f_eq_start','f_eq_end','f_sat_end','tail_slope','f_E')

%%

k_new = logspace(log10(k_rad_m_Pyxis(1)),log10(100),64)';
Sk_new = interp1(k_rad_m_Pyxis,Sk,k_new);

log_k_Pyxis = log10(k_new);

indices_per_window = 24;

k_holder_struc = struct();

hard_code_inds = [45 46 51];
hard_code_k = k_rad_m_Pyxis([5 4 3]);

for particular_ind = 1:size(Sf,2)

    try

        log_Sk_Pyxis = log10(Sk_new(:,particular_ind));

        k_slope_Pyxis = NaN*log_Sk_Pyxis;
        R2_Pyxis = k_slope_Pyxis;
        SE_Pyxis = k_slope_Pyxis;

        % i_high = find(k_new>70,1,'first');

        for i = indices_per_window/2:length(log_Sk_Pyxis)-indices_per_window/2+1

            istart = i - indices_per_window/2 + 1;
            iend = i + indices_per_window/2 - 1;

            x = log_k_Pyxis(istart:iend);
            y = log_Sk_Pyxis(istart:iend);

            lm = fitlm(x,y);

            k_slope_Pyxis(i) = lm.Coefficients.Estimate(2);
            R2_Pyxis(i) = lm.Rsquared.Ordinary;
            SE_Pyxis(i) = lm.Coefficients.SE(2);

        end

        i_high = length(log_Sk_Pyxis)-indices_per_window/2+1;

        k_slope_Pyxis = k_slope_Pyxis + 1;

        % ind_peak

        ind_sat_end = find(k_slope_Pyxis+2*SE_Pyxis<-0.25,1,'first');
        if isempty(ind_sat_end) | ind_sat_end == i_high
            k_sat_end = NaN;
        else
            k_sat_end = k_new(ind_sat_end);% + indices_per_window/2*median(diff(k_rad_m_Pyxis));
        end
        

        tail_slope = mean(k_slope_Pyxis(k_new>k_sat_end),'omitnan');

        k_holder_struc(particular_ind).k_sat_end = k_sat_end;
        k_holder_struc(particular_ind).tail_slope = tail_slope;

        colors = flipud(spectral(7));
        violet = colors(1,:);
        teal = colors(2,:);

        figure(11);clf
        tlayout = tiledlayout(2,1);
        ax_struc = struct();

        nexttile(1)
        hold on
        plot([1e0 1e3],[1 1]*0.5,'r--','linewidth',2)
        plot([1e0 1e3],[1 1]*0,'r:','linewidth',2)
        plot([1 1]*k_sat_end,[-1 1]*3,'k:','linewidth',3)
        plot(k_new,k_slope_Pyxis,'.-','Color',violet,'markersize',15)
        for i = 1:length(k_slope_Pyxis)
            plot(k_new(i)*[1 1],k_slope_Pyxis(i)*[1 1]+[-1 1]*SE_Pyxis(i)*2,'-','linewidth',1.5,'Color',violet)
        end
        hold off
        box on
        ax_struc(1).ax = gca;
        ax_struc(1).ax.XScale = 'log';
        xlim([1e0 1e3])
        ylim([-2 2])
        xlabel('k [rad m^{-1}]')
        ylabel('n, s.t. B(k)\proptok^{n}')

        nexttile(2)
        plot(k_rad_m_Pyxis,k_rad_m_Pyxis.*Sk(:,particular_ind),'k-','linewidth',2)
        hold on
        plot([1 1]*k_sat_end,[1e-5 1e0],'k:','linewidth',3)
        s_Pyxis = scatter(k_new,k_new.*Sk_new(:,particular_ind),60,k_slope_Pyxis);colorbar
        hold off
        box on
        ax_struc(2).ax = gca;
        ax_struc(2).ax.XScale = 'log';
        ax_struc(2).ax.YScale = 'log';
        xlim([1e0 1e3])
        ylim([1e-4 1e-1])
        clim([-1 1]*1)
        colormap(coolwarm)
        xlabel('k [rad m^{-1}]')
        ylabel('B(k)')

        s_Pyxis.Marker = 'o';
        s_Pyxis.LineWidth = 1.5;

        tile_cleaner(ax_struc,tlayout)

    end

    pause(0.1)
    % pause()

end

k_sat_end = [k_holder_struc.k_sat_end];
tail_slope = [k_holder_struc.tail_slope]-3;
k_sat_end(hard_code_inds) = hard_code_k;

% save('../data/wavenumber_spect_range_limits.mat','k_sat_end','tail_slope')
