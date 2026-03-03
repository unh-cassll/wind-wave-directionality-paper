function freq_spect_integral_moments(example_run_ind,fignum,fsize)

clc

short_wave_spectra_nc_name = 'data/ASIT2019_wave_spectra_stats_timeseries_empirical_gain.nc';
long_wave_spectra_nc_name = 'data/ASIT2019_EPSS_directional_spectra.nc';

g = 9.806;
water_depth_m = 15;

f_Hz_EPSS = ncread(long_wave_spectra_nc_name,'frequency');
theta_rad_EPSS = double(ncread(long_wave_spectra_nc_name,'direction'))*pi/180;

f_Hz_Pyxis = ncread(short_wave_spectra_nc_name,'f_Hz');
theta_rad_Pyxis = ncread(short_wave_spectra_nc_name,'theta_rad');

SFTHETA = double(ncread(short_wave_spectra_nc_name,'S_f_theta'));
FFTHETA = double(ncread(long_wave_spectra_nc_name,'F_f_d'))*180/pi;

SFTHETA(SFTHETA==0) = NaN;
FFTHETA(FFTHETA==0) = NaN;

cmap = flipud(spectral(7));
violet = cmap(1,:);
teal = cmap(2,:);
crimson = cmap(7,:);

Ff_EPSS = squeeze(sum(FFTHETA(example_run_ind,:,:),2,'omitnan'))*median(diff(theta_rad_EPSS))';

dtheta = median(diff(theta_rad_Pyxis));

f_cut = 0.4;

Sf_Pyxis = smoothdata(squeeze(sum(SFTHETA(:,:,example_run_ind),2,'omitnan'))*dtheta,'movmean',5);

F_f_Pyxis = ((2*pi*f_Hz_Pyxis).^4/g^2).^-1.*Sf_Pyxis;
ind_trim_Pyxis = find(f_Hz_Pyxis>f_cut*1.1,1,'first');
ind_trim_EPSS = find(f_Hz_EPSS<f_cut*1,1,'last');

f_Hz_combined = [f_Hz_EPSS(2:ind_trim_EPSS); f_Hz_Pyxis(ind_trim_Pyxis:end)];
F_f_combined = [Ff_EPSS(2:ind_trim_EPSS); F_f_Pyxis(ind_trim_Pyxis:end)];

ind_p = find(F_f_combined==max(F_f_combined,[],'all','omitnan'),1,'first');
f_p = f_Hz_combined(ind_p);

f_m01 = trapz(f_Hz_combined,f_Hz_combined.*F_f_combined)/trapz(f_Hz_combined,F_f_combined);
f_m02 = sqrt(trapz(f_Hz_combined,f_Hz_combined.^2.*F_f_combined)/trapz(f_Hz_combined,F_f_combined));
f_E = trapz(f_Hz_combined,F_f_combined)/trapz(f_Hz_combined,f_Hz_combined.^-1.*F_f_combined);

f_ticks = 10.^(-2:1:1);

figure(fignum);clf
hold on
h_p = plot(f_p*[1 1],[1e-10 1e10],'-','Color',crimson,'linewidth',2);
h_E = plot(f_E*[1 1],[1e-10 1e10],'-.','Color',crimson,'linewidth',2);
h_m01 = plot(f_m01*[1 1],[1e-10 1e10],'--','Color',crimson,'linewidth',2);
h_m02 = plot(f_m02*[1 1],[1e-10 1e10],':','Color',crimson,'linewidth',2);
h_Pyxis = plot(f_Hz_Pyxis(ind_trim_Pyxis:end),F_f_Pyxis(ind_trim_Pyxis:end),'-','Color',teal,'linewidth',3);
h_Nortek = plot(f_Hz_EPSS(1:ind_trim_EPSS),Ff_EPSS(1:ind_trim_EPSS),'-','Color',[violet 1],'linewidth',3);
plot([0.3 0.8],1e-2*[0.3 0.8].^-4,'k--','linewidth',2)
plot([0.8 2],0.8e-2*[0.8 2].^-5,'k:','linewidth',2)
hold off
box on
xlim([1e-2 2e1])
ylim([1e-10 1e1])
ax = gca;
ax.XTick = f_ticks;
ax.XScale = 'log';
ax.YScale = 'log';
H = [h_Nortek h_Pyxis h_p h_E h_m01 h_m02];
L = {'E-PSS','direct','f_p','f_E','f_{m01}','f_{m02}'};
legend(H,L,'Location','southwest')
xlabel('f [Hz]')
ylabel('F(f) [m^2 Hz^{-1}]')
text(0.8,0.6,'$\propto f^{-4}$','HorizontalAlignment','center','FontSize',fsize,'Interpreter','latex')
text(2,1e-2,'$\propto f^{-5}$','HorizontalAlignment','center','FontSize',fsize,'Interpreter','latex')
text(8,1e-5,'k^{-2}\timesS(f)','Color',teal,'FontWeight','bold','HorizontalAlignment','center','FontSize',16)

freq_holder.f_p = f_p;
freq_holder.f_E = f_E;
freq_holder.f_m01 = f_m01;
freq_holder.f_m02 = f_m02;
freq_table = struct2table(freq_holder);
freq_table.Properties.VariableNames = {'f_p', 'f_E','f_{m01}','f_{m02}'};
disp(freq_table);

[c,~] = lindisp_with_current(2*pi*[f_p f_E f_m01 f_m02],water_depth_m,0);

celerity_holder.c_p = c(1);
celerity_holder.c_E = c(2);
celerity_holder.c_m01 = c(3);
celerity_holder.c_m02 = c(4);
freq_table = struct2table(celerity_holder);
freq_table.Properties.VariableNames = {'c_p', 'c_E','c_{m01}','c_{m02}'};
disp(freq_table);
