function stacked_timeseries(fignum,fsize)

supporting_data_nc_name = 'data/ASIT2019_supporting_data_compilation.nc';

cmap = flipud(spectral(21));
violet = cmap(1,:);
bluish = cmap(2,:)*1.1;
teal = cmap(4,:);
greenish = cmap(6,:)*0.7;
goldenrod = cmap(16,:)*0.8;
crimson = cmap(21,:);

DN_start = datenum(2019,10,16,0,0,0);
time = double(ncread(supporting_data_nc_name,'time'));
DN = DN_start + time/60/24;
DTime = datetime(datevec(DN),'TimeZone','UTC');

COARE_U10 = ncread(supporting_data_nc_name,'COARE_U10');
COARE_ustar_m_s = ncread(supporting_data_nc_name,'COARE_ustar');
EC_ustar_m_s = ncread(supporting_data_nc_name,'EC_ustar');

Tpeak = ncread(supporting_data_nc_name,'COARE_peak_wave_period');
Tmean = ncread(supporting_data_nc_name,'COARE_mean_wave_period');

f_Hz = ncread(supporting_data_nc_name,'frequency_omni');
F_f_m2_Hz = ncread(supporting_data_nc_name,'Omni_F_f');
Hm0 = 4*sqrt(trapz(f_Hz,F_f_m2_Hz));

DT_lim = [datetime(2019,10,7,'TimeZone','UTC') datetime(2020,2,14,'TimeZone','UTC')];

lw = 1;

num_stacks = 4;
num_wide = 4;

nbins = 40;

U10_lims = [0 20];
ustar_lims = [5e-3 5e0];
Hs_lims = [0 4];
T_lims = [0 15];
Plims = [0 1.25];

text_x_left = 0.03;
text_x_right = 0.92;
text_y = 0.9;

figure(fignum);clf
tlayout = tiledlayout(num_stacks,num_wide);
ax_struc = struct();

% time series

nexttile((1-1)*num_wide+1,[1 num_wide-1])
plot(DTime,COARE_U10,'Color',bluish,'linewidth',lw)
xlim(DT_lim)
ylim(U10_lims)
ax_struc((1-1)*num_wide+1).ax = gca;
ylabel('$\mathrm{U_{10}\ [m\ s^{-1}]}$','Interpreter','LaTeX')
text(text_x_left,text_y,'(a)','Units','Normalized','FontSize',fsize,'HorizontalAlignment','center')

nexttile((2-1)*num_wide+1,[1 num_wide-1])
hold on
h_EC = plot(DTime,EC_ustar_m_s,'-','Color',violet,'markersize',3,'linewidth',lw);
h_COARE = plot(DTime,COARE_ustar_m_s,'-','Color',teal,'markersize',3,'linewidth',lw);
hold off
box on
xlim(DT_lim)
ylim(ustar_lims)
ax_struc((2-1)*num_wide+1).ax = gca;
ax_struc((2-1)*num_wide+1).ax.YScale='log';
ax_struc((2-1)*num_wide+1).ax.YTick = 10.^(-2:1:0);
ax_struc((2-1)*num_wide+1).ax.YTickLabel = {'0.01','0.1','1'};
ylabel('$\mathrm{u_{*}\ [m\ s^{-1}]}$','Interpreter','LaTeX')

text(text_x_left,text_y,'(b)','Units','Normalized','FontSize',fsize,'HorizontalAlignment','center')

nexttile((3-1)*num_wide+1,[1 num_wide-1])
plot(DTime,Hm0,'-','Color',crimson,'markersize',3,'linewidth',lw)
xlim(DT_lim)
ylim(Hs_lims)
ax_struc((3-1)*num_wide+1).ax = gca;
ylabel('$\mathrm{H_s\ [m]}$','Interpreter','LaTeX')
text(text_x_left,text_y,'(c)','Units','Normalized','FontSize',fsize,'HorizontalAlignment','center')

nexttile((4-1)*num_wide+1,[1 num_wide-1])
hold on
h_peak = plot(DTime,Tpeak,'-','Color',goldenrod,'markersize',3,'linewidth',lw);
h_avg = plot(DTime,Tmean,'-','Color',greenish,'markersize',3,'linewidth',lw);
hold off
box on
xlim(DT_lim)
ylim(T_lims)
ax_struc((4-1)*num_wide+1).ax = gca;
ylabel('$\mathrm{T\ [s]}$','Interpreter','LaTeX')

text(text_x_left,text_y,'(d)','Units','Normalized','FontSize',fsize,'HorizontalAlignment','center')

% histograms
[counts_U10,centers_U10] = hist(COARE_U10,nbins);
P_U10 = counts_U10/max(counts_U10);

nexttile(1*num_wide)
plot(P_U10,centers_U10,'-','Color',bluish,'LineWidth',lw*1.5)
xlim(Plims)
ylim(U10_lims)
ax_struc(1*num_wide).ax = gca;
text(text_x_right,text_y,'(e)','Units','Normalized','FontSize',fsize,'HorizontalAlignment','center')

ustar_m_s_EC = EC_ustar_m_s;
ustar_m_s_EC(ustar_m_s_EC>1) = NaN;

[counts_ustar_EC,centers_ustar_EC] = hist(ustar_m_s_EC,nbins);
P_ustar_EC = counts_ustar_EC/max(counts_ustar_EC);
[counts_ustar_bulk,centers_ustar_bulk] = hist(COARE_ustar_m_s,nbins);
P_ustar_bulk = counts_ustar_bulk/max(counts_ustar_bulk);

nexttile(2*num_wide)
hold on
h_EC = plot(P_ustar_EC,centers_ustar_EC,'-','Color',violet,'LineWidth',lw*1.5);
h_COARE = plot(P_ustar_bulk,centers_ustar_bulk,'-','Color',teal,'LineWidth',lw*1.5);
hold off
box on
xlim(Plims)
ylim(ustar_lims)
ax_struc(2*num_wide).ax = gca;
ax_struc(2*num_wide).ax.YScale='log';

text(0.9,0.4,'EC','Color',violet,'FontSize',fsize,'HorizontalAlignment','center')
text(0.85,0.015,'COARE 3.5','Color',teal,'FontSize',fsize,'HorizontalAlignment','center')
text(text_x_right,text_y,'(f)','Units','Normalized','FontSize',fsize,'HorizontalAlignment','center')

[counts_Hs,centers_Hs] = hist(Hm0,nbins);
P_Hs = counts_Hs/max(counts_Hs);

nexttile(3*num_wide)
plot(P_Hs,centers_Hs,'-','Color',crimson,'LineWidth',lw*1.5)
xlim(Plims)
ylim(Hs_lims)
ax_struc(3*num_wide).ax = gca;
text(text_x_right,text_y,'(g)','Units','Normalized','FontSize',fsize,'HorizontalAlignment','center')

[counts_Tpeak,centers_Tpeak] = hist(Tpeak,nbins);
P_Tpeak = counts_Tpeak/max(counts_Tpeak);
[counts_Tave,centers_Tave] = hist(Tmean,nbins);
P_Tave = counts_Tave/max(counts_Tave);

nexttile(4*num_wide)
hold on
h_peak = plot(P_Tpeak,centers_Tpeak,'-','Color',goldenrod,'LineWidth',lw*1.5);
h_avg = plot(P_Tave,centers_Tave,'-','Color',greenish,'LineWidth',lw*1.5);
hold off
box on
xlim(Plims)
ylim(T_lims)
ax_struc(4*num_wide).ax = gca;
xlabel('$\mathrm{P(X)\ [normalized]}$','Interpreter','latex')

text(0.75,11.5,'$\mathrm{T_p}$','Color',goldenrod,'FontSize',fsize,'HorizontalAlignment','center','Interpreter','latex')
text(1.1,2,'$\mathrm{T_m}$','Color',greenish,'FontSize',fsize,'HorizontalAlignment','center','Interpreter','latex')
text(text_x_right,text_y,'(h)','Units','Normalized','FontSize',fsize,'HorizontalAlignment','center')

tile_cleaner(ax_struc,tlayout)
