
function generate_MV_ASIT_map(fignum,fsize)

nc_name = 'ASIT2019_supporting_environmental_observations.nc';
LON = ncread(nc_name,'bathy_lon_deg');
LAT = ncread(nc_name,'bathy_lat_deg');
DEPTH = ncread(nc_name,'bathy_depth_m');

depthlevels = -19:2:-9;

asit_position = [41.333150 -70.573135];
asit_box_lon = asit_position(2) + [-1 1]*1/60;
asit_box_lat = asit_position(1) + [-1 1]*1/60;

lonlims = [-71.0 -69.9]+0.01;
latlims = [41 42] + 1/6;
depth_lims = [depthlevels(1) depthlevels(end)];

markercolor = [200 0 0]/255;
patchcolor = [253 190 110]/255;

cmap = deep(length(depthlevels)+1);
cmap = cmap(1:end-1,:);

figure(fignum);clf
tlayout = tiledlayout(1,2);
ax_struc = struct();

nexttile()

m_proj('miller','long',lonlims,'lat',latlims);
m_gshhs_f('patch',patchcolor);
m_gshhs_f('speckle','color','k');
hold on
m_plot(asit_position(2),asit_position(1),'o','markeredgecolor',markercolor,'markerfacecolor',markercolor,'markersize',3)
m_plot([asit_box_lon(1) asit_box_lon(2) asit_box_lon(2) asit_box_lon(1) asit_box_lon(1)],[asit_box_lat(2) asit_box_lat(2) asit_box_lat(1) asit_box_lat(1) asit_box_lat(2)],'Color',markercolor,'linewidth',2)
hold off
m_grid('box','fancy','fontsize',fsize,'linewidth',2,'tickdir','out','xaxisloc','bottom','yaxisloc','left');
text(-0.1,0.95,'(a)','Units','Normalized','FontSize',fsize)
ax_struc(1).ax = gca;

nexttile(2)

m_proj('miller','long',asit_box_lon,'lat',asit_box_lat);
m_gshhs_f('patch',patchcolor);
m_gshhs_f('speckle','color','k');
hold on
[C,h] = m_contour(LON,LAT,DEPTH,depthlevels,'ShowText','on','LineWidth',3,"LabelFormat","%0.0f m");
m_plot(asit_position(2),asit_position(1),'o','markeredgecolor',markercolor,'markerfacecolor',markercolor,'markersize',5)
hold off
m_grid('box','fancy','fontsize',fsize,'linewidth',2,'tickdir','out','xaxisloc','bottom','yaxisloc','left');
clabel(C,h,'FontName','Latin Modern Math','FontSize',fsize*0.75)
text(-0.1,0.95,'(b)','Units','Normalized','FontSize',fsize)

colormap(cmap)
clim(depth_lims)

tlayout.TileSpacing = 'tight';

% top line
x = [0.237 0.545];
y = [0.257 0.925];
annotation('line',x,y,'Color',markercolor,'LineWidth',2);

% bottom line
x = [0.237 0.545];
y = [0.23 0.115];
annotation('line',x,y,'Color',markercolor,'LineWidth',2);


