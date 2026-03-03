%
function wind_wave_schematic_topdown(fignum,fsize)

wse_data = load('data/example_wse.mat');

x_offsets = [0.25 0.5 0.75 0.3]+0.02;
y_offsets = [0.6 0.5 0.4 0.2];
angles = [20 15 0 10];

veclength = 0.2;
gain = 1.03;

x_bot = [0 0 0 0];
y_bot = [0 0 0 0];

x_top = veclength*sind(angles)+x_bot;
y_top = veclength*cosd(angles)+y_bot;

x = [(x_bot+x_offsets)' (x_top+x_offsets)'];
y = [(y_bot+y_offsets)' (y_top+y_offsets)'];

x_bigger = [(x_bot+x_offsets)' (x_top*gain+x_offsets)'];
y_bigger = [(y_bot+y_offsets)' (y_top*gain+y_offsets)'];

cmap = flipud(spectral(7));
crimson = cmap(7,:);

cmap = [[230 120 60]/255;[75 150 115]/255;[18 60 128]/255;crimson];
labels = {'wind velocity','short waves','long waves','breakers'};

figure(fignum);clf

imagesc(wse_data.x,wse_data.y,wse_data.wse_m)
colormap('gray')
grid off
clim([-1 1]*2.25)

holder_struc = struct();

hold on

for n = 1:length(angles)

    a = annotation('arrow',x_bigger(n,:),y_bigger(n,:));
    a.LineWidth = 6;
    a.HeadSize = 26;
    a.Color = 'k';

    holder_struc(n).a = annotation('textarrow',x(n,:),y(n,:),'String',labels{n},'FontSize',fsize,'HorizontalAlignment','center','FontWeight','bold');

    holder_struc(n).a.LineWidth = 4;
    holder_struc(n).a.HeadSize = 20;
    holder_struc(n).a.Color = cmap(n,:);
    holder_struc(n).a.TextEdgeColor = 'k';
    holder_struc(n).a.TextBackgroundColor = [1 1 1 0.2];
    holder_struc(n).a.TextLineWidth = 1;

end

hold off

ax = gca;
ax.XTick = [];
ax.YTick = [];