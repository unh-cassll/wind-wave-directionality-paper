function figure_style(fsize)

set(0,'DefaultAxesFontName','Liberation Serif')
set(0,'DefaultTextFontName','Liberation Serif')
set(0,'DefaultAxesFontSize',fsize)
set(0,'DefaultFigureColormap',viridis)
set(0,'DefaultAxesLineWidth',1)
set(0,'DefaultLineMarkerSize',20)
set(0,'DefaultAxesXGrid','on','DefaultAxesYGrid','on')
set(groot,'defaultAxesTickDir','out');
set(groot,'defaultAxesTickDirMode','manual');
set(groot,'defaultFigureRenderer','painters')

v_order = [7 6 3 1 5 2];
cmap = spectral(7);
set(0,'DefaultAxesColorOrder',cmap(v_order,:));
