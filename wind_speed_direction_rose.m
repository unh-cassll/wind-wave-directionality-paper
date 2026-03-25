function wind_speed_direction_rose(fignum,fsize)

supporting_data_nc_name = 'data/ASIT2019_supporting_data_compilation.nc';

COARE_U10 = ncread(supporting_data_nc_name,'COARE_U10');
wdir_deg_all = mod(ncread(supporting_data_nc_name,'COARE_winddir'),360);

Options = struct();
Options.AngleNorth = 0;
Options.AngleEast = 90;
Options.axesfontname = 'Liberation Serif';
Options.textfontname = 'Liberation Serif';
Options.frequencyFontSize = fsize*3/4;
Options.axesfontsize = fsize;
Options.titlefontsize = fsize;
Options.legendfontsize = fsize*3/4;
Options.TitleString = '';
Options.legendvariable = 'U_{10}';
Options.lablegend = [];
Options.legendposition = 'southeast';
% Options.NSpeeds = 7;

U_high = 13;
U_low = 1;
inds_keep = COARE_U10 >= U_low & COARE_U10 <= U_high;

dU = 2;

Options.vwinds = U_low:dU:11;

% Options.cmap = flipud(spectral);
Options.cmap = (magma);

figure(fignum);clf
ax = gca;
Options.axes = ax;
WindRose(wdir_deg_all(inds_keep),COARE_U10(inds_keep),Options);
