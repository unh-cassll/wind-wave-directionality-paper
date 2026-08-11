function wind_speed_direction_rose(fignum,fsize)

compilation_nc_name = 'data/ASIT2019_supporting_data_compilation.nc';

% Wind forcing per the source selected in the aa_step01 preamble; the hybrid
% *_best variables exist only on the production grid, so this falls back. The
% ustar gate is NOT applied: this figure documents the wind climatology, so
% filtering out the light-wind sectors would misrepresent the record.
wind = wind_forcing(compilation_nc_name,false);
COARE_U10 = wind.U10;
wdir_deg_all = mod(ncread(compilation_nc_name,'COARE_winddir'),360);

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
