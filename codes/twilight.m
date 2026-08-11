% Cyclic colormap in the style of matplotlib's twilight, for signed angles that
% wrap: black in the middle, saturated blue and red on the two flanks, white at
% both ends, so -180 and +180 meet without a seam.
%
% Built in CIE L*a*b* rather than tabulated, so the ramp is perceptually even by
% construction: lightness rises linearly from black at the centre to white at
% either end, chroma follows a half-sine that vanishes at both achromatic
% extremes, and hue is fixed at one value per flank (blue below, red above).
%
% Used with clim([-180 180]), a value maps to:
%   0 deg    black       +/-90 deg  blue (-) / red (+)      +/-180 deg  white
%
% You provide:
% *n - number of colors [default 256]
%
% Returns:
% *cm - n-by-3 sRGB colormap, first row = -180 deg, last row = +180 deg
%
% N. Laxague 2026
%
function cm = twilight(n)

if nargin < 1 || isempty(n); n = 256; end

s = linspace(-1,1,n)';   % normalized angle: -1 = -180 deg, 0 = 0 deg, +1 = +180 deg

L_min = 8;               % centre: near-black
L_max = 92;              % ends: near-white
C_max = 45;              % flanks, at +/-90 deg
h_blue = 285;            % hue of the negative flank [deg]
h_red = 20;              % hue of the positive flank [deg]

L = L_min + (L_max-L_min)*abs(s);
C = C_max*sin(pi*abs(s));

h = h_blue*ones(size(s));
h(s > 0) = h_red;

cm = lab2srgb([L C.*cosd(h) C.*sind(h)]);


% CIE L*a*b* (D65) to gamma-encoded sRGB, clipped to gamut
function rgb = lab2srgb(lab)

L = lab(:,1); a = lab(:,2); b = lab(:,3);

fy = (L+16)/116;
fx = fy + a/500;
fz = fy - b/200;

finv = @(t) (t > 6/29).*t.^3 + (t <= 6/29).*(3*(6/29)^2*(t - 4/29));

white = [0.95047 1 1.08883];
X = white(1)*finv(fx);
Y = white(2)*finv(fy);
Z = white(3)*finv(fz);

M = [ 3.2406 -1.5372 -0.4986
     -0.9689  1.8758  0.0415
      0.0557 -0.2040  1.0570];

lin = [X Y Z]*M';

rgb = (lin <= 0.0031308).*(12.92*lin) + (lin > 0.0031308).*(1.055*max(lin,0).^(1/2.4) - 0.055);
rgb = min(max(rgb,0),1);
