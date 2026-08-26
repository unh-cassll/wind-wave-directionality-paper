% Wavenumber from frequency via the linear dispersion relation, no current
%
% lindisp_with_current carries surface tension and the finite-depth correction,
% so this is valid from the peak through the gravity-capillary range
%
% You provide:
% *f_Hz    - frequency [Hz], any shape
% *depth_m - water depth [m, default 15]
%
% Returns:
% *k - wavenumber [rad/m], same shape as f_Hz, NaN where f is not finite
%
% Nathan Laxague 2026
%
function k = dispersion_wavenumber(f_Hz,depth_m)

if nargin < 2 || isempty(depth_m); depth_m = 15; end

f = f_Hz(:);
k = NaN(size(f));

ok = isfinite(f) & f > 0;
if any(ok)
    c = lindisp_with_current(2*pi*f(ok),depth_m,0);
    k(ok) = 2*pi*f(ok)./c(:);
end

k = reshape(k,size(f_Hz));
