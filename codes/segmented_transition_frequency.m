% Equilibrium-saturation transition frequency from a fitted broken power law
%
% Model, continuous at the break, on x = log10(f), y = log10(E):
% y = a + b1*(x-xb) + (b2-b1)*max(0,x-xb)
%
% Unlike a threshold-crossing scan, this tests whether two slopes are justified
% at all, so spectra with no f^-4 range are rejected rather than assigned a
% transition on the peak flank
%
% You provide:
% *f, E   - frequency [Hz] and elevation spectrum, vectors
% *f_lo   - lower edge of the fit band [Hz]
% *f_hi   - upper edge of the fit band [Hz]
% *slope_lims - acceptable equilibrium-side slope range [default [-4.8 -3.2]]
%
% Returns:
% *f_break   - breakpoint frequency [Hz], NaN if rejected
% *ci        - approximate 95% profile interval on f_break
% *slope_eq  - fitted slope below the break (expected near -4)
% *slope_sat - fitted slope above the break (expected near -5)
% *p_two     - F-test p-value, two slopes against one
% *rmse      - residual RMS of the two-slope fit [log10 units]
% *n         - points in the fit band
% *accepted  - passed every acceptance test
% *reason    - why it was rejected, '' if accepted
%
% N. Laxague 2026
%
function out = segmented_transition_frequency(f,E,f_lo,f_hi,slope_lims)

if nargin < 5 || isempty(slope_lims); slope_lims = [-4.8 -3.2]; end

out = struct('f_break',NaN,'ci',[NaN NaN],'slope_eq',NaN,'slope_sat',NaN, ...
    'p_two',NaN,'rmse',NaN,'n',0,'accepted',false,'reason','');

band = f >= f_lo & f <= f_hi & isfinite(f) & isfinite(E) & E > 0;
x = log10(f(band));
y = log10(E(band));
[x,ord] = sort(x(:)); y = y(ord);
n = numel(x);
out.n = n;

if n < 12
    out.reason = 'too few points in the fit band';
    return
end

% Candidate breakpoints: interior only, leaving at least 4 points per segment
cand = x(5:end-4);
sse = NaN(numel(cand),1);
pars = NaN(numel(cand),3);
for i = 1:numel(cand)
    X = [ones(n,1) (x-cand(i)) max(0,x-cand(i))];
    b = X\y;
    r = y - X*b;
    sse(i) = sum(r.^2);
    pars(i,:) = b(:)';
end

[sse_min,ib] = min(sse);
xb = cand(ib);
b = pars(ib,:);

out.slope_eq = b(2);
out.slope_sat = b(2)+b(3);
out.f_break = 10^xb;
out.rmse = sqrt(sse_min/(n-4));

% Two slopes against one
X1 = [ones(n,1) x];
r1 = y - X1*(X1\y);
sse1 = sum(r1.^2);
F = ((sse1-sse_min)/2)/(sse_min/(n-4));
out.p_two = 1-fcdf(F,2,n-4);

% Approximate profile interval: breakpoints whose SSE is not significantly worse
thresh = sse_min*(1 + finv(0.95,1,n-4)/(n-4));
in_ci = cand(sse <= thresh);
out.ci = [10^min(in_ci) 10^max(in_ci)];

% Acceptance
if out.p_two > 0.05
    out.reason = 'single slope fits as well';
elseif out.slope_eq < slope_lims(1) || out.slope_eq > slope_lims(2)
    out.reason = sprintf('equilibrium slope %.2f outside [%.1f %.1f]', ...
        out.slope_eq,slope_lims(1),slope_lims(2));
elseif out.slope_sat >= out.slope_eq
    out.reason = 'spectrum does not steepen across the break';
elseif ib <= 2 || ib >= numel(cand)-1
    out.reason = 'breakpoint at the edge of the fit band';
else
    out.accepted = true;
end

if ~out.accepted
    out.f_break = NaN;
    out.ci = [NaN NaN];
end
