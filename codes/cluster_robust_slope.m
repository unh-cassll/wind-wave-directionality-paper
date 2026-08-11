% Regression slope with standard errors that tolerate within-block correlation,
% keeping every observation rather than deflating the sample size
%
% Two estimators: a CR1 cluster-robust standard error on G-1 degrees of freedom,
% and a wild cluster bootstrap with the null imposed. Quote the bootstrap when
% G is of order ten, where CR1 is unreliable
%
% You provide:
% *x,y      - predictor and response
% *block_id - acquisition block label per observation
% *n_boot   - wild cluster bootstrap draws [default 4999]
% *weights  - 'webb' (default) or 'rademacher'
%
% Returns:
%   slope, se_ols, se_CR, G, p_ols, p_CR, p_wild, CI_wild
%
% N. Laxague 2026
%
function s = cluster_robust_slope(x,y,block_id,n_boot,weights)

if nargin < 4 || isempty(n_boot);  n_boot = 4999; end
if nargin < 5 || isempty(weights); weights = 'webb'; end

% Webb six-point weights by default: Rademacher admits only 2^(G-1) distinct
% bootstrap statistics, so at G = 9 it cannot resolve p below about 1/256
webb_support = [-sqrt(1.5) -1 -sqrt(0.5) sqrt(0.5) 1 sqrt(1.5)];

ok = isfinite(x) & isfinite(y) & isfinite(block_id);
x = x(ok); y = y(ok); g = block_id(ok);

N = numel(y);
X = [ones(N,1) x(:)];
b = X\y(:);
u = y(:) - X*b;

XtXi = inv(X'*X);
blocks = unique(g);
G = numel(blocks);

% Cluster-robust meat matrix with the CR1 finite-sample correction
meat = zeros(2);
for i = 1:G
    j = g == blocks(i);
    Xg = X(j,:);
    ug = u(j);
    meat = meat + (Xg'*ug)*(ug'*Xg);
end
c1 = G/(G-1)*(N-1)/(N-2);
V_CR = c1*XtXi*meat*XtXi;

se_ols = sqrt(sum(u.^2)/(N-2)*XtXi(2,2));
se_CR = sqrt(V_CR(2,2));

t_ols = b(2)/se_ols;
t_CR = b(2)/se_CR;

% Wild cluster bootstrap with the null slope imposed
Xr = ones(N,1);
br = Xr\y(:);
ur = y(:) - Xr*br;
t_boot = NaN(n_boot,1);
for k = 1:n_boot
    w = ones(N,1);
    if strcmpi(weights,'rademacher')
        draw = 2*(rand(G,1) > 0.5) - 1;          % two-point, one per block
    else
        draw = webb_support(randi(6,G,1))';      % six-point, one per block
    end
    for i = 1:G
        w(g == blocks(i)) = draw(i);
    end
    yb = Xr*br + ur.*w;
    bb = X\yb;
    ub = yb - X*bb;
    meat_b = zeros(2);
    for i = 1:G
        j = g == blocks(i);
        meat_b = meat_b + (X(j,:)'*ub(j))*(ub(j)'*X(j,:));
    end
    Vb = c1*XtXi*meat_b*XtXi;
    t_boot(k) = bb(2)/sqrt(Vb(2,2));
end

s.slope = b(2);
s.se_ols = se_ols;
s.se_CR = se_CR;
s.G = G;
s.p_ols = 2*(1-tcdf(abs(t_ols),N-2));
s.p_CR = 2*(1-tcdf(abs(t_CR),G-1));
s.p_wild = (1+sum(abs(t_boot) >= abs(t_CR)))/(n_boot+1);
s.CI_wild = b(2) + [-1 1]*prctile(abs(t_boot),95)*se_CR;
