% Regression slope with standard errors that tolerate within-block correlation,
% keeping every observation rather than deflating the sample size
%
% Two estimators on the same block structure:
% * sandwich standard error for grouped observations [Liang & Zeger 1986], with
%   a finite-sample correction, on G-1 degrees of freedom
% * wild bootstrap [Wu 1986] with the null imposed, one weight per block, so
%   within-block dependence survives into each replicate
%
% Quote the bootstrap when G is of order ten; the sandwich estimator alone is
% unreliable there
%
% Inputs:
% * x,y      - predictor and response
% * block_id - acquisition block label per observation
% * n_boot   - wild bootstrap draws [default 4999]
% * weights  - 'webb' (default, six-point) or 'rademacher' (two-point)
%
% Returns:
% * slope, se_ols, se_block, G, p_ols, p_block, p_wild, CI_wild
%
% Nathan Laxague 2026
%
function s = cluster_robust_slope(x,y,block_id,n_boot,weights)

if nargin < 4 || isempty(n_boot)
    n_boot = 4999;
end

if nargin < 5 || isempty(weights)
    weights = 'webb';
end

% Six-point weights by default: Rademacher admits only 2^(G-1) distinct
% bootstrap statistics, so at G = 9 it cannot resolve p below about 1/256
webb_support = [-sqrt(1.5) -1 -sqrt(0.5) sqrt(0.5) 1 sqrt(1.5)];

inds_keep = isfinite(x) & isfinite(y) & isfinite(block_id);
x = x(inds_keep);
y = y(inds_keep);
g = block_id(inds_keep);

N = numel(y);

X = [ones(N,1) x(:)];
b = X\y(:);
u = y(:) - X*b;

XtXi = inv(X'*X);

blocks = unique(g);
G = numel(blocks);

% Sandwich meat matrix, summed over blocks
meat = zeros(2);
for i = 1:G

    inds_block = g == blocks(i);
    Xg = X(inds_block,:);
    ug = u(inds_block);

    meat = meat + (Xg'*ug)*(ug'*Xg);

end

% Finite-sample correction
c1 = G/(G-1)*(N-1)/(N-2);

V_block = c1*XtXi*meat*XtXi;

se_ols = sqrt(sum(u.^2)/(N-2)*XtXi(2,2));
se_block = sqrt(V_block(2,2));

t_ols = b(2)/se_ols;
t_block = b(2)/se_block;

% Wild bootstrap with the null slope imposed, one weight per block
Xr = ones(N,1);
br = Xr\y(:);
ur = y(:) - Xr*br;

t_boot = NaN*ones(n_boot,1);

for k = 1:n_boot

    if strcmpi(weights,'rademacher')
        draw = 2*(rand(G,1) > 0.5) - 1;      % two-point, one per block
    else
        draw = webb_support(randi(6,G,1))';  % six-point, one per block
    end

    w = ones(N,1);
    for i = 1:G
        w(g == blocks(i)) = draw(i);
    end

    yb = Xr*br + ur.*w;
    bb = X\yb;
    ub = yb - X*bb;

    meat_b = zeros(2);
    for i = 1:G
        inds_block = g == blocks(i);
        meat_b = meat_b + (X(inds_block,:)'*ub(inds_block))*(ub(inds_block)'*X(inds_block,:));
    end

    Vb = c1*XtXi*meat_b*XtXi;

    t_boot(k) = bb(2)/sqrt(Vb(2,2));

end

s.slope = b(2);
s.se_ols = se_ols;
s.se_block = se_block;
s.G = G;
s.p_ols = 2*(1-tcdf(abs(t_ols),N-2));
s.p_block = 2*(1-tcdf(abs(t_block),G-1));
s.p_wild = (1+sum(abs(t_boot) >= abs(t_block)))/(n_boot+1);
s.CI_wild = b(2) + [-1 1]*prctile(abs(t_boot),95)*se_block;
