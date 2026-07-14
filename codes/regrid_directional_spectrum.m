% Resample a directional spectrum from one direction grid onto another.
% Periodic (2*pi) resampling via 3x tiling and per-frequency interp1 over
% finite samples; target points outside the source support are left NaN.
%
% S:         nf x ndir_src spectral density on directions theta_src (rad)
% theta_src: source directions (rad)
% theta_tgt: target directions (rad)
% S_out:     nf x ndir_tgt density resampled onto theta_tgt
%
% N. Laxague 2026
%
function S_out = regrid_directional_spectrum(S,theta_src,theta_tgt)

theta_src = theta_src(:)';
theta_tgt = theta_tgt(:)';

% Orient S as nf x ndir (direction along columns)
[s1,s2] = size(S);
if length(theta_src) == s1 && length(theta_src) ~= s2
    S = S';
end

[nf,~] = size(S);
ntgt = length(theta_tgt);

% Periodic tiling of source directions across +/- 2*pi
theta_src_big = [theta_src-2*pi theta_src theta_src+2*pi];

S_out = NaN*ones(nf,ntgt);

for i = 1:nf

    row_big = [S(i,:) S(i,:) S(i,:)];
    finite_mask = isfinite(theta_src_big) & isfinite(row_big);

    if nnz(finite_mask) < 2
        continue
    end

    xs = theta_src_big(finite_mask);
    ys = row_big(finite_mask);
    [xs,order] = sort(xs);
    ys = ys(order);
    [xs,ia] = unique(xs);
    ys = ys(ia);

    S_out(i,:) = interp1(xs,ys,theta_tgt,'linear',NaN);

end
