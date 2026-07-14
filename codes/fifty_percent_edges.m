% Per-row angles of the OUTERMOST 0.5 crossing of the peak-normalized directional
% distribution D (nrow x ntheta) on each side of its peak, linearly interpolated
% between bins. Outermost (not first) so bimodal rows trace the full lobe envelope
% rather than hooking inward at a central sub-0.5 dip. theta is monotonic.
%
% theta_halfwidth = (right - left)/2, invariant to the direction frame of D.
%
function [left,right] = fifty_percent_edges(D,theta)

theta = theta(:);
nrow = size(D,1);
left = NaN(nrow,1);
right = NaN(nrow,1);
for i = 1:nrow
    d = smoothdata(D(i,:),'movmean',7,'omitnan');   % suppress angular noise
    [mx,ip] = max(d);
    if ~isfinite(mx); continue; end
    jr = find(d(ip:end)>=0.5,1,'last');             % outermost >=0.5 to the right
    if ~isempty(jr)
        b = ip+jr-1;
        if b<numel(d) && d(b)~=d(b+1); right(i) = theta(b)+(0.5-d(b))/(d(b+1)-d(b))*(theta(b+1)-theta(b)); else; right(i) = theta(b); end
    end
    jl = find(d(1:ip)>=0.5,1,'first');              % outermost >=0.5 to the left
    if ~isempty(jl)
        b = jl;
        if b>1 && d(b)~=d(b-1); left(i) = theta(b)+(0.5-d(b))/(d(b-1)-d(b))*(theta(b-1)-theta(b)); else; left(i) = theta(b); end
    end
end
