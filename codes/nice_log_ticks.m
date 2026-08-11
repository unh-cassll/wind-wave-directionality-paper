% Readable integer tick values for a logarithmic axis spanning lims
%
% Uses the 1-2-3-5-7 decade progression and always includes the endpoints, so
% ticks stay round rather than falling on whatever bin edges are in play
%
% You provide:
% *lims - [min max], both positive
%
% Returns:
% *t - integer tick values, ascending
%
% N. Laxague 2026
%
function t = nice_log_ticks(lims)

base = [1 2 3 5 7];
decades = floor(log10(lims(1))):ceil(log10(lims(2)));

cand = reshape(base(:)*10.^decades(:)',1,[]);
cand = cand(cand > lims(1) & cand < lims(2));

t = unique(round([lims(1) cand lims(2)]));
