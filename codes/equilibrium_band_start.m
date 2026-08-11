% Frequency at which the spectrum leaves the peak flank, per run
%
% The equilibrium-range validity test fits a broken power law above the peak,
% and a fixed band start of 2.5 f_p sits above the f^-4 range in young seas:
% the spectrum there has already steepened to the saturation slope, so the
% lower segment reports saturation and the run is rejected for the wrong
% reason. This returns the frequency at which the local slope first steepens
% past slope_enter, which is where the equilibrium range actually begins
%
% Local slopes are taken on a 0.2-decade window stepped through f/f_p, so the
% answer is insensitive to individual grid points
%
% You provide:
% *f, F        - frequency [Hz] and elevation spectrum [n_f x n_run]
% *f_p         - peak frequency per run [Hz]
% *slope_enter - slope marking the end of the peak flank [default -3.5]
%
% Returns:
% *f_start - band start per run [Hz], NaN where the spectrum never steepens
%
% N. Laxague 2026
%
function f_start = equilibrium_band_start(f,F,f_p,slope_enter)

if nargin < 4 || isempty(slope_enter); slope_enter = -3.5; end

half_decade = 0.10;
ratios = 10.^(log10(1.1):0.04:log10(12));

n_run = numel(f_p);
f_start = NaN(n_run,1);

for particular_ind = 1:n_run

    if ~isfinite(f_p(particular_ind)); continue; end

    for j = 1:numel(ratios)

        f_centre = ratios(j)*f_p(particular_ind);
        in_window = f >= f_centre/10^half_decade & f <= f_centre*10^half_decade & ...
            isfinite(F(:,particular_ind)) & F(:,particular_ind) > 0;

        if sum(in_window) < 3; continue; end

        p = polyfit(log10(f(in_window)),log10(F(in_window,particular_ind)),1);

        if p(1) <= slope_enter
            f_start(particular_ind) = f_centre;
            break
        end

    end

end
