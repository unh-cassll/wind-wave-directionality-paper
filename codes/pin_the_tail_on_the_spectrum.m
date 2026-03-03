% Given wavenumber water surface elevation spectrum, pins on a piecewise
% tail:
% * k^-2.5 in the equilibrium range
% * k^-3 in the saturation range
%
% Transition from equilibrium to saturation is defined as k^3F(k) = 8e-3
%
% Nathan Laxague 2024
%
function out_struc = pin_the_tail_on_the_spectrum(in_struc,k_max)

% form of 'in_struc':
% * wavenumber in rad/m -> obs_wave_spectrum.k_rad_m
% ** Mx1 array
% * direction going-to in rad -> obs_wave_spectrum.theta_rad
% ** 1xL array
% * spectral energy density in m^4 -> obs_wave_spectrum.F_k_theta
% ** MxL array

degree_of_saturation = 8e-3;

k_rad_m = in_struc.k_rad_m;
theta_rad = in_struc.theta_rad;
F_k_theta = in_struc.F_k_theta;

[k_rad_m,inds_unique] = unique(k_rad_m);
F_k_theta = F_k_theta(inds_unique,:);

F_k = trapz(theta_rad,k_rad_m.*F_k_theta,2);

dk = min(k_rad_m);
k_high_old = max(k_rad_m);

D_k = F_k_theta./(F_k./k_rad_m);

B_k = k_rad_m.^3.*F_k;

ind_start = floor(0.75*length(k_rad_m));
inds_consider = ind_start:length(k_rad_m);

D_k_add = mean(D_k(inds_consider,:),1);

p = polyfit(log10(k_rad_m(inds_consider)),log10(B_k(inds_consider)),1);

if p(1) >= 0

    B_high = 10.^(p(1)*log10(k_high_old)+p(2));

    if B_high < degree_of_saturation

        k_transition = 10.^((log10(degree_of_saturation)-log10(B_high))/p(1)+log10(k_high_old));

        if k_transition < k_high_old - dk

            eq_k = [k_high_old+dk;k_transition];
            sat_k = [k_transition+dk;k_max];

            k_add = [eq_k; sat_k];
            B_k_add = [B_high; degree_of_saturation; degree_of_saturation*[1; 1]];

            k_new = logspace(log10(dk),log10(k_max),2^nextpow2(length(k_rad_m)))';
            B_k_new = 10.^interp1(log10([k_rad_m; k_add]),log10([B_k; B_k_add]),log10(k_new));

        else

            k_add = [k_high_old+dk; k_max];
            B_k_add = [B_high; degree_of_saturation];

            k_new = logspace(log10(dk),log10(k_max),2^nextpow2(length(k_rad_m)))';
            B_k_new = 10.^interp1(log10([k_rad_m; k_add]),log10([B_k; B_k_add]),log10(k_new));

        end

        D_k_new = interp2(theta_rad(:)',[k_rad_m; k_add],[D_k; ones(length(k_add),1)*D_k_add],theta_rad(:)',k_new);

    else

        k_new = logspace(log10(dk),log10(k_max),2^nextpow2(length(k_rad_m)))';

        B_k(end-16:end) = mean(B_k(end-16:end),'omitnan');
        D_k(end-16:end,:) = repmat(mean(D_k(end-16:end,:),'omitnan'),[17 1]);

        B_k_new = 10.^interp1(log10(k_rad_m),log10(B_k),log10(k_new),'linear','extrap');
        D_k_new = 10.^interp1(log10(k_rad_m),log10(D_k),log10(k_new),'linear','extrap');

    end

else

    k_add = [k_high_old + dk; k_max];
    B_k_add = 10.^(p(1)*log10(k_add)+p(2));

    k_new = logspace(log10(dk),log10(k_max),2^nextpow2(length(k_rad_m)))';
    B_k_new = 10.^interp1(log10([k_rad_m; k_add]),log10([B_k; B_k_add]),log10(k_new));

    D_k_new = 10.^interp2(theta_rad(:)',log10([k_rad_m; k_add]),log10([D_k; ones(length(k_add),1)*D_k_add]),theta_rad(:)',log10(k_new));

end

F_k_new = k_new.^-3.*B_k_new;
F_k_theta_new = k_new.^-1.*F_k_new.*D_k_new;
F_k_theta_new(isnan(F_k_theta_new)) = 0;

out_struc.k_rad_m = k_new;
out_struc.F_k_theta_new = F_k_theta_new;


