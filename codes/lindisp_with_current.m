function [c,cg] = lindisp_with_current(omega,h,current_m_s)
%
% Nathan Laxague's brute-force attempt at numerically inverting the LDR
%
% You provide:
% *omega - wave radian frequency in rad/s
% *h - water depth in m
% *current_m_s - current in m/s
%
% Doppler shift calculation assumes waves and currents are all aligned
%
% 2025

omega = omega(:);
h = h(:);
current_m_s = current_m_s(:);

g = 9.806;
rho_w = 1020;
sigma = 0.072;

k_vec = logspace(-3,log10(300),100);
omega_disp = sqrt((g*k_vec+sigma/rho_w*k_vec.^3).*tanh(k_vec*h))+k_vec*current_m_s;

k_from_omega = fit(omega_disp(:),k_vec(:),'spline');
k = k_from_omega(omega);

c = omega./k;
cg = c/2.*(1+(2*k*h)./(sinh(2*k*h)));