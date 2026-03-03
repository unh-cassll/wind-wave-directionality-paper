% Given current profile and frequency-directional spectrum, computs
% wavenumber-directional spectrum via Doppler shift
%
% Nathan Laxague 2020-2024
%
function [wave_F_k_theta,k_rad_m,U_current_k,D_current_k] = directional_Doppler_shift_spectrum(Umag_m_s,Udir_deg,z_m,water_depth_m,wave_f_Hz,wave_F_f_theta,wave_theta_rad,tail_flag,k_max)

% Constants
g = 9.81;       % acceleration due to gravity in m/s^2
sigma = 0.072;  % surface tension in N/m
rho_w = 1030;   % water density in kg/m^3

% Compute U(k) from U(z) - Equation 18 from Kirby & Chen [1989], JGR
% Input arrays of size NxM (depths by frequencies)
% Output array U(k) of size Mx1 (number of frequencies/wavenumbers)
kmax_Doppler = 1000;
k_reference = logspace(-4,log10(kmax_Doppler),1024)';
wave_theta_rad = wave_theta_rad(:)';
wave_f_Hz = wave_f_Hz(:);
[s1,~] = size(wave_F_f_theta);
if length(wave_theta_rad) == s1
    wave_F_f_theta = wave_F_f_theta';
end

num_current_depths = length(z_m);

U_block = repmat(Umag_m_s(:),[1 length(k_reference)]);
D_block = repmat(Udir_deg(:),[1 length(k_reference)]);
z_block = repmat(z_m(:),[1 length(k_reference)]);

% Compute U(k) given U(z)
k_block = repmat(k_reference(:)',[num_current_depths 1]);
front_part = 2*k_block./(sinh(2*k_block*water_depth_m));
front_part = 2*k_block;%./(sinh(2*k_block*water_depth_m));

integrand_downwind = U_block.*cosd(D_block).*cosh(2*k_block.*(water_depth_m+z_block));
integrand_downwind = U_block.*cosd(D_block).*exp(2*k_block.*(z_block));
U_current_k_downwind = -1*squeeze(trapz(z_m,front_part.*integrand_downwind))';
U_current_k_downwind(isnan(U_current_k_downwind)) = 0;

integrand_crosswind = U_block.*sind(D_block).*cosh(2*k_block.*(water_depth_m+z_block));
integrand_crosswind = U_block.*sind(D_block).*exp(2*k_block.*(z_block));
U_current_k_crosswind = -1*squeeze(trapz(z_m,front_part.*integrand_crosswind))';
U_current_k_crosswind(isnan(U_current_k_crosswind)) = 0;

U_current_k = sqrt(U_current_k_downwind.^2+U_current_k_crosswind.^2);
D_current_k = atan2(U_current_k_crosswind,U_current_k_downwind);

% Compute reference wave frequency array from a test wavenumber and a
% known current velocity
f_mat = (sqrt((g*k_reference+sigma/rho_w*k_reference.^3).*tanh(k_reference*water_depth_m))+k_reference.*U_current_k.*cos(D_current_k-wave_theta_rad))/(2*pi);

% Find matching wavenumber (numerical inversion of LDR)
k_rad_m = NaN*wave_F_f_theta;
for i = 1:length(wave_f_Hz)
    for j = 1:length(wave_theta_rad)
        f_diff = abs(f_mat(:,j)-wave_f_Hz(i));
        ind = find(f_diff==min(f_diff),1,'first');
        k_rad_m(i,j) = k_reference(ind);
    end
end

% Convert from F(f,theta) to F(k,theta), conserving energy
Cg = 2*pi*wave_f_Hz./k_rad_m/2.*(1+2*k_rad_m*water_depth_m./sinh(2*k_rad_m*water_depth_m));
wave_F_k_theta_converted = wave_F_f_theta.*Cg./(pi*k_rad_m);
% Interpolates onto uniform wavenumber grid
warning('off','MATLAB:scatteredInterpolant:DupPtsAvValuesWarnId');
k_vec = linspace(max(k_rad_m(1,:)),min(k_rad_m(end,:)),1024)';
F = scatteredInterpolant(double(reshape(k_rad_m,[],1)),double(reshape(repmat(wave_theta_rad,length(wave_f_Hz),1),[],1)),double(reshape(wave_F_k_theta_converted,[],1)));
wave_F_k_theta_interp_vec = F(double(reshape(repmat(k_vec(:),1,length(wave_theta_rad)),[],1)),double(reshape(repmat(wave_theta_rad(:)',length(k_vec),1),[],1)));
wave_F_k_theta_interp_mat = reshape(wave_F_k_theta_interp_vec,length(k_vec),length(wave_theta_rad));
wave_F_k_theta_interp_mat = abs(wave_F_k_theta_interp_mat);
k_rad_m = k_vec;

if tail_flag

    % Attach k^-3 tail (including to directional spreading function)
    in_struc.k_rad_m = k_vec;
    in_struc.theta_rad = wave_theta_rad;
    in_struc.F_k_theta = wave_F_k_theta_interp_mat;
    tail_struc = pin_the_tail_on_the_spectrum(in_struc,k_max);
    k_rad_m = tail_struc.k_rad_m;
    wave_F_k_theta = tail_struc.F_k_theta_new;
    U_current_k = interp1(k_reference,U_current_k,k_vec);
    D_current_k = interp1(k_reference,D_current_k,k_vec);
else
    wave_F_k_theta = wave_F_k_theta_interp_mat;
    if Umag_m_s > 0
        U_current_k = interp1(k_reference,U_current_k,mean(k_rad_m,2,'omitnan'));
        D_current_k = interp1(k_reference,D_current_k,mean(k_rad_m,2,'omitnan'));
    else
        U_current_k = 0;
        D_current_k = 0;
    end
end