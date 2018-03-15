function [V_r_eq, V_t, T, M, rho] = reduced_eq_airspeed(V_c, hp, TAT, W, W_s)
%This function transforms the calibrated airspeeds (which is the IAS from
%the flight test) into the reduced, equivalent airspeeds. All input
%variables are expected to be in SI-units. V_c, hp, TAT, and W should be
%one-dimensional arrays of the same size.

%First some basic ISA-values:
p0     = 101325;     % [Pa]
rho0   = 1.225;      % [kg/m^3]
T0     = 288.15;     % [degK]
a0     = 340.294;    % [m/s]
g0     = 9.80665;    % [m/s^2]
gamma  = 1.401;      % [-]
R      = 287.04;     % [m^2/(degK*sec^2)]
lambda = -0.0065;    % [degK/m]

%First find the ISA pressures and use it to find the corresponding mach 
%numbers: 
p = p0*(1+lambda*hp/T0).^(-g0/(lambda*R));

M = sqrt(2/(gamma-1)*((1+p0./p.*((1+(gamma-1)/(2*gamma)*rho0/p0*V_c.^2).^(gamma/(gamma-1)) -...
    1)).^((gamma-1)/gamma) - 1));

%Now find the static air temperatures and corresponding speeds of sound,
%along with the air densities.=:
T = TAT./(1 + ((gamma-1)/2).*M.^2);

a = (gamma*R.*T).^0.5;

rho = p./(R.*T);

%Now calculate the true and equivalent airspeeds:
V_t = M.*a;

V_eq = V_t.*sqrt(rho./rho0);

%Now convert into reduced equivalent airspeeds:
V_r_eq = V_eq.*sqrt(W_s./W);

