%This function calculates the pressure and density at a given height and
%outputs all atmospheric conditions into SI-units.

function [hp, Vc, Temp, p, rho] = ISA_converted(hp, Vc, TAT_K)
%Basic ISA-values:
p0     = 101325;     % [Pa]
rho0   = 1.225;      % [kg/m^3]
T0     = 288.15;     % [degK]
a0     = 340.294;    % [m/s]
g0     = 9.80665;    % [m/s^2]
gamma  = 1.401;      % [-]
R      = 287.04;     % [m^2/(degK*sec^2)]
lambda = -0.0065;  % [degK/m]

%First convert hp, IAS, and TAT into SI-units:
% h  = hp_ft.*unitsratio('meter', 'feet');
% Vc = convvel(IAS_kts, 'kts', 'm/s');
% TAT = convtemp(TAT_C, 'C', 'K');

%Find the ISA pressures and use it to find the corresponding mach numbers:
p = p0.*((1+(lambda.*hp)./T0).^(-g0./(lambda.*R)));

M = sqrt(2/(gamma-1)*((1+p0./p.*((1+(gamma-1)/(2*gamma)*rho0/p0*Vc.^2).^(gamma/(gamma-1)) -...
    1)).^((gamma-1)/gamma) - 1));

Temp = TAT_K./(1 + ((gamma-1)/2).*M.^2);

rho = p./(R.*Temp);
end