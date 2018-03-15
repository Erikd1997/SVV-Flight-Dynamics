%This function calculates the pressure and density at a given height and
%outputs all atmospheric conditions into SI-units.

function [h, V_c, T, p, rho] = ISA_converted(hp_ft, IAS_kts, TAT_C)
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
h  = hp_ft.*unitsratio('meter', 'feet');
V_c = convvel(IAS_kts, 'kts', 'm/s');
TAT = convtemp(TAT_C, 'C', 'K');

%Find the ISA pressures and use it to find the corresponding mach numbers:
p = p0.*((1+(lambda.*h)./T0).^(-g0./(lambda.*R)));

M = ((2/(gamma-1))*((1+p0./p.*(((1+((gamma-1)/(2*gamma).*rho0/p0*V_c.^2)).^...
    (gamma/(gamma-1)) - 1)).^((gamma-1)/gamma) - 1))).^0.5

T = TAT./(1 + ((gamma-1)/2).*M.^2);

rho = p./(R.*T);
end