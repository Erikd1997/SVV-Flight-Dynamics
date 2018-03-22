%This function calculates the pressure and density at a given height and
%outputs all atmospheric conditions into SI-units.

<<<<<<< HEAD
function [hp, Vc, Temp, p, rho] = ISA_converted(hp, Vc, TAT_K)
=======
function [h, V_c, T, p, rho] = ISA_converted(h, V_c, TAT)
>>>>>>> 1513f994b678d122aca6dedc06d43c1968494faa
%Basic ISA-values:
p0     = 101325;     % [Pa]
rho0   = 1.225;      % [kg/m^3]
T0     = 288.15;     % [degK]
g0     = 9.80665;    % [m/s^2]
gamma  = 1.401;      % [-]
R      = 287.04;     % [m^2/(degK*sec^2)]
lambda = -0.0065;  % [degK/m]

<<<<<<< HEAD
%First convert hp, IAS, and TAT into SI-units:
% h  = hp_ft.*unitsratio('meter', 'feet');
% Vc = convvel(IAS_kts, 'kts', 'm/s');
% TAT = convtemp(TAT_C, 'C', 'K');

=======
>>>>>>> 1513f994b678d122aca6dedc06d43c1968494faa
%Find the ISA pressures and use it to find the corresponding mach numbers:
p = p0.*((1+(lambda.*hp)./T0).^(-g0./(lambda.*R)));

<<<<<<< HEAD
M = sqrt(2/(gamma-1)*((1+p0./p.*((1+(gamma-1)/(2*gamma)*rho0/p0*Vc.^2).^(gamma/(gamma-1)) -...
    1)).^((gamma-1)/gamma) - 1));

Temp = TAT_K./(1 + ((gamma-1)/2).*M.^2);
=======
M = sqrt(2/(gamma-1)*((1+p0./p.*((1+(gamma-1)/(2*gamma)*rho0/p0*V_c.^2).^(gamma/(gamma-1)) -...
    1)).^((gamma-1)/gamma) - 1));

T = TAT./(1 + ((gamma-1)/2).*M.^2);
>>>>>>> 1513f994b678d122aca6dedc06d43c1968494faa

rho = p./(R.*Temp);
end