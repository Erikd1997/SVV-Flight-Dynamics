function Mass_CG_Moment_array = Mass_CG_Moment(BEM, BEM_arm,datasheet)
[F_used1, F_used2, F_used3, Payload, Fuel0] = FuelReader(datasheet);

BEM = BEM*0.45359237;                                           %[kg]
BEM_arm = BEM_arm*0.0254;                                       %[m]

F_used = [F_used1;F_used2;F_used3];                             %[lbs]

Payload_arm1_14 = [131;131;170;214;214;251;251;288;288]*0.0254; %[m]
Payload_arm15 = [131;131;170;214;214;251;251;288;131]*0.0254;   %[m]

Fuel_mass = Fuel0-F_used;                                       %[lbs]

cg_fuel = 0.0254*Fuel_cg_as_function_of_fuel_mass(Fuel_mass);   %[m]

Mass = sum(Payload) + BEM + Fuel_mass*0.45359237;               %[kg]

Fuel_mass = Fuel_mass*0.45359237;                               %[kg]

cg1_14 = (sum(Payload.*Payload_arm1_14) + Fuel_mass(1:14).*cg_fuel(1:14)' + BEM*BEM_arm)./Mass(1:14); %[m]
cg15 = (sum(Payload.*Payload_arm15) + Fuel_mass(end)*cg_fuel(end) + BEM*BEM_arm)./Mass(end);         %[m]

cg = [cg1_14; cg15];                                             %[m]

Moment = Mass.*9.81.*cg;                                          %[N*m]

Mass_CG_Moment_array = [Mass cg Moment];
end