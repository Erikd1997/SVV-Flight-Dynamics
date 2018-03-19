function Mass_CG_Moment_array = Mass_CG_Moment(BEM, BEM_arm,datasheet)
[F_used1, F_used2, F_used3, Payload, Fuel0] = FuelReader(datasheet);

F_used = [F_used1;F_used2;F_used3];

Payload_arm1_14 = [131;131;170;214;214;251;251;288;288]*0.0254;
Payload_arm15 = [131;131;170;214;214;251;251;288;131]*0.0254;

Fuel_mass = Fuel0-F_used;

cg_fuel = Fuel_cg_as_function_of_fuel_mass(Fuel_mass);

Mass = sum(Payload) + BEM + Fuel_mass;

cg1_14 = (sum(Payload.*Payload_arm1_14) + Fuel_mass(1:14).*cg_fuel(1:14) + BEM*BEM_arm)./Mass(1:14);
cg15 = (sum(Payload.*Payload_arm15) + Fuel_mass(end)*cg_fuel(end) + BEM*BEM_arm)./Mass(end);

cg = [cg1_14 cg15];

Moment = Mass./cg;

Mass_CG_Moment_array = [(Mass*0.45359237)' (cg*)' Moment'];