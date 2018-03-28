%This program handles the calculation of Cm_alpha and Cm_delta from our
%flight test data.

%First, import the data from the excel file and convert them into SI-units:
datasheet = 'Post_Flight_Datasheet_Flight_2_DD_6_3_2018_for_test.xlsx';

%Load relevant sub-programs
Cit_par
ISA_converted

ET          = xlsread(datasheet, 'C59:C65');
hp          = xlsread(datasheet, 'D59:D65').*unitsratio('meter', 'feet');
IAS         = convvel(xlsread(datasheet, 'E59:E65'), 'kts', 'm/s');
alpha       = convang(xlsread(datasheet, 'F59:F65'), 'deg', 'rad');
delta_e     = convang(xlsread(datasheet, 'G59:G65'), 'deg', 'rad');
delta_e_tr  = convang(xlsread(datasheet, 'H59:H65'), 'deg', 'rad');
F_e         = xlsread(datasheet, 'I59:I65');
FF_l        = convvel((convmass(xlsread(datasheet, 'J59:J65'), 'lbm', 'kg')), 'km/h', 'km/s');
FF_r        = convvel((convmass(xlsread(datasheet, 'K59:K65'), 'lbm', 'kg')), 'km/h', 'km/s');
F_used      = convmass(xlsread(datasheet, 'L59:L65'), 'lbm', 'kg');
TAT         = convtemp(xlsread(datasheet, 'M59:M65'), 'C', 'K');

%Start with computing Cm_d:
V_c_Cm_d  = mean(xlsread(datasheet, 'E75:E76'));
hp_Cm_d   = mean(xlsread(datasheet, 'D75:D76'))*unitsratio('meter', 'feet');
Delta_de_Cm_d  = convang((xlsread(datasheet, 'G75')-xlsread(datasheet, 'G76')), 'deg', 'rad');

T_Cm_d = mean(xlsread(datasheet, 'M75:M76'))/(1 + ((gamma-1)/2)*M^2);

p_Cm_d   = p0*((1+(lambda*(hp_Cm_d)/Temp0)^(-g0/(lambda*R)));
rho_Cm_d = rho0*((1+(lambda*hp_Cm_d/Temp0)))^(-((g/(lambda*R))+1));
C_N      = W/(0.5*rho_Cm_d*V_c_Cm_d^2*S);
%Delta_cg = cg[13]-cg[14]
Delta_cg = -0.1

Cm_d = -(1/Delta_de_Cm_d)*C_N*Delta_cg/c

%Now compute Cm_alpha:


