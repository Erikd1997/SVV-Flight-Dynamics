%This program handles the calculation of Cm_alpha and Cm_delta from our
%flight test data.

%First, import the data from the excel file and convert them into SI-units:
datasheet = 'Post_Flight_Datasheet_Flight_2_DD_6_3_2018_for_test.xlsx';

%Load relevant sub-programs
Cit_par
ISA_converted

ET_s1          = xlsread(datasheet, 'C59:C65');
hp_s1          = xlsread(datasheet, 'D59:D65').*unitsratio('meter', 'feet');
IAS_s1         = convvel(xlsread(datasheet, 'E59:E65'), 'kts', 'm/s');
alpha_s1       = convang(xlsread(datasheet, 'F59:F65'), 'deg', 'rad');
delta_e_s1     = convang(xlsread(datasheet, 'G59:G65'), 'deg', 'rad');
delta_e_tr_s1  = convang(xlsread(datasheet, 'H59:H65'), 'deg', 'rad');
Fe_s1          = xlsread(datasheet, 'I59:I65');
FFl_s1         = convvel((convmass(xlsread(datasheet, 'J59:J65'), 'lbm', 'kg')), 'km/h', 'km/s');
FFr_s1         = convvel((convmass(xlsread(datasheet, 'K59:K65'), 'lbm', 'kg')), 'km/h', 'km/s');
F_used_s1      = convmass(xlsread(datasheet, 'L59:L65'), 'lbm', 'kg');
TAT_s1         = convtemp(xlsread(datasheet, 'M59:M65'), 'C', 'K');

%Start with computing Cm_d:
IAS_kts_s3   = mean(xlsread(datasheet, 'E75:E76'));
hp_ft_s3     = mean(xlsread(datasheet, 'D75:D76'));
TAT_C_s3     = mean(xlsread(datasheet, 'M75:M76'));
Delta_de_s3  = convang((xlsread(datasheet, 'G75')-xlsread(datasheet, 'G76')), 'deg', 'rad');

[h_Cm_d, V_Cm_d, T_Cm_d, p_Cm_d, rho_Cm_d] = ISA_converted(hp_ft_s3, ...
    IAS_kts_s3, TAT_C_s3);

C_N      = W/(0.5*rho_Cm_d*V_c_Cm_d^2*S);

%Delta_cg = cg[13]-cg[14]
Delta_cg = -0.1;

Cm_d = -(1/Delta_de_s3)*C_N*Delta_cg/c

%Now compute Cm_alpha:


