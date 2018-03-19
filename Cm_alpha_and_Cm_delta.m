close all
clear all

%This program handles the calculation of Cm_alpha and Cm_delta from our
%flight test data.

%First, import the data from the excel file and convert them into SI-units:
datasheet = 'Post_Flight_Datasheet_Flight_2_DD_6_3_2018_for_test.xlsx';

%Load relevant sub-programs
Cit_par

%Save these values in a struct
save('Cit_par.mat')

%Load in values in a struct 'c'
c = load('Cit_par.mat');

ET_s2          = xlsread(datasheet, 'C59:C65');
hp_s2          = xlsread(datasheet, 'D59:D65').*unitsratio('meter', 'feet');
IAS_s2         = convvel(xlsread(datasheet, 'E59:E65'), 'kts', 'm/s');
alpha_s2       = convang(xlsread(datasheet, 'F59:F65'), 'deg', 'rad');
delta_e_s2     = convang(xlsread(datasheet, 'G59:G65'), 'deg', 'rad');
delta_e_tr_s2  = convang(xlsread(datasheet, 'H59:H65'), 'deg', 'rad');
Fe_s2          = xlsread(datasheet, 'I59:I65');
FFl_s2         = convvel((convmass(xlsread(datasheet, 'J59:J65'), 'lbm', 'kg')), 'km/h', 'km/s');
FFr_s2         = convvel((convmass(xlsread(datasheet, 'K59:K65'), 'lbm', 'kg')), 'km/h', 'km/s');
F_used_s2      = convmass(xlsread(datasheet, 'L59:L65'), 'lbm', 'kg');
TAT_s2         = convtemp(xlsread(datasheet, 'M59:M65'), 'C', 'K');

%Start with computing Cm_d:
IAS_kts_s3   = mean(xlsread(datasheet, 'E75:E76'));
hp_ft_s3     = mean(xlsread(datasheet, 'D75:D76'));
TAT_C_s3     = mean(xlsread(datasheet, 'M75:M76'));
Delta_de_s3  = convang((xlsread(datasheet, 'G76:G76')-xlsread(datasheet, 'G75:G75')), 'deg', 'rad');

[h_Cm_d, V_Cm_d, T_Cm_d, p_Cm_d, rho_Cm_d] = ISA_converted(hp_ft_s3, ...
    IAS_kts_s3, TAT_C_s3);

C_N      = W/(0.5*rho_Cm_d*V_Cm_d^2*S);

%Delta_cg = cg[13]-cg[14]
Delta_cg = -0.1;

Cm_d = -(1/Delta_de_s3)*C_N*Delta_cg/c.c;

%Now start computing Cm_alpha, by obtaining the reduced airspeed.

W_s2 = -F_used_s2 + Wi;
[V_r_eq, V_t_s2, Temp_s2, M_s2, rho_s2] = reduced_eq_airspeed(IAS_s2, hp_s2, TAT_s2, W_s2, Ws);

%for which the thrust at all measurement times is required. Therefor create
%columns for the altitude, Mach number, temperature difference from ISA 
%conditions and fuel flows to the left and right engine to use these as 
%input for the thrust calculation program.

T = ThrustFile(c, hp_s2, M_s2, TAT_s2, FFl_s2, FFr_s2);
