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
Vc_s2          = convvel(xlsread(datasheet, 'E59:E65'), 'kts', 'm/s');
alpha_s2       = convang(xlsread(datasheet, 'F59:F65'), 'deg', 'rad');
de_s2          = convang(xlsread(datasheet, 'G59:G65'), 'deg', 'rad');
de_tr_s2       = convang(xlsread(datasheet, 'H59:H65'), 'deg', 'rad');
Fe_s2          = xlsread(datasheet, 'I59:I65');
FFl_s2         = convvel((convmass(xlsread(datasheet, 'J59:J65'), 'lbm', 'kg')), 'km/h', 'km/s');
FFr_s2         = convvel((convmass(xlsread(datasheet, 'K59:K65'), 'lbm', 'kg')), 'km/h', 'km/s');
F_used_s2      = convmass(xlsread(datasheet, 'L59:L65'), 'lbm', 'kg');
TAT_s2         = convtemp(xlsread(datasheet, 'M59:M65'), 'C', 'K');

%Start with computing Cm_d:
Vc_s3        = convvel(mean(xlsread(datasheet, 'E75:E76')), 'kts', 'm/s');
hp_s3        = mean(xlsread(datasheet, 'D75:D76')).*unitsratio('meter', 'feet');
TAT_s3       = convtemp(mean(xlsread(datasheet, 'M75:M76')), 'C', 'K');
Delta_de_s3  = convang((xlsread(datasheet, 'G76:G76')-xlsread(datasheet, 'G75:G75')), 'deg', 'rad');

[h_Cm_d, V_Cm_d, T_Cm_d, p_Cm_d, rho_Cm_d] = ISA_converted(hp_s3, ...
    Vc_s3, TAT_s3);

C_N      = (Wi-g*convmass(mean(xlsread(datasheet, 'L75:L76')), 'lbm', 'kg'))...
    /(0.5*rho_Cm_d*V_Cm_d^2*S);

%Delta_cg = cg[13]-cg[14]
Delta_cg = -0.0537;

Cm_d = -(1/Delta_de_s3)*C_N*Delta_cg/c.c;

%--------------------------------------------------------------------------
%Now start computing Cm_alpha, by first obtaining the reduced airspeed. 
%These calculations use measurement series 2.

W_s2 = -F_used_s2 + Wi;
<<<<<<< HEAD
[V_r_eq_s2, V_t_s2, Temp_s2, M_s2, rho_s2] = reduced_eq_airspeed(Vc_s2, hp_s2, TAT_s2, W_s2, Ws);
=======
[V_r_eq, V_t_s2, Temp_s2, M_s2, rho_s2] = reduced_eq_airspeed(IAS_s2, hp_s2, TAT_s2, W_s2, Ws);
>>>>>>> 1513f994b678d122aca6dedc06d43c1968494faa

[hp_s2, Vc_s2, Temp_s2, p_s2, rho_s2] = ISA_converted(hp_s2, Vc_s2, TAT_s2);

%Then the thrust coefficient and standard thrust coefficient at all 
%measurement times is required. Therefor create columns for the altitude, 
%Mach number, temperature difference from ISA conditions and fuel flows to 
%the left and right engine to use these as input for the thrust calculation
%program.

T_s2 = ThrustFile(c, hp_s2, M_s2, TAT_s2, FFl_s2, FFr_s2);

mf_s_s2 = ones(7,1)*mf_s;
Ts_s2 = ThrustFile(c, hp_s2, M_s2, TAT_s2, mf_s_s2, mf_s_s2);

T_tot_s2 = T_s2(:,1)+T_s2(:,2);
Ts_tot_s2 = Ts_s2(:,1)+Ts_s2(:,2);

T_c_s2 = T_tot_s2./(0.5.*rho_s2.*Vc_s2.^2.*S);
Ts_c_s2 = Ts_tot_s2./(0.5.*rho_s2.*V_r_eq_s2.^2.*S);

%Now find the reduced equivalent elevator deflection.
de_r_eq_s2 = de_s2 - (1/Cm_d).*Cm_Tc.*(Ts_c_s2-T_c_s2);
de_s2;

%Now plot it against alpha and against V_r_eq, and find the derrivative of
%the reduced equivalent elevator deflection w.r.t. alpha to obtain Cm_alpha
de_alpha_s2_matrix = [alpha_s2, de_r_eq_s2];
[~,idx] = sort(de_alpha_s2_matrix(:,1));            % sort just the first column
de_alpha_s2_sorted = de_alpha_s2_matrix(idx,:);     % sort the whole matrix using the sort indices

de_alpha_fit_f = fittype('a*x + b');
de_alpha_fit = fit(de_alpha_s2_sorted(:,1), de_alpha_s2_sorted(:,2),...
    de_alpha_fit_f, 'StartPoint', [0, 0]);

% plot(sort(V_r_eq_s2), sort(de_r_eq_s2))
plot(de_alpha_s2_sorted(:,1), de_alpha_s2_sorted(:,2))

dde_dalpha = de_alpha_fit.a;

Cm_alpha = -Cm_d*dde_dalpha
Cm_d