close all
clear all
clc

%Load in values in a struct 'c'
c = load('Cit_par.mat');
c.Temp0 = 19+273.15;


%Aircraft empty weight
BEM = 9165;      % [lbs]

%% Stationary measurements series 1 processing
datasheet = 'Post_Flight_Datasheet_Flight_2_DD_14_3_2018.xlsx'; %name
[~, hp, IAS_SI, alpha, FFl, FFr, F_used, TAT_K, weight] = ExcelReader_S1(datasheet, BEM);

%Convert velocity
[~, V_t, TAT, M, rho] = reduced_eq_airspeed(IAS_SI, hp, TAT_K, weight, 0);

%Find thrust
T = ThrustFile(c, hp, M, TAT, FFl, FFr);

%Construct graphs
Cl = 9.81*weight./(0.5.*rho.*c.S.*V_t.^(2));
Cd = (T(:,1)+T(:,2))./(0.5.*rho.*c.S.*V_t.^(2));

%Fit Cl and Cd
Cl_fit_f = fittype('b*(x-a)');
Cd_fit_f = fittype('a + b*x^2');

Cl_fit = fit(alpha, Cl, Cl_fit_f, 'StartPoint', [0, 0]);
Cd_fit = fit(Cl_fit(alpha), Cd, Cd_fit_f, 'StartPoint', [0,1/(pi*c.A*c.e)]);

%Determine coefficients from functions
coeff_Cd = coeffvalues(Cd_fit);
coeff_Cl = coeffvalues(Cl_fit);

Cd0 = coeff_Cd(1);
e = 1/(pi*c.A*coeff_Cd(2));
Cla = coeff_Cl(2);
alpha0 = coeff_Cl(1);

%Plot solutions
figure(1)
plot(rad2deg(alpha), Cl_fit(alpha), 'xk-')
hold on
plot(rad2deg(alpha), Cl, 'o')
ylabel('C_{L} [-]')
xlabel('\alpha [degree]')

figure(2)
plot(Cl_fit(alpha), Cd_fit(Cl_fit(alpha)), 'xk-')
hold on
plot(Cl_fit(alpha), Cd, 'o')
ylabel('C_{D} [-]')
xlabel('C_{L} [-]')

figure(3)
plot(alpha, Cd_fit(Cl_fit(alpha)), 'xk-')
hold on
plot(alpha, Cd, 'o')
ylabel('C_{D} [-]')
xlabel('\alpha')

figure(4)
plot(Cl_fit(alpha).^2, Cd_fit(Cl_fit(alpha)), 'xk-')
hold on
plot(Cl_fit(alpha).^2, Cd, 'o')
ylabel('C_{D}')
xlabel('C_{L}^{2}')

disp('')
disp('-------These plots are obtained in the following conditions-------')
disp('')
%% Aircraft configuration
disp('Aircraft has gears and flaps up')

%% Mach number range
M_range = string(num2str(min(M),'%.4f'))+' - '+string(num2str(max(M),'%.4f'));
disp('Mach number range is '+M_range)

%% Reynolds number range
% Function for dynamic viscosity
mu0 = 18.27e-6;     %[Pa*s]
T0 = 291.15;        %[K]
C = 120;            %[K]

mu= @(T) mu0*(T0+C)./(T+C).*(T/T0).^(3/2);

Re = rho.*V_t.*c.c./mu(TAT);
Re_range = string(num2str(min(Re),'%.3d'))+' - '+string(num2str(max(Re),'%.3d'));
disp('Reynolds number range is '+Re_range)

%% Display coefficients on screen
disp('')
disp('-------Next the coefficients derived from the flight test data are-------')
disp('')
disp('alpha_{0} = '+string(alpha0))
disp('C_{L_{alpha}} = '+string(Cla))
disp('e = '+string(e))
disp('Cd_{0} = '+string(Cd0))