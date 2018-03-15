close all
clear all

%Load in values in a struct 'c'
c = load('Cit_par.mat');
c.Temp0 = 19+273.15;

%% Stationary measurements series 1 processing
hp_meas = [6040;6035;6030;6040;6030;6040];  %[ft]
IAS_meas = [251;222;192;165;130;115];       %[kts]
alpha_meas = [1.3;2.0;3.0;4.5;8.3;10.1];    %[deg]
FFl_meas = [786;637;515;458;414;373];       %[lbs/hr]
FFr_meas = [776;650;548;470;435;385];       %[lbs/hr]
Fused_meas = [410;450;487;540;574;604];     %[lbs]
TAT_meas = [10.5;8.2;6.1;4.8;3.5;3.0];      %[degree C]

%Initial parameters
Wi = 2800*0.45359237+ 92+95+76+61+59+66+77+84 +3655;              %[lbs]

%Convert all measurements to SI units
hp = hp_meas * 0.3048;
alpha = deg2rad(alpha_meas);
weight = Wi - Fused_meas*0.45359237;
IAS_SI = IAS_meas*0.514444444;
TAT_K = TAT_meas + 273.15;
[~, V_t, TAT, M, rho] = reduced_eq_airspeed(IAS_SI, hp_meas, TAT_K, weight, 0);
FFl = FFl_meas*0.3048/3600;
FFr = FFr_meas*0.3048/3600;
T = ThrustFile(c, hp, M, TAT, FFl, FFr);

%Construct graphs
Cl = weight./(0.5.*rho.*c.S.*V_t.^(2));
Cd = 2*T./(0.5.*rho.*c.S.*V_t.^(2));

Cl_fit_f = fittype('b*(x-a)');
Cd_fit_f = fittype('a + b*x^2');

Cl_fit = fit(alpha, Cl, Cl_fit_f, 'StartPoint', [0, 0]);
Cd_fit = fit(Cl_fit(alpha), Cd, Cd_fit_f, 'StartPoint', [0,0]);

%Plot solutions
figure(1)
plot(rad2deg(alpha), Cl_fit(alpha), 'xr-')
hold on
plot(rad2deg(alpha), Cl, 'ok-')

figure(2)
plot(Cl_fit(alpha), Cd_fit(Cl_fit(alpha)), 'xr-')
hold on
plot(Cl_fit(alpha), Cd, 'ok-')