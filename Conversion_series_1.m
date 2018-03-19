close all
clear all

%Load in values in a struct 'c'
c = load('Cit_par.mat');
c.Temp0 = 19+273.15;

%% Stationary measurements series 1 processing
hp_meas = [5064;5063;5070;5060;5065;5050];  %[ft]
IAS_meas = [248.3;220.3;193.3;160.3;143.5;113];       %[kts]
alpha_meas = [1.47;2.1;3.1;5.1;6.5;11];    %[deg]
FFl_meas = [770;674;553.8;453.8;429.8;424];       %[lbs/hr]
FFr_meas = [783.3;684.2;573.6;481.3;451.7;445];       %[lbs/hr]
Fused_meas = [395.83;426.5;457.2;480.8;500.8;539];     %[lbs]
TAT_meas = [9;7.2;5.2;3.8;2.65;1.8];      %[degree C]

%Initial parameters
Wi = 2000*0.45359237+ 92+82+130+78+73+75+89+74.5+92 +3655;              %[lbs]

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
Cd = (T(:,1)+T(:,2))./(0.5.*rho.*c.S.*V_t.^(2));

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