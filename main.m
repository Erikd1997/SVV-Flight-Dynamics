%% Opening regular Cit_par
clc
clear
%Call the file containing all variables
Cit_par

%Save these values in a struct
save('Cit_par.mat')

%Delete all saved variables from workspace
clear
clc
%% Opening updated Cit_par
%Call the file containing all variables
Cit_par_updated

%Save these values in a struct
save('Cit_par_updated.mat')

%Delete all saved variables from workspace
clear
clc

%Load in necessary values in a struct
c_init = load('Cit_par.mat');
c_update = load('Cit_par_updated.mat');
flight = load('FTISxprt-20180314_101817.mat');
datasheet = 'Post_Flight_Datasheet_Flight_2_DD_14_3_2018.xlsx';

%Which plots should appear?
plot_DutchRoll = false;
plot_Phugoid = false;
plot_ShortPeriod = false;
plot_AperiodicRoll = false;
plot_Spiral = false;

%Do you wanna plot all the necessary figures for the report to save them?
PlotAllForSaving = false;

%Do you wanna plot the result plots for results in numerical model?
ResultPlots = true;

%% Check each manoeuvre
%First save necessary data
da = flight.flightdata.delta_a.data;
de = flight.flightdata.delta_e.data;
dr = flight.flightdata.delta_r.data;
time = flight.flightdata.time.data;

Roll = (flight.flightdata.Ahrs1_Roll.data);
Pitch = (flight.flightdata.Ahrs1_Pitch.data);
RollRate = (flight.flightdata.Ahrs1_bRollRate.data);
PitchRate = (flight.flightdata.Ahrs1_bPitchRate.data);
YawRate = (flight.flightdata.Ahrs1_bYawRate.data);
AOA = (flight.flightdata.vane_AOA.data);
TAS = flight.flightdata.Dadc1_tas.data*0.514444444;
hp = flight.flightdata.Dadc1_alt.data;
FU_left = flight.flightdata.lh_engine_FU.data;
FU_right = flight.flightdata.rh_engine_FU.data;
Wi = (9165+2800)*0.45359237 + 82+92+67+81+80+75+83+76+95;

%Read times of each manoeuvre from datasheet
t_DutchRoll = (44*60+43.7)  *10;
t_DutchRoll_YD = 44*60      *10;
t_Phugoid = (49*60+10)      *10;
t_ShortPeriod = (47*60+50.5)*10;
t_AperiodicRoll = (46*60+57)*10;
t_Spiral = (54*60+58)       *10;

%Specify how long each manoeuvre takes
DR_length = 20              *10;
DR_YD_length = 70           *10;
Phugoid_length = 200        *10;
ShortPeriod_length = 10     *10;
AperiodicRoll_length = 15   *10;
Spiral_length = 80         *10;

%Dutch Roll
dr_DutchRoll = dr(t_DutchRoll:t_DutchRoll+DR_length);
da_DutchRoll = da(t_DutchRoll:t_DutchRoll+DR_length);
Roll_DutchRoll = Roll(t_DutchRoll:t_DutchRoll+DR_length);
RollRate_DutchRoll = RollRate(t_DutchRoll:t_DutchRoll+DR_length);
YawRate_DutchRoll = YawRate(t_DutchRoll:t_DutchRoll+DR_length);
time_DutchRoll = time(t_DutchRoll:t_DutchRoll+DR_length);
time_DutchRoll = time_DutchRoll - time_DutchRoll(1);

%Phugoid            %DONE
de_Phugoid = de(t_Phugoid:t_Phugoid+Phugoid_length);
AOA_Phugoid = AOA(t_Phugoid:t_Phugoid+Phugoid_length);
Pitch_Phugoid = Pitch(t_Phugoid:t_Phugoid+Phugoid_length);
PitchRate_Phugoid = PitchRate(t_Phugoid:t_Phugoid+Phugoid_length);
TAS_Phugoid = TAS(t_Phugoid:t_Phugoid+Phugoid_length);
time_Phugoid = time(t_Phugoid:t_Phugoid+Phugoid_length);
time_Phugoid = time_Phugoid - time_Phugoid(1);

%Short Period       %DONE
de_ShortPeriod = de(t_ShortPeriod:t_ShortPeriod+ShortPeriod_length);
AOA_ShortPeriod = AOA(t_ShortPeriod:t_ShortPeriod+ShortPeriod_length);
Pitch_ShortPeriod = Pitch(t_ShortPeriod:t_ShortPeriod+ShortPeriod_length);
PitchRate_ShortPeriod = PitchRate(t_ShortPeriod:t_ShortPeriod+ShortPeriod_length);
TAS_ShortPeriod = TAS(t_ShortPeriod:t_ShortPeriod+ShortPeriod_length);
time_ShortPeriod = time(t_ShortPeriod:t_ShortPeriod+ShortPeriod_length);
time_ShortPeriod = time_ShortPeriod - time_ShortPeriod(1);

%Aperiodic Roll     %DONE
dr_AperiodicRoll = dr(t_AperiodicRoll:t_AperiodicRoll+AperiodicRoll_length);
da_AperiodicRoll = da(t_AperiodicRoll:t_AperiodicRoll+AperiodicRoll_length);
Roll_AperiodicRoll = Roll(t_AperiodicRoll:t_AperiodicRoll+AperiodicRoll_length);
RollRate_AperiodicRoll = RollRate(t_AperiodicRoll:t_AperiodicRoll+AperiodicRoll_length);
YawRate_AperiodicRoll = YawRate(t_AperiodicRoll:t_AperiodicRoll+AperiodicRoll_length);
time_AperiodicRoll = time(t_AperiodicRoll:t_AperiodicRoll+AperiodicRoll_length);
time_AperiodicRoll = time_AperiodicRoll - time_AperiodicRoll(1);

%Spiral             %DONE
dr_Spiral = dr(t_Spiral:t_Spiral+Spiral_length);
da_Spiral = da(t_Spiral:t_Spiral+Spiral_length);
Roll_Spiral = Roll(t_Spiral:t_Spiral+Spiral_length);
RollRate_Spiral = RollRate(t_Spiral:t_Spiral+Spiral_length);
YawRate_Spiral = YawRate(t_Spiral:t_Spiral+Spiral_length);
time_Spiral = time(t_Spiral:t_Spiral+Spiral_length);
time_Spiral = time_Spiral - time_Spiral(1);

%% Solve simulations of each manouevre from initial cit_par file
%------Dutch Roll------
c_init.hp0    = hp(t_DutchRoll)*0.3048;
c_init.V0     = TAS(t_DutchRoll);
c_init.alpha0 = deg2rad(AOA(t_DutchRoll));
c_init.th0    = deg2rad(Pitch(t_DutchRoll));
c_init.rho    = c_init.rho0*((1+(c_init.lambda*c_init.hp0/c_init.Temp0)))^(-((c_init.g/(c_init.lambda*c_init.R))+1));
c_init.m      = Wi - FU_left(t_DutchRoll) - FU_right(t_DutchRoll);
c_init.W      = c_init.m*c_init.g;
c_init.muc    = c_init.m/(c_init.rho*c_init.S*c_init.c);
c_init.mub    = c_init.m/(c_init.rho*c_init.S*c_init.b);
c_init.CL     = 2*c_init.W/(c_init.rho*c_init.V0^2*c_init.S);
c_init.CD     = c_init.CD0 + (c_init.alpha0*c_init.CLa)^2/(pi*c_init.A*c_init.e);
c_init.CX0    = c_init.W*sin(c_init.th0)/(0.5*c_init.rho*c_init.V0^2*c_init.S);
c_init.CZ0    = -c_init.W*cos(c_init.th0)/(0.5*c_init.rho*c_init.V0^2*c_init.S);

%Create state-space-system
[~, sysA_DutchRoll_init] = state_space_system(c_init);

%Save and non-dimensionalize eigenvalues
eig_DutchRoll_init = eig(sysA_DutchRoll_init)*c_init.b/c_init.V0;

%Initial condition
x0 = deg2rad([0, Roll_DutchRoll(1), RollRate_DutchRoll(1), YawRate_DutchRoll(1)]);
%Simulate using state space system
y_DutchRoll_init = lsim(sysA_DutchRoll_init, deg2rad([da_DutchRoll dr_DutchRoll]), time_DutchRoll, x0);

% Transform
y_DutchRoll_init(:,1) = rad2deg(y_DutchRoll_init(:,1));
y_DutchRoll_init(:,2) = rad2deg(y_DutchRoll_init(:,2));
y_DutchRoll_init(:,3) = rad2deg(y_DutchRoll_init(:,3));
y_DutchRoll_init(:,4) = rad2deg(y_DutchRoll_init(:,4));

%------Phugoid------
c_init.hp0 = hp(t_Phugoid)*0.3048;
c_init.V0 = TAS(t_Phugoid);
c_init.alpha0 = deg2rad(AOA(t_Phugoid));
c_init.th0 = deg2rad(Pitch(t_Phugoid));
c_init.rho = c_init.rho0*((1+(c_init.lambda*c_init.hp0/c_init.Temp0)))^(-((c_init.g/(c_init.lambda*c_init.R))+1));
c_init.m = Wi - FU_left(t_Phugoid) - FU_right(t_Phugoid);
c_init.W = c_init.m*c_init.g;
c_init.muc    = c_init.m/(c_init.rho*c_init.S*c_init.c);
c_init.mub    = c_init.m/(c_init.rho*c_init.S*c_init.b);
c_init.CL = 2*c_init.W/(c_init.rho*c_init.V0^2*c_init.S);
c_init.CD = c_init.CD0 + (c_init.alpha0*c_init.CLa)^2/(pi*c_init.A*c_init.e);
c_init.CX0    = c_init.W*sin(c_init.th0)/(0.5*c_init.rho*c_init.V0^2*c_init.S);
c_init.CZ0    = -c_init.W*cos(c_init.th0)/(0.5*c_init.rho*c_init.V0^2*c_init.S);

%Create state-space-system
[sysS_Phugoid_init, ~] = state_space_system(c_init);

%Save and non-dimensionalize eigenvalues
eig_Phugoid_init = eig(sysS_Phugoid_init)*c_init.c/c_init.V0;

%Simulate using state space system
y_Phugoid_init = lsim(sysS_Phugoid_init, deg2rad(de_Phugoid), time_Phugoid);
% Transform
y_Phugoid_init(:,1) = y_Phugoid_init(:,1)+c_init.V0;
y_Phugoid_init(:,2) = rad2deg(y_Phugoid_init(:,2) + c_init.alpha0);
y_Phugoid_init(:,3) = rad2deg(y_Phugoid_init(:,3) + c_init.th0);
y_Phugoid_init(:,4) = rad2deg(y_Phugoid_init(:,4));

%------Short Period------
c_init.hp0 = hp(t_ShortPeriod)*0.3048;
c_init.V0 = TAS(t_ShortPeriod);
c_init.alpha0 = deg2rad(AOA(t_ShortPeriod));
c_init.th0 = deg2rad(Pitch(t_ShortPeriod));
c_init.rho = c_init.rho0*((1+(c_init.lambda*c_init.hp0/c_init.Temp0)))^(-((c_init.g/(c_init.lambda*c_init.R))+1));
c_init.m = Wi - FU_left(t_ShortPeriod) - FU_right(t_ShortPeriod);
c_init.W = c_init.m*c_init.g;
c_init.muc    = c_init.m/(c_init.rho*c_init.S*c_init.c);
c_init.mub    = c_init.m/(c_init.rho*c_init.S*c_init.b);
c_init.CL = 2*c_init.W/(c_init.rho*c_init.V0^2*c_init.S);
c_init.CD = c_init.CD0 + (c_init.alpha0*c_init.CLa)^2/(pi*c_init.A*c_init.e);
c_init.CX0    = c_init.W*sin(c_init.th0)/(0.5*c_init.rho*c_init.V0^2*c_init.S);
c_init.CZ0    = -c_init.W*cos(c_init.th0)/(0.5*c_init.rho*c_init.V0^2*c_init.S);

%Create state-space-system
[sysS_ShortPeriod_init, ~] = state_space_system(c_init);

%Save and non-dimensionalize eigenvalues
eig_ShortPeriod_init = eig(sysS_ShortPeriod_init)*c_init.c/c_init.V0;

%Simulate using state space system
y_ShortPeriod_init = lsim(sysS_ShortPeriod_init, deg2rad(de_ShortPeriod), time_ShortPeriod);

% Transform
y_ShortPeriod_init(:,1) = y_ShortPeriod_init(:,1)+c_init.V0;
y_ShortPeriod_init(:,2) = rad2deg(y_ShortPeriod_init(:,2) + c_init.alpha0);
y_ShortPeriod_init(:,3) = rad2deg(y_ShortPeriod_init(:,3) + c_init.th0);
y_ShortPeriod_init(:,4) = rad2deg(y_ShortPeriod_init(:,4));

%------Aperiodic Roll------
c_init.hp0 = hp(t_AperiodicRoll)*0.3048;
c_init.V0 = TAS(t_AperiodicRoll);
c_init.alpha0 = deg2rad(AOA(t_AperiodicRoll));
c_init.th0 = deg2rad(Pitch(t_AperiodicRoll));
c_init.rho = c_init.rho0*((1+(c_init.lambda*c_init.hp0/c_init.Temp0)))^(-((c_init.g/(c_init.lambda*c_init.R))+1));
c_init.m = Wi - FU_left(t_AperiodicRoll) - FU_right(t_AperiodicRoll);
c_init.W = c_init.m*c_init.g;
c_init.muc    = c_init.m/(c_init.rho*c_init.S*c_init.c);
c_init.mub    = c_init.m/(c_init.rho*c_init.S*c_init.b);
c_init.CL = 2*c_init.W/(c_init.rho*c_init.V0^2*c_init.S);
c_init.CD = c_init.CD0 + (c_init.alpha0*c_init.CLa)^2/(pi*c_init.A*c_init.e);
c_init.CX0    = c_init.W*sin(c_init.th0)/(0.5*c_init.rho*c_init.V0^2*c_init.S);
c_init.CZ0    = -c_init.W*cos(c_init.th0)/(0.5*c_init.rho*c_init.V0^2*c_init.S);

%Create state-space-system
[~, sysA_AperiodicRoll_init] = state_space_system(c_init);

%Initial condition
x0 = [0, deg2rad(Roll_AperiodicRoll(1)), deg2rad(RollRate_AperiodicRoll(1)), deg2rad(YawRate_AperiodicRoll(1))];

%Save and non-dimensionalize eigenvalues
eig_AperiodicRoll_init = eig(sysA_AperiodicRoll_init)*c_init.b/c_init.V0;

%Simulate using state space system
y_AperiodicRoll_init = lsim(sysA_AperiodicRoll_init, deg2rad([da_AperiodicRoll-da_AperiodicRoll(1) dr_AperiodicRoll-dr_AperiodicRoll(1)]), time_AperiodicRoll, x0);

% Transform
y_AperiodicRoll_init(:,1) = rad2deg(y_AperiodicRoll_init(:,1));
y_AperiodicRoll_init(:,2) = rad2deg(y_AperiodicRoll_init(:,2));
y_AperiodicRoll_init(:,3) = rad2deg(y_AperiodicRoll_init(:,3));
y_AperiodicRoll_init(:,4) = rad2deg(y_AperiodicRoll_init(:,4));

%------Spiral------
c_init.hp0 = hp(t_Spiral)*0.3048;
c_init.V0 = TAS(t_Spiral);
c_init.alpha0 = deg2rad(AOA(t_Spiral));
c_init.th0 = deg2rad(Pitch(t_Spiral));
c_init.rho = c_init.rho0*((1+(c_init.lambda*c_init.hp0/c_init.Temp0)))^(-((c_init.g/(c_init.lambda*c_init.R))+1));
c_init.m = Wi - FU_left(t_Spiral) - FU_right(t_Spiral);
c_init.W = c_init.m*c_init.g;
c_init.muc    = c_init.m/(c_init.rho*c_init.S*c_init.c);
c_init.mub    = c_init.m/(c_init.rho*c_init.S*c_init.b);
c_init.CL = 2*c_init.W/(c_init.rho*c_init.V0^2*c_init.S);
c_init.CD = c_init.CD0 + (c_init.alpha0*c_init.CLa)^2/(pi*c_init.A*c_init.e);
c_init.CX0    = c_init.W*sin(c_init.th0)/(0.5*c_init.rho*c_init.V0^2*c_init.S);
c_init.CZ0    = -c_init.W*cos(c_init.th0)/(0.5*c_init.rho*c_init.V0^2*c_init.S);

%Create state-space-system
[~, sysA_Spiral_init] = state_space_system(c_init);

%Initial state
x0_Spiral = [0, deg2rad(Roll(t_Spiral)), deg2rad(RollRate(t_Spiral)), deg2rad(YawRate(t_Spiral))];

%Save and non-dimensionalize eigenvalues
eig_Spiral_init = eig(sysA_Spiral_init)*c_init.b/c_init.V0;

%Simulate using state space system
y_Spiral_init = lsim(sysA_Spiral_init, deg2rad([da_Spiral dr_Spiral]), time_Spiral, x0_Spiral);

% Transform
y_Spiral_init(:,1) = rad2deg(y_Spiral_init(:,1));
y_Spiral_init(:,2) = rad2deg(y_Spiral_init(:,2));
y_Spiral_init(:,3) = rad2deg(y_Spiral_init(:,3));
y_Spiral_init(:,4) = rad2deg(y_Spiral_init(:,4));

%% Solve simulations of each manouevre from updated cit_par file
%------Dutch Roll------
c_update.hp0    = hp(t_DutchRoll)*0.3048;
c_update.V0     = TAS(t_DutchRoll);
c_update.alpha0 = deg2rad(AOA(t_DutchRoll));
c_update.th0    = deg2rad(Pitch(t_DutchRoll));
c_update.rho    = c_update.rho0*((1+(c_update.lambda*c_update.hp0/c_update.Temp0)))^(-((c_update.g/(c_update.lambda*c_update.R))+1));
c_update.m      = Wi - FU_left(t_DutchRoll) - FU_right(t_DutchRoll);
c_update.W      = c_update.m*c_update.g;
c_update.muc    = c_update.m/(c_update.rho*c_update.S*c_update.c);
c_update.mub    = c_update.m/(c_update.rho*c_update.S*c_update.b);
c_update.CL     = 2*c_update.W/(c_update.rho*c_update.V0^2*c_update.S);
c_update.CD     = c_update.CD0 + (c_update.alpha0*c_update.CLa)^2/(pi*c_update.A*c_update.e);
c_update.CX0    = c_update.W*sin(c_update.th0)/(0.5*c_update.rho*c_update.V0^2*c_update.S);
c_update.CZ0    = -c_update.W*cos(c_update.th0)/(0.5*c_update.rho*c_update.V0^2*c_update.S);

%Create state-space-system
[~, sysA_DutchRoll_update] = state_space_system(c_update);

%Save and non-dimensionalize eigenvalues
eig_DutchRoll_update = eig(sysA_DutchRoll_update)*c_update.b/c_update.V0;

%Initial condition
x0 = deg2rad([0, Roll_DutchRoll(1), RollRate_DutchRoll(1), YawRate_DutchRoll(1)]);
%Simulate using state space system
y_DutchRoll_update = lsim(sysA_DutchRoll_update, deg2rad([da_DutchRoll dr_DutchRoll]), time_DutchRoll, x0);

% Transform
y_DutchRoll_update(:,1) = rad2deg(y_DutchRoll_update(:,1));
y_DutchRoll_update(:,2) = rad2deg(y_DutchRoll_update(:,2));
y_DutchRoll_update(:,3) = rad2deg(y_DutchRoll_update(:,3));
y_DutchRoll_update(:,4) = rad2deg(y_DutchRoll_update(:,4));

%------Phugoid------
c_update.hp0 = hp(t_Phugoid)*0.3048;
c_update.V0 = TAS(t_Phugoid);
c_update.alpha0 = deg2rad(AOA(t_Phugoid));
c_update.th0 = deg2rad(Pitch(t_Phugoid));
c_update.rho = c_update.rho0*((1+(c_update.lambda*c_update.hp0/c_update.Temp0)))^(-((c_update.g/(c_update.lambda*c_update.R))+1));
c_update.m = Wi - FU_left(t_Phugoid) - FU_right(t_Phugoid);
c_update.W = c_update.m*c_update.g;
c_update.muc    = c_update.m/(c_update.rho*c_update.S*c_update.c);
c_update.mub    = c_update.m/(c_update.rho*c_update.S*c_update.b);
c_update.CL = 2*c_update.W/(c_update.rho*c_update.V0^2*c_update.S);
c_update.CD = c_update.CD0 + (c_update.alpha0*c_update.CLa)^2/(pi*c_update.A*c_update.e);
c_update.CX0    = c_update.W*sin(c_update.th0)/(0.5*c_update.rho*c_update.V0^2*c_update.S);
c_update.CZ0    = -c_update.W*cos(c_update.th0)/(0.5*c_update.rho*c_update.V0^2*c_update.S);

%Create state-space-system
[sysS_Phugoid_update, ~] = state_space_system(c_update);

%Save and non-dimensionalize eigenvalues
eig_Phugoid_update = eig(sysS_Phugoid_update)*c_update.c/c_update.V0;

%Simulate using state space system
y_Phugoid_update = lsim(sysS_Phugoid_update, deg2rad(de_Phugoid), time_Phugoid);
% Transform
y_Phugoid_update(:,1) = y_Phugoid_update(:,1)+c_update.V0;
y_Phugoid_update(:,2) = rad2deg(y_Phugoid_update(:,2) + c_update.alpha0);
y_Phugoid_update(:,3) = rad2deg(y_Phugoid_update(:,3) + c_update.th0);
y_Phugoid_update(:,4) = rad2deg(y_Phugoid_update(:,4));

%------Short Period------
c_update.hp0 = hp(t_ShortPeriod)*0.3048;
c_update.V0 = TAS(t_ShortPeriod);
c_update.alpha0 = deg2rad(AOA(t_ShortPeriod));
c_update.th0 = deg2rad(Pitch(t_ShortPeriod));
c_update.rho = c_update.rho0*((1+(c_update.lambda*c_update.hp0/c_update.Temp0)))^(-((c_update.g/(c_update.lambda*c_update.R))+1));
c_update.m = Wi - FU_left(t_ShortPeriod) - FU_right(t_ShortPeriod);
c_update.W = c_update.m*c_update.g;
c_update.muc    = c_update.m/(c_update.rho*c_update.S*c_update.c);
c_update.mub    = c_update.m/(c_update.rho*c_update.S*c_update.b);
c_update.CL = 2*c_update.W/(c_update.rho*c_update.V0^2*c_update.S);
c_update.CD = c_update.CD0 + (c_update.alpha0*c_update.CLa)^2/(pi*c_update.A*c_update.e);
c_update.CX0    = c_update.W*sin(c_update.th0)/(0.5*c_update.rho*c_update.V0^2*c_update.S);
c_update.CZ0    = -c_update.W*cos(c_update.th0)/(0.5*c_update.rho*c_update.V0^2*c_update.S);

%Create state-space-system
[sysS_ShortPeriod_update, ~] = state_space_system(c_update);

%Save and non-dimensionalize eigenvalues
eig_ShortPeriod_update = eig(sysS_ShortPeriod_update)*c_update.c/c_update.V0;

%Simulate using state space system
y_ShortPeriod_update = lsim(sysS_ShortPeriod_update, deg2rad(de_ShortPeriod), time_ShortPeriod);

% Transform
y_ShortPeriod_update(:,1) = y_ShortPeriod_update(:,1)+c_update.V0;
y_ShortPeriod_update(:,2) = rad2deg(y_ShortPeriod_update(:,2) + c_update.alpha0);
y_ShortPeriod_update(:,3) = rad2deg(y_ShortPeriod_update(:,3) + c_update.th0);
y_ShortPeriod_update(:,4) = rad2deg(y_ShortPeriod_update(:,4));

%------Aperiodic Roll------
c_update.hp0 = hp(t_AperiodicRoll)*0.3048;
c_update.V0 = TAS(t_AperiodicRoll);
c_update.alpha0 = deg2rad(AOA(t_AperiodicRoll));
c_update.th0 = deg2rad(Pitch(t_AperiodicRoll));
c_update.rho = c_update.rho0*((1+(c_update.lambda*c_update.hp0/c_update.Temp0)))^(-((c_update.g/(c_update.lambda*c_update.R))+1));
c_update.m = Wi - FU_left(t_AperiodicRoll) - FU_right(t_AperiodicRoll);
c_update.W = c_update.m*c_update.g;
c_update.muc    = c_update.m/(c_update.rho*c_update.S*c_update.c);
c_update.mub    = c_update.m/(c_update.rho*c_update.S*c_update.b);
c_update.CL = 2*c_update.W/(c_update.rho*c_update.V0^2*c_update.S);
c_update.CD = c_update.CD0 + (c_update.alpha0*c_update.CLa)^2/(pi*c_update.A*c_update.e);
c_update.CX0    = c_update.W*sin(c_update.th0)/(0.5*c_update.rho*c_update.V0^2*c_update.S);
c_update.CZ0    = -c_update.W*cos(c_update.th0)/(0.5*c_update.rho*c_update.V0^2*c_update.S);

%Create state-space-system
[~, sysA_AperiodicRoll_update] = state_space_system(c_update);

%Initial condition
x0 = [0, deg2rad(Roll_AperiodicRoll(1)), deg2rad(RollRate_AperiodicRoll(1)), deg2rad(YawRate_AperiodicRoll(1))];

%Save and non-dimensionalize eigenvalues
eig_AperiodicRoll_update = eig(sysA_AperiodicRoll_update)*c_update.b/c_update.V0;

%Simulate using state space system
y_AperiodicRoll_update = lsim(sysA_AperiodicRoll_update, deg2rad([da_AperiodicRoll-da_AperiodicRoll(1) dr_AperiodicRoll-dr_AperiodicRoll(1)]), time_AperiodicRoll, x0);

% Transform
y_AperiodicRoll_update(:,1) = rad2deg(y_AperiodicRoll_update(:,1));
y_AperiodicRoll_update(:,2) = rad2deg(y_AperiodicRoll_update(:,2));
y_AperiodicRoll_update(:,3) = rad2deg(y_AperiodicRoll_update(:,3));
y_AperiodicRoll_update(:,4) = rad2deg(y_AperiodicRoll_update(:,4));

%------Spiral------
c_update.hp0 = hp(t_Spiral)*0.3048;
c_update.V0 = TAS(t_Spiral);
c_update.alpha0 = deg2rad(AOA(t_Spiral));
c_update.th0 = deg2rad(Pitch(t_Spiral));
c_update.rho = c_update.rho0*((1+(c_update.lambda*c_update.hp0/c_update.Temp0)))^(-((c_update.g/(c_update.lambda*c_update.R))+1));
c_update.m = Wi - FU_left(t_Spiral) - FU_right(t_Spiral);
c_update.W = c_update.m*c_update.g;
c_update.muc    = c_update.m/(c_update.rho*c_update.S*c_update.c);
c_update.mub    = c_update.m/(c_update.rho*c_update.S*c_update.b);
c_update.CL = 2*c_update.W/(c_update.rho*c_update.V0^2*c_update.S);
c_update.CD = c_update.CD0 + (c_update.alpha0*c_update.CLa)^2/(pi*c_update.A*c_update.e);
c_update.CX0    = c_update.W*sin(c_update.th0)/(0.5*c_update.rho*c_update.V0^2*c_update.S);
c_update.CZ0    = -c_update.W*cos(c_update.th0)/(0.5*c_update.rho*c_update.V0^2*c_update.S);

%Create state-space-system
[~, sysA_Spiral_update] = state_space_system(c_update);

%Initial state
x0_Spiral = [0, deg2rad(Roll(t_Spiral)), deg2rad(RollRate(t_Spiral)), deg2rad(YawRate(t_Spiral))];

%Save and non-dimensionalize eigenvalues
eig_Spiral_update = eig(sysA_Spiral_update)*c_update.b/c_update.V0;

%Simulate using state space system
y_Spiral_update = lsim(sysA_Spiral_update, deg2rad([da_Spiral dr_Spiral]), time_Spiral, x0_Spiral);

% Transform
y_Spiral_update(:,1) = rad2deg(y_Spiral_update(:,1));
y_Spiral_update(:,2) = rad2deg(y_Spiral_update(:,2));
y_Spiral_update(:,3) = rad2deg(y_Spiral_update(:,3));
y_Spiral_update(:,4) = rad2deg(y_Spiral_update(:,4));

%% Plot Measurements
%Dutch Roll
if plot_DutchRoll
    figure(1)
    subplot(4,1,1)
    plot(time_DutchRoll, [RollRate_DutchRoll y_DutchRoll_init(:,3) y_DutchRoll_update(:,3)])
    legend({'Measurement', 'Initial values', 'Updated values'})
    ylabel('p (deg/s)')
    xlabel('t (s)')
    title('Roll Rate - Dutch Roll')
    subplot(4,1,2)
    plot(time_DutchRoll, [YawRate_DutchRoll y_DutchRoll_init(:,4) y_DutchRoll_update(:,4)])
    legend({'Measurement', 'Initial values', 'Updated values'})
    ylabel('r (deg/s)')
    xlabel('t (s)')
    title('Yaw Rate - Dutch Roll')
    subplot(4,1,3)
    plot([RollRate_DutchRoll y_DutchRoll_init(:,3) y_DutchRoll_update(:,3)], [YawRate_DutchRoll y_DutchRoll_init(:,4) y_DutchRoll_update(:,4)])
    legend({'Measurement', 'Initial values', 'Updated values'})
    ylabel('r (deg/s)')
    xlabel('p (deg/s)')
    title('Roll rate vs Yaw rate - Dutch Roll')
    subplot(4,1,4)
    plot(time_DutchRoll, [da_DutchRoll dr_DutchRoll])
    xlabel('t (s)')
    ylabel('deflection (deg)')
    xlabel('t (s)')
    legend({'Aileron','Rudder'})
    title('Inputs - Dutch Roll')
end

%Phugoid
if plot_Phugoid
    figure(2)
    subplot(5,1,1)
    plot(time_Phugoid, [TAS_Phugoid y_Phugoid_init(:,1) y_Phugoid_update(:,1)])
    legend({'Measurement', 'Initial values', 'Updated values'})
    ylabel('V (m/s)')
    xlabel('t (s)')
    title('True Airspeed - Phugoid')
    subplot(5,1,2)
    plot(time_Phugoid, [AOA_Phugoid y_Phugoid_init(:,2) y_Phugoid_update(:,2)])
    legend({'Measurement', 'Initial values', 'Updated values'})
    ylabel('\alpha (deg)')
    xlabel('t (s)')
    title('Angle of Attack - Phugoid')
    subplot(5,1,3)
    plot(time_Phugoid, [Pitch_Phugoid y_Phugoid_init(:,3) y_Phugoid_update(:,3)])
    legend({'Measurement', 'Initial values', 'Updated values'})
    ylabel('\theta (deg)')
    xlabel('t (s)')
    title('Pitch angle - Phugoid')
    subplot(5,1,4)
    plot(time_Phugoid, [PitchRate_Phugoid y_Phugoid_init(:,4) y_Phugoid_update(:,4)])
    legend({'Measurement', 'Initial values', 'Updated values'})
    ylabel('q (deg/s)')
    xlabel('t (s)')
    title('Pitch Rate - Phugoid')
    subplot(5,1,5)
    plot(time_Phugoid, de_Phugoid)
    ylabel('de (deg)')
    xlabel('t (s)')
    title('Elevator Deflection - Phugoid')
end

%Short Period
if plot_ShortPeriod
    figure(3)
    subplot(5,1,1)
    plot(time_ShortPeriod, [TAS_ShortPeriod y_ShortPeriod_init(:,1) y_ShortPeriod_update(:,1)])
    legend({'Measurement', 'Initial values', 'Updated values'})
    ylabel('V (m/s)')
    xlabel('t (s)')
    title('True Airspeed - Short Period')
    subplot(5,1,2)
    plot(time_ShortPeriod, [AOA_ShortPeriod y_ShortPeriod_init(:,2) y_ShortPeriod_update(:,2)])
    legend({'Measurement', 'Initial values', 'Updated values'})
    ylabel('\alpha (deg)')
    xlabel('t (s)')
    title('Angle of Attack - Short Period')
    subplot(5,1,3)
    plot(time_ShortPeriod, [Pitch_ShortPeriod y_ShortPeriod_init(:,3) y_ShortPeriod_update(:,3)])
    legend({'Measurement', 'Initial values', 'Updated values'})
    ylabel('\theta (deg)')
    xlabel('t (s)')
    title('Pitch Angle - Short Period')
    subplot(5,1,4)
    plot(time_ShortPeriod, [PitchRate_ShortPeriod y_ShortPeriod_init(:,4) y_ShortPeriod_update(:,4)])
    legend({'Measurement', 'Initial values', 'Updated values'})
    ylabel('q (deg/s)')
    xlabel('t (s)')
    title('Pitch Rate - Short Period')
    subplot(5,1,5)
    plot(time_ShortPeriod, de_ShortPeriod)
    ylabel('de (deg)')
    xlabel('t (s)')
    title('Elevator Deflection - Short Period')
end

%Aperiodic Roll
if plot_AperiodicRoll
    figure(4)
    subplot(3,1,1)
    plot(time_AperiodicRoll, [Roll_AperiodicRoll y_AperiodicRoll_init(:,2) y_AperiodicRoll_update(:,2)])
    legend({'Measurement', 'Initial values', 'Updated values'})
    ylabel('\phi (deg/s)')
    xlabel('t (s)')
    title('Roll Angle - Aperiodic Roll')
    subplot(3,1,2)
    plot(time_AperiodicRoll, [RollRate_AperiodicRoll y_AperiodicRoll_init(:,3) y_AperiodicRoll_update(:,3)])
    legend({'Measurement', 'Initial values', 'Updated values'})
    ylabel('p (deg/s)')
    xlabel('t (s)')
    title('Roll Rate - Aperiodic Roll')
    subplot(3,1,3)
    plot(time_AperiodicRoll, [da_AperiodicRoll dr_AperiodicRoll])
    ylabel('deflection (deg)')
    xlabel('t (s)')
    legend({'Aileron','Rudder'})
    title('Inputs - Aperiodic Roll')
end

%Spiral
if plot_Spiral
    figure(5)
    subplot(4,1,1)
    plot(time_Spiral, [Roll_Spiral y_Spiral_init(:,2) y_Spiral_update(:,2)])
    legend({'Measurement', 'Initial values', 'Updated values'})
    ylabel('\phi (deg)')
    xlabel('t (s)')
    title('Roll Angle - Spiral')
    subplot(4,1,2)
    plot(time_Spiral, [YawRate_Spiral y_Spiral_init(:,4) y_Spiral_update(:,4)])
    legend({'Measurement', 'Initial values', 'Updated values'})
    ylabel('r (deg/s)')
    xlabel('t (s)')
    title('Yaw Rate - Spiral')
    subplot(4,1,3)
    plot(time_Spiral, [RollRate_Spiral y_Spiral_init(:,3) y_Spiral_update(:,3)])
    legend({'Measurement', 'Initial values', 'Updated values'})
    ylabel('p (deg/s)')
    xlabel('t (s)')
    title('Roll Rate - Spiral')
    subplot(4,1,4)
    plot(time_Spiral, [da_Spiral dr_Spiral])
    ylabel('deflection (deg)')
    xlabel('t (s)')
    legend({'Aileron','Rudder'})
    title('Inputs - Spiral')
end

if PlotAllForSaving
    %Short Period
    figure
    plot(time_ShortPeriod, [TAS_ShortPeriod y_ShortPeriod_init(:,1) y_ShortPeriod_update(:,1)])
    legend({'Measurement', 'Initial values', 'Updated values'})
    ylabel('V [m/s]')
    xlabel('t [s]')
%     title('True Airspeed - Short Period')
    
    figure
    plot(time_ShortPeriod, [AOA_ShortPeriod y_ShortPeriod_init(:,2) y_ShortPeriod_update(:,2)])
    legend({'Measurement', 'Initial values', 'Updated values'})
    ylabel('\alpha [deg]')
    xlabel('t [s]')
%     title('Angle of Attack - Short Period')
    
    figure
    plot(time_ShortPeriod, [Pitch_ShortPeriod y_ShortPeriod_init(:,3) y_ShortPeriod_update(:,3)])
    legend({'Measurement', 'Initial values', 'Updated values'})
    ylabel('\theta [deg]')
    xlabel('t [s]')
%     title('Pitch Angle - Short Period')
    
    figure
    plot(time_ShortPeriod, [PitchRate_ShortPeriod y_ShortPeriod_init(:,4) y_ShortPeriod_update(:,4)])
    legend({'Measurement', 'Initial values', 'Updated values'})
    ylabel('q [deg/s]')
    xlabel('t [s]')
%     title('Pitch Rate - Short Period')
    
    figure
    plot(time_ShortPeriod, de_ShortPeriod)
    ylabel('de [deg]')
    xlabel('t [s]')
%     title('Elevator Deflection - Short Period')
    
    %Phugoid    
    figure
    plot(time_Phugoid, [TAS_Phugoid y_Phugoid_init(:,1) y_Phugoid_update(:,1)])
    legend({'Measurement', 'Initial values', 'Updated values'})
    ylabel('V [m/s]')
    xlabel('t [s]')
%     title('True Airspeed - Phugoid')
    
    figure
    plot(time_Phugoid, [AOA_Phugoid y_Phugoid_init(:,2) y_Phugoid_update(:,2)])
    legend({'Measurement', 'Initial values', 'Updated values'})
    ylabel('\alpha [deg]')
    xlabel('t [s]')
%     title('Angle of Attack - Phugoid')
    
    figure
    plot(time_Phugoid, [Pitch_Phugoid y_Phugoid_init(:,3) y_Phugoid_update(:,3)])
    legend({'Measurement', 'Initial values', 'Updated values'})
    ylabel('\theta [deg]')
    xlabel('t [s]')
%     title('Pitch angle - Phugoid')
    
    figure
    plot(time_Phugoid, [PitchRate_Phugoid y_Phugoid_init(:,4) y_Phugoid_update(:,4)])
    legend({'Measurement', 'Initial values', 'Updated values'})
    ylabel('q [deg/s]')
    xlabel('t [s]')
%     title('Pitch Rate - Phugoid')
    
    figure
    plot(time_Phugoid, de_Phugoid)
    ylabel('de [deg]')
    xlabel('t [s]')
%     title('Elevator Deflection - Phugoid')
    
    %Dutch Roll
    figure
    plot(time_DutchRoll, [RollRate_DutchRoll y_DutchRoll_init(:,3) y_DutchRoll_update(:,3)])
    legend({'Measurement', 'Initial values', 'Updated values'})
    ylabel('p [deg/s]')
    xlabel('t [s]')
%     title('Roll Rate - Dutch Roll')
    
    figure
    plot(time_DutchRoll, [YawRate_DutchRoll y_DutchRoll_init(:,4) y_DutchRoll_update(:,4)])
    legend({'Measurement', 'Initial values', 'Updated values'})
    ylabel('r [deg/s]')
    xlabel('t [s]')
%     title('Yaw Rate - Dutch Roll')

    figure
    plot([RollRate_DutchRoll y_DutchRoll_init(:,3) y_DutchRoll_update(:,3)], [YawRate_DutchRoll y_DutchRoll_init(:,4) y_DutchRoll_update(:,4)])
    legend({'Measurement', 'Initial values', 'Updated values'})
    ylabel('r [deg/s]')
    xlabel('p [deg/s]')
%     title('Roll rate vs Yaw rate - Dutch Roll')
    
    figure   
    plot(time_DutchRoll, [da_DutchRoll dr_DutchRoll])
    ylabel('deflection [deg]')
    xlabel('t [s]')
    legend({'Aileron','Rudder'})
%     title('Inputs - Dutch Roll')

    %Aperiodic Roll
    figure
    plot(time_AperiodicRoll, [Roll_AperiodicRoll y_AperiodicRoll_init(:,2) y_AperiodicRoll_update(:,2)])
    legend({'Measurement', 'Initial values', 'Updated values'})
    ylabel('\phi [deg/s]')
    xlabel('t [s]')
%     title('Roll Angle - Aperiodic Roll')
    
    figure
    plot(time_AperiodicRoll, [RollRate_AperiodicRoll y_AperiodicRoll_init(:,3) y_AperiodicRoll_update(:,3)])
    legend({'Measurement', 'Initial values', 'Updated values'})
    ylabel('p [deg/s]')
    xlabel('t [s]')
%     title('Roll Rate - Aperiodic Roll')
    
    figure
    plot(time_AperiodicRoll, [da_AperiodicRoll dr_AperiodicRoll])
    ylabel('deflection [deg]')
    xlabel('t [s]')
    legend({'Aileron','Rudder'})
%     title('Inputs - Aperiodic Roll')
    
    %Spiral
    figure
    plot(time_Spiral, [Roll_Spiral y_Spiral_init(:,2) y_Spiral_update(:,2)])
    legend({'Measurement', 'Initial values', 'Updated values'})
    ylabel('\phi [deg]')
    xlabel('t [s]')
%     title('Roll Angle - Spiral')
    
    figure
    plot(time_Spiral, [YawRate_Spiral y_Spiral_init(:,4) y_Spiral_update(:,4)])
    legend({'Measurement', 'Initial values', 'Updated values'})
    ylabel('r [deg/s]')
    xlabel('t [s]')
%     title('Yaw Rate - Spiral')
    
    figure
    plot(time_Spiral, [RollRate_Spiral y_Spiral_init(:,3) y_Spiral_update(:,3)])
    legend({'Measurement', 'Initial values', 'Updated values'})
    ylabel('p [deg/s]')
    xlabel('t [s]')
%     title('Roll Rate - Spiral')
    
    figure
    plot(time_Spiral, [da_Spiral dr_Spiral])
    ylabel('deflection [deg]')
    xlabel('t [s]')
    legend({'Aileron','Rudder'})
%     title('Inputs - Spiral')
end

if ResultPlots
    %Short Period
    figure
    plot(time_ShortPeriod, [y_ShortPeriod_init(:,1)])
    ylabel('V [m/s]')
    xlabel('t [s]')
    
    figure
    plot(time_ShortPeriod, [y_ShortPeriod_init(:,2)])
    ylabel('\alpha [deg]')
    xlabel('t [s]')
    
    figure
    plot(time_ShortPeriod, [y_ShortPeriod_init(:,3)])
    ylabel('\theta [deg]')
    xlabel('t [s]')
    
    figure
    plot(time_ShortPeriod, [y_ShortPeriod_init(:,4)])
    ylabel('q [deg/s]')
    xlabel('t [s]')

    %Phugoid    
    figure
    plot(time_Phugoid, [y_Phugoid_init(:,1)])
    ylabel('V [m/s]')
    xlabel('t [s]')
    
    figure
    plot(time_Phugoid, [y_Phugoid_init(:,2)])
    ylabel('\alpha [deg]')
    xlabel('t [s]')
    
    figure
    plot(time_Phugoid, [y_Phugoid_init(:,3)])
    ylabel('\theta [deg]')
    xlabel('t [s]')
    
    figure
    plot(time_Phugoid, [y_Phugoid_init(:,4)])
    ylabel('q [deg/s]')
    xlabel('t [s]')
    
    %Dutch Roll
    figure
    plot(time_DutchRoll, [y_DutchRoll_init(:,3)])
    ylabel('p [deg/s]')
    xlabel('t [s]')
    
    figure
    plot(time_DutchRoll, [y_DutchRoll_init(:,4)])
    ylabel('r [deg/s]')
    xlabel('t [s]')

    figure
    plot([y_DutchRoll_init(:,3)], [y_DutchRoll_init(:,4)])
    ylabel('r [deg/s]')
    xlabel('p [deg/s]')
    
    %Aperiodic Roll
    figure
    plot(time_AperiodicRoll, [y_AperiodicRoll_init(:,2)])
    ylabel('\phi [deg/s]')
    xlabel('t [s]')
    
    figure
    plot(time_AperiodicRoll, [y_AperiodicRoll_init(:,3)])
    ylabel('p [deg/s]')
    xlabel('t [s]')
    
    %Spiral
    figure
    plot(time_Spiral, [y_Spiral_init(:,2)])
    ylabel('\phi [deg]')
    xlabel('t [s]')
    
    figure
    plot(time_Spiral, [y_Spiral_init(:,4)])
    ylabel('r [deg/s]')
    xlabel('t [s]')
    
    figure
    plot(time_Spiral, [y_Spiral_init(:,3)])
    ylabel('p [deg/s]')
    xlabel('t [s]')
end