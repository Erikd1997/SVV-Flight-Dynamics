%clc
close all
clear

%Call the file containing all variables
Cit_par

%Save these values in a struct
save('Cit_par.mat')

%Delete all saved variables from workspace
clear

%Load in necessary values in a struct
c = load('Cit_par.mat');
flight = load('FTISxprt-20180314_101817.mat');
datasheet = 'Post_Flight_Datasheet_Flight_2_DD_14_3_2018.xlsx';

%Which plots should appear?
plot_DutchRoll = false;
plot_Phugoid = false;
plot_ShortPeriod = true;
plot_AperiodicRoll = false;
plot_Spiral = false;

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

%% Solve simulations of each manouevre
%------Dutch Roll------
c.hp0    = hp(t_DutchRoll)*0.3048;
c.V0     = TAS(t_DutchRoll);
c.alpha0 = deg2rad(AOA(t_DutchRoll));
c.th0    = deg2rad(Pitch(t_DutchRoll));
c.rho    = c.rho0*((1+(c.lambda*c.hp0/c.Temp0)))^(-((c.g/(c.lambda*c.R))+1));
c.m      = Wi - FU_left(t_DutchRoll) - FU_right(t_DutchRoll);
c.W      = c.m*c.g;
c.muc    = c.m/(c.rho*c.S*c.c);
c.mub    = c.m/(c.rho*c.S*c.b);
c.CL     = 2*c.W/(c.rho*c.V0^2*c.S);
c.CD     = c.CD0 + (c.alpha0*c.CLa)^2/(pi*c.A*c.e);
c.CX0    = c.W*sin(c.th0)/(0.5*c.rho*c.V0^2*c.S);
c.CZ0    = -c.W*cos(c.th0)/(0.5*c.rho*c.V0^2*c.S);

%Create state-space-system
[~, sysA_DutchRoll] = state_space_system(c);

%Save and non-dimensionalize eigenvalues
eig_DutchRoll = eig(sysA_DutchRoll)*c.b/c.V0;

%Initial condition
x0 = deg2rad([0, Roll_DutchRoll(1), RollRate_DutchRoll(1), YawRate_DutchRoll(1)]);
%Simulate using state space system
y_DutchRoll = lsim(sysA_DutchRoll, deg2rad([da_DutchRoll -dr_DutchRoll]), time_DutchRoll, x0);

% Transform
y_DutchRoll(:,1) = rad2deg(y_DutchRoll(:,1));
y_DutchRoll(:,2) = rad2deg(y_DutchRoll(:,2));
y_DutchRoll(:,3) = rad2deg(y_DutchRoll(:,3));
y_DutchRoll(:,4) = rad2deg(y_DutchRoll(:,4));

%------Phugoid------
c.hp0 = hp(t_Phugoid)*0.3048;
c.V0 = TAS(t_Phugoid);
c.alpha0 = deg2rad(AOA(t_Phugoid));
c.th0 = deg2rad(Pitch(t_Phugoid));
c.rho = c.rho0*((1+(c.lambda*c.hp0/c.Temp0)))^(-((c.g/(c.lambda*c.R))+1));
c.m = Wi - FU_left(t_Phugoid) - FU_right(t_Phugoid);
c.W = c.m*c.g;
c.muc    = c.m/(c.rho*c.S*c.c);
c.mub    = c.m/(c.rho*c.S*c.b);
c.CL = 2*c.W/(c.rho*c.V0^2*c.S);
c.CD = c.CD0 + (c.alpha0*c.CLa)^2/(pi*c.A*c.e);
c.CX0    = c.W*sin(c.th0)/(0.5*c.rho*c.V0^2*c.S);
c.CZ0    = -c.W*cos(c.th0)/(0.5*c.rho*c.V0^2*c.S);

%Create state-space-system
[sysS_Phugoid, ~] = state_space_system(c);

%Save and non-dimensionalize eigenvalues
eig_Phugoid = eig(sysS_Phugoid)*c.c/c.V0;

%Simulate using state space system
y_Phugoid = lsim(sysS_Phugoid, deg2rad(de_Phugoid), time_Phugoid);
% Transform
y_Phugoid(:,1) = y_Phugoid(:,1)+c.V0;
y_Phugoid(:,2) = rad2deg(y_Phugoid(:,2) + c.alpha0);
y_Phugoid(:,3) = rad2deg(y_Phugoid(:,3) + c.th0);
y_Phugoid(:,4) = rad2deg(y_Phugoid(:,4));

%------Short Period------
c.hp0 = hp(t_ShortPeriod)*0.3048;
c.V0 = TAS(t_ShortPeriod);
c.alpha0 = deg2rad(AOA(t_ShortPeriod));
c.th0 = deg2rad(Pitch(t_ShortPeriod));
c.rho = c.rho0*((1+(c.lambda*c.hp0/c.Temp0)))^(-((c.g/(c.lambda*c.R))+1));
c.m = Wi - FU_left(t_ShortPeriod) - FU_right(t_ShortPeriod);
c.W = c.m*c.g;
c.muc    = c.m/(c.rho*c.S*c.c);
c.mub    = c.m/(c.rho*c.S*c.b);
c.CL = 2*c.W/(c.rho*c.V0^2*c.S);
c.CD = c.CD0 + (c.alpha0*c.CLa)^2/(pi*c.A*c.e);
c.CX0    = c.W*sin(c.th0)/(0.5*c.rho*c.V0^2*c.S);
c.CZ0    = -c.W*cos(c.th0)/(0.5*c.rho*c.V0^2*c.S);

%Create state-space-system
[sysS_ShortPeriod, ~] = state_space_system(c);

%Save and non-dimensionalize eigenvalues
eig_ShortPeriod = eig(sysS_ShortPeriod)*c.c/c.V0;

%Simulate using state space system
y_ShortPeriod = lsim(sysS_ShortPeriod, deg2rad(de_ShortPeriod), time_ShortPeriod);

% Transform
y_ShortPeriod(:,1) = y_ShortPeriod(:,1)+c.V0;
y_ShortPeriod(:,2) = rad2deg(y_ShortPeriod(:,2) + c.alpha0);
y_ShortPeriod(:,3) = rad2deg(y_ShortPeriod(:,3) + c.th0);
y_ShortPeriod(:,4) = rad2deg(y_ShortPeriod(:,4));

%------Aperiodic Roll------
c.hp0 = hp(t_AperiodicRoll)*0.3048;
c.V0 = TAS(t_AperiodicRoll);
c.alpha0 = deg2rad(AOA(t_AperiodicRoll));
c.th0 = deg2rad(Pitch(t_AperiodicRoll));
c.rho = c.rho0*((1+(c.lambda*c.hp0/c.Temp0)))^(-((c.g/(c.lambda*c.R))+1));
c.m = Wi - FU_left(t_AperiodicRoll) - FU_right(t_AperiodicRoll);
c.W = c.m*c.g;
c.muc    = c.m/(c.rho*c.S*c.c);
c.mub    = c.m/(c.rho*c.S*c.b);
c.CL = 2*c.W/(c.rho*c.V0^2*c.S);
c.CD = c.CD0 + (c.alpha0*c.CLa)^2/(pi*c.A*c.e);
c.CX0    = c.W*sin(c.th0)/(0.5*c.rho*c.V0^2*c.S);
c.CZ0    = -c.W*cos(c.th0)/(0.5*c.rho*c.V0^2*c.S);

%Create state-space-system
[~, sysA_AperiodicRoll] = state_space_system(c);

%Initial condition
x0 = [0, deg2rad(Roll_AperiodicRoll(1)), deg2rad(RollRate_AperiodicRoll(1)), deg2rad(YawRate_AperiodicRoll(1))];

%Save and non-dimensionalize eigenvalues
eig_AperiodicRoll = eig(sysA_AperiodicRoll)*c.b/c.V0;

%Simulate using state space system
y_AperiodicRoll = lsim(sysA_AperiodicRoll, deg2rad([da_AperiodicRoll-da_AperiodicRoll(1) dr_AperiodicRoll-dr_AperiodicRoll(1)]), time_AperiodicRoll, x0);

% Transform
y_AperiodicRoll(:,1) = rad2deg(y_AperiodicRoll(:,1));
y_AperiodicRoll(:,2) = rad2deg(y_AperiodicRoll(:,2));
y_AperiodicRoll(:,3) = rad2deg(y_AperiodicRoll(:,3));
y_AperiodicRoll(:,4) = rad2deg(y_AperiodicRoll(:,4));

%------Spiral------
c.hp0 = hp(t_Spiral)*0.3048;
c.V0 = TAS(t_Spiral);
c.alpha0 = deg2rad(AOA(t_Spiral));
c.th0 = deg2rad(Pitch(t_Spiral));
c.rho = c.rho0*((1+(c.lambda*c.hp0/c.Temp0)))^(-((c.g/(c.lambda*c.R))+1));
c.m = Wi - FU_left(t_Spiral) - FU_right(t_Spiral);
c.W = c.m*c.g;
c.muc    = c.m/(c.rho*c.S*c.c);
c.mub    = c.m/(c.rho*c.S*c.b);
c.CL = 2*c.W/(c.rho*c.V0^2*c.S);
c.CD = c.CD0 + (c.alpha0*c.CLa)^2/(pi*c.A*c.e);
c.CX0    = c.W*sin(c.th0)/(0.5*c.rho*c.V0^2*c.S);
c.CZ0    = -c.W*cos(c.th0)/(0.5*c.rho*c.V0^2*c.S);

%Create state-space-system
[~, sysA_Spiral] = state_space_system(c);

%Initial state
x0_Spiral = [0, deg2rad(Roll(t_Spiral)), deg2rad(RollRate(t_Spiral)), deg2rad(YawRate(t_Spiral))];

%Save and non-dimensionalize eigenvalues
eig_Spiral = eig(sysA_Spiral)*c.b/c.V0;

%Simulate using state space system
y_Spiral = lsim(sysA_Spiral, deg2rad([da_Spiral dr_Spiral]), time_Spiral, x0_Spiral);

% Transform
y_Spiral(:,1) = rad2deg(y_Spiral(:,1));
y_Spiral(:,2) = rad2deg(y_Spiral(:,2));
y_Spiral(:,3) = rad2deg(y_Spiral(:,3));
y_Spiral(:,4) = rad2deg(y_Spiral(:,4));

%% Plot Measurements
%Dutch Roll
if plot_DutchRoll
    figure
    subplot(4,1,1)
    plot(time_DutchRoll, [RollRate_DutchRoll y_DutchRoll(:,3)])
    ylabel('r (deg/s)')
    xlabel('t (s)')
    title('Roll Rate - Dutch Roll')
    subplot(4,1,2)
    plot(time_DutchRoll, [YawRate_DutchRoll y_DutchRoll(:,4)])
    ylabel('p (deg/s)')
    xlabel('t (s)')
    title('Yaw Rate - Dutch Roll')
    subplot(4,1,3)
    plot([RollRate_DutchRoll y_DutchRoll(:,3)], [YawRate_DutchRoll y_DutchRoll(:,4)])
    ylabel('p (deg/s)')
    xlabel('r (deg/s)')
    title('Roll rate vs Yaw rate - Dutch Roll')
    subplot(4,1,4)
    plot(time_DutchRoll, [da_DutchRoll dr_DutchRoll])
    xlabel('t (s)')
    ylabel('da (deg)')
    xlabel('t (s)')
    title('Aileron Deflection (blue), Rudder Deflection (red) - Dutch Roll')
end

%Phugoid
if plot_Phugoid
    figure
    subplot(5,1,1)
    plot(time_Phugoid, [TAS_Phugoid y_Phugoid(:,1)])
    ylabel('V (m/s)')
    xlabel('t (s)')
    title('True Airspeed - Phugoid')
    subplot(5,1,2)
    plot(time_Phugoid, [AOA_Phugoid y_Phugoid(:,2)])
    ylabel('\alpha (deg)')
    xlabel('t (s)')
    title('Angle of Attack - Phugoid')
    subplot(5,1,3)
    plot(time_Phugoid, [Pitch_Phugoid y_Phugoid(:,3)])
    ylabel('\theta (deg)')
    xlabel('t (s)')
    title('Pitch angle - Phugoid')
    subplot(5,1,4)
    plot(time_Phugoid, [PitchRate_Phugoid y_Phugoid(:,4)])
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
    figure
    subplot(5,1,1)
    plot(time_ShortPeriod, [TAS_ShortPeriod y_ShortPeriod(:,1)])
    ylabel('V (m/s)')
    xlabel('t (s)')
    title('True Airspeed - Short Period')
    subplot(5,1,2)
    plot(time_ShortPeriod, [AOA_ShortPeriod y_ShortPeriod(:,2)])
    ylabel('\alpha (deg)')
    xlabel('t (s)')
    title('Angle of Attack - Short Period')
    subplot(5,1,3)
    plot(time_ShortPeriod, [Pitch_ShortPeriod y_ShortPeriod(:,3)])
    ylabel('\theta (deg)')
    xlabel('t (s)')
    title('Pitch Angle - Short Period')
    subplot(5,1,4)
    plot(time_ShortPeriod, [PitchRate_ShortPeriod y_ShortPeriod(:,4)])
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
    figure
    subplot(3,1,1)
    plot(time_AperiodicRoll, [Roll_AperiodicRoll y_AperiodicRoll(:,2)])
    ylabel('\phi (deg/s)')
    xlabel('t (s)')
    title('Roll Angle - Aperiodic Roll')
    subplot(3,1,2)
    plot(time_AperiodicRoll, [RollRate_AperiodicRoll y_AperiodicRoll(:,3)])
    ylabel('r (deg/s)')
    xlabel('t (s)')
    title('Roll Rate - Aperiodic Roll')
    subplot(3,1,3)
    plot(time_AperiodicRoll, da_AperiodicRoll)
    ylabel('da (deg)')
    xlabel('t (s)')
    title('Aileron Deflection - Aperiodic Roll')
end

%Spiral
if plot_Spiral
    figure
    subplot(4,1,1)
    plot(time_Spiral, [Roll_Spiral y_Spiral(:,2)])
    ylabel('\phi (deg)')
    xlabel('t (s)')
    title('Roll Angle - Spiral')
    subplot(4,1,2)
    plot(time_Spiral, [YawRate_Spiral y_Spiral(:,4)])
    ylabel('p (deg/s)')
    xlabel('t (s)')
    title('Yaw Rate - Spiral')
    subplot(4,1,3)
    plot(time_Spiral, [RollRate_Spiral y_Spiral(:,3)])
    ylabel('r (deg/s)')
    xlabel('t (s)')
    title('Roll Rate - Spiral')
    subplot(4,1,4)
    plot(time_Spiral, [da_Spiral dr_Spiral])
    ylabel('da (deg)')
    xlabel('t (s)')
    title('Aileron Deflection (blue), Rudder Deflection (red) - Spiral')
end