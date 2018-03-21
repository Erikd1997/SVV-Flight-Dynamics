clc
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

%% Create state-space-system
[sysS, sysA, A_S, A_A] = state_space_system(c);

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
AOA = deg2rad(flight.flightdata.vane_AOA.data);
TAS = flight.flightdata.Dadc1_tas.data*0.514444444;

%Read times of each manoeuvre from datasheet
t_DutchRoll = (44*60+43.7)  *10;
t_DutchRoll_YD = 44*60      *10;
t_Phugoid = (49*60+10)      *10;
t_ShortPeriod = (47*60+50.5)*10;
t_AperiodicRoll = (46*60+57)*10;
t_Spiral = (54*60+45)       *10;

%Specify how long each manoeuvre takes
DR_length = 20              *10;
DR_YD_length = 70           *10;
Phugoid_length = 200        *10;
ShortPeriod_length = 10     *10;
AperiodicRoll_length = 30   *10;
Spiral_length = 100         *10;

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

%% Plot Measurements
%Dutch Roll
figure
subplot(3,1,1)
plot(time_DutchRoll, RollRate_DutchRoll)
ylabel('r (deg/s)')
xlabel('t (s)')
title('Roll Rate - Dutch Roll')
subplot(3,1,2)
plot(time_DutchRoll, YawRate_DutchRoll)
ylabel('p (deg/s)')
xlabel('t (s)')
title('Yaw Rate - Dutch Roll')
subplot(3,1,3)
plot(time_DutchRoll, da_DutchRoll)
xlabel('t (s)')
ylabel('da (deg)')
xlabel('t (s)')
title('Aileron Deflection - Dutch Roll')

%Phugoid
figure
subplot(5,1,1)
plot(time_Phugoid, TAS_Phugoid)
ylabel('V (m/s)')
xlabel('t (s)')
title('True Airspeed - Phugoid')
subplot(5,1,2)
plot(time_Phugoid, AOA_Phugoid)
ylabel('\alpha (deg)')
xlabel('t (s)')
title('Angle of Attack - Phugoid')
subplot(5,1,3)
plot(time_Phugoid, Pitch_Phugoid)
ylabel('\theta (deg)')
xlabel('t (s)')
title('Pitch angle - Phugoid')
subplot(5,1,4)
plot(time_Phugoid, PitchRate_Phugoid)
ylabel('q (deg/s)')
xlabel('t (s)')
title('Pitch Rate - Phugoid')
subplot(5,1,5)
plot(time_Phugoid, de_Phugoid)
ylabel('de (deg)')
xlabel('t (s)')
title('Elevator Deflection - Phugoid')

%Short Period
figure
subplot(5,1,1)
plot(time_ShortPeriod, TAS_ShortPeriod)
ylabel('V (m/s)')
xlabel('t (s)')
title('True Airspeed - Short Period')
subplot(5,1,2)
plot(time_ShortPeriod, AOA_ShortPeriod)
ylabel('\alpha (deg)')
xlabel('t (s)')
title('Angle of Attack - Short Period')
subplot(5,1,3)
plot(time_ShortPeriod, Pitch_ShortPeriod)
ylabel('\theta (deg)')
xlabel('t (s)')
title('Pitch Angle - Short Period')
subplot(5,1,4)
plot(time_ShortPeriod, PitchRate_ShortPeriod)
ylabel('q (deg/s)')
xlabel('t (s)')
title('Pitch Rate - Short Period')
subplot(5,1,5)
plot(time_ShortPeriod, de_ShortPeriod)
ylabel('de (deg)')
xlabel('t (s)')
title('Elevator Deflection - Short Period')

%Aperiodic Roll
figure
subplot(3,1,1)
plot(time_AperiodicRoll, Roll_AperiodicRoll)
ylabel('r (deg/s)')
xlabel('t (s)')
title('Roll Angle - Aperiodic Roll')
subplot(3,1,2)
plot(time_AperiodicRoll, RollRate_AperiodicRoll)
ylabel('p (deg/s)')
xlabel('t (s)')
title('Roll Rate - Aperiodic Roll')
subplot(3,1,3)
plot(time_AperiodicRoll, da_AperiodicRoll)
ylabel('da (deg)')
xlabel('t (s)')
title('Aileron Deflection - Aperiodic Roll')

%% Transform
V  = yS(:,1)*c.V0+c.V0;
alpha = yS(:,2) + c.alpha0;
theta = yS(:,3) + c.th0;
q = yS(:,4)*c.V0/c.c;

beta = yA(:,1);
phi = yA(:,2);
p = yA(:,3)*2*c.V0/c.c;
r = yA(:,4)*2*c.V0/c.c;

% %Plot everything!!!
% %Start with symmetric case
% figure(1)
% subplot(4,1,1)
% plot(tS, V)
% title('Velocity')
% subplot(4,1,2)
% plot(tS, alpha)
% title('Angle of attack')
% subplot(4,1,3)
% plot(tS,theta)
% title('Pitch angle')
% subplot(4,1,4)
% plot(tS,q)
% title('Pitch rate')
% 
% %Asymmetric case
% figure(2)
% subplot(4,1,1)
% plot(tA, beta)
% title('Yaw angle')
% subplot(4,1,2)
% plot(tA, phi)
% title('Roll angle')
% subplot(4,1,3)
% plot(tA,p)
% title('Yaw rate')
% subplot(4,1,4)
% plot(tA,r)
% title('Roll rate')