%clc
%close all

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

Roll = deg2rad(flight.flightdata.Ahrs1_Roll.data);
Pitch = deg2rad(flight.flightdata.Ahrs1_Pitch.data);
RollRate = deg2rad(flight.flightdata.Ahrs1_bRollRate.data);
PitchRate = deg2rad(flight.flightdata.Ahrs1_bPitchRate.data);
YawRate = deg2rad(flight.flightdata.Ahrs1_bYawRate.data);
AOA = deg2rad(flight.flightdata.vane_AOA.data);
V = flight.flightdata.Dadc1_tas.data*0.514444444;

%Read times of each manoeuvre from datasheet
t_DutchRoll = 43*60         *10;
t_DutchRoll_YD = 44*60      *10;
t_Phugoid = 48*60           *10;
t_ShortPeriod = 47*60       *10;
t_AperiodicRoll = 48*60     *10;
t_Spiral = 54*60            *10;

%Specify how long each manoeuvre takes
DR_length = 70              *10;
DR_YD_length = 70           *10;
Phugoid_length = 210        *10;
ShortPeriod_length = 65     *10;
AperiodicRoll_length = 90   *10;
Spiral_length = 120         *10;

%Dutch roll
Roll_DutchRoll = Roll(t_DutchRoll:t_DutchRoll+DR_length);
RollRate_DutchRoll = RollRate(t_DutchRoll:t_DutchRoll+DR_length);
YawRate_DutchRoll = YawRate(t_DutchRoll:t_DutchRoll+DR_length);
time_DutchRoll = time(t_DutchRoll:t_DutchRoll+DR_length);

%Phugoid
Roll_DutchRoll = Roll(t_DutchRoll:t_DutchRoll+DR_length);
RollRate_DutchRoll = RollRate(t_DutchRoll:t_DutchRoll+DR_length);
YawRate_DutchRoll = YawRate(t_DutchRoll:t_DutchRoll+DR_length);
time_DutchRoll = time(t_DutchRoll:t_DutchRoll+DR_length);

%
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