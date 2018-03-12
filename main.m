%clc
%close all

%Call the file containing all variables
Cit_par

%Save these values in a struct
save('Cit_par.mat')

%Delete all saved variables from workspace
clear

%Load in values in a struct 'c'
c = load('Cit_par.mat');

%Create state-space-system
[sysS, sysA] = state_space_system(c);

%Create initial state vector
x0S = [c.V0,  c.alpha0, c.th0, 0];
x0A = [  0,       0,      0,   0];

tS = 0:0.1:500;
uS = deg2rad(sin(tS));                    %de = 15 degrees

tA = 0:0.1:500;
uA = [deg2rad(sin(tA)); deg2rad(sin(tA))]; %da = 15degrees, dr = 26 degrees
%uA(:,1:20) = ones(2,20)*deg2rad(1);

figure(1)
lsim(sysS,uS,tS,x0S)
figure(2)
lsim(sysA,uA,tA,x0A)