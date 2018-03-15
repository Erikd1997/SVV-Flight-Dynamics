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
[sysS, sysA, A_S, A_A] = state_space_system(c);

%Create initial state vector
x0S = [ 0, 0, 0, 0];
x0A = [ 0, 0, 0, 0];

%Create time arrays
tS = 0:0.1:500;
tA = 0:0.1:500;

%Create input data
uS = deg2rad(sin(tS));                    %de = 15 degrees
uA = [deg2rad(sin(tA)); deg2rad(sin(tA))]; %da = 15degrees, dr = 26 degrees
%uA(:,1:20) = ones(2,20)*deg2rad(1);

%Evaluate state-space model
yS = lsim(sysS,uS,tS,x0S);
yA = lsim(sysA,uA,tA,x0A);

%[yS, tS] = step(sysS);
%[yA, tA] = step(sysA);

%Save and transform, if necessary
V  = yS(:,1)*c.V0+c.V0;
alpha = yS(:,2) + c.alpha0;
theta = yS(:,3) + c.th0;
q = yS(:,4)*c.V0/c.c;

beta = yA(:,1);
phi = yA(:,2);
p = yA(:,3)*2*c.V0/c.c;
r = yA(:,4)*2*c.V0/c.c;

%Plot everything!!!
%Start with symmetric case
figure(1)
subplot(4,1,1)
plot(tS, V)
title('Velocity')
subplot(4,1,2)
plot(tS, alpha)
title('Angle of attack')
subplot(4,1,3)
plot(tS,theta)
title('Pitch angle')
subplot(4,1,4)
plot(tS,q)
title('Pitch rate')

%Asymmetric case
figure(2)
subplot(4,1,1)
plot(tA, beta)
title('Yaw angle')
subplot(4,1,2)
plot(tA, phi)
title('Roll angle')
subplot(4,1,3)
plot(tA,p)
title('Yaw rate')
subplot(4,1,4)
plot(tA,r)
title('Roll rate')

