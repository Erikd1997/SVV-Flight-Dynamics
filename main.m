clc
clear all
close all

%Call the file containing all variables
Cit_par

%Save these values in a struct
save('Cit_par.mat')
c = load('Cit_par.mat');

%Create state-space-system
[sysS, sysA, A_S, A_A] = state_space_system(c);

%Create initial state vector
x0 = [c.V0, c.alpha0, c.th0, 0];

tS = 0:1:150;
uS = ones(1,length(tS))*-0.005;

tA = 0:1:150;
uA = [ones(1,length(tA))*-0.005; ones(1,length(tA))*-0.005];

figure(1)
lsim(sysS,uS,tS,x0)
figure(2)
lsim(sysA,uA,tA,x0)