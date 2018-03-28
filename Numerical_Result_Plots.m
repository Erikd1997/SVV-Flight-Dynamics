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

[sysA, sysS] = state_space_system(c);