%Call the file containing all variables
Cit_par

%Save these values in a struct
save('Cit_par.mat')
c = load('Cit_par.mat');

%Create state-space-system
[sysS, sysA, A_S, A_A] = state_space_system(c);