%Call the file containing all variables
Cit_par

%Save these values in a struct
save('Cit_par.mat')
c = load('Cit_par.mat');

%Create state-space-system
[sysS, sysA] = state_space_system(c);

%Create initial state vector
x0 = [c.V0, c.alpha0, c.th0, 0];

t = 0:1:150;
u = ones(1,length(t))*-0.005;