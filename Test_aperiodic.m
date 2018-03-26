clear
Cit_par_Aperiodic_Roll
%save('Cit_par.mat')
%clear
%load('Cit_par.mat')

%% Symmetric
%Phugoid 

A_p = 2*muc * (CZa * Cmq - 2 * muc * Cma);
B_p = -2*muc * (-CXu * Cma + Cmu * CXa) + Cmq * (CZu * CXa - CXu * CZa);
C_p = CZ0 * (Cmu * CZa - CZu * Cma);

lambda_p_plus = (-B_p + sqrt((B_p)^2 - 4* A_p * C_p)) / (2 * A_p);
lambda_p_min = (-B_p - sqrt((B_p)^2 - 4* A_p * C_p)) / (2 * A_p);


xi_p_plus = real(lambda_p_plus);
eta_p_plus = imag(lambda_p_plus);

xi_p_min = real(lambda_p_min);
eta_p_min = imag(lambda_p_min);

%Short Period
% A_sp = -2 * muc * KY2;
% B_sp = Cmadot + Cmq;
% C_sp = Cma;

A_sp = -2*muc*KY2*(CZadot-2*muc);
B_sp = (CZadot-2*muc)*Cmq -2*muc*KY2*CZa - Cmadot*(CZq+2*muc);
C_sp = CZa*Cmq-Cma*CZq-2*muc*Cma;

lambda_sp_plus = (-B_sp + sqrt((B_sp)^2 - 4*A_sp * C_sp)) / (2 * A_sp);
lambda_sp_min = (-B_sp - sqrt((B_sp)^2 - 4*A_sp * C_sp)) / (2 * A_sp);

xi_sp_plus = real(lambda_sp_plus);
eta_sp_plus = imag(lambda_sp_plus);

xi_sp_min = real(lambda_sp_min);
eta_sp_min = imag(lambda_sp_min);
%% Assymetric

%Aperiodic Roll

lambda_ar = Clp / (4 * mub * KX2)

xi_ar = real(lambda_ar);

% Dutch Roll
%lecture
% % A_dr=-2*mub*KZ2;
% % B_dr=Cnr/2;
% % C_dr=-Cnb;

%report
A_dr=8*mub*mub*KZ2;
B_dr=-2*mub*(Cnr+2*KZ2*CYb);
C_dr=4*mub*Cnb+CYb*Cnr;
lambda_dr_plus = (-B_dr + sqrt((B_dr)^2 - 4*A_dr * C_dr)) / (2 * A_dr);
lambda_dr_min = (-B_dr - sqrt((B_dr)^2 - 4*A_dr * C_dr)) / (2 * A_dr);
%lambda_dr =(Cnr + 2 * KZ2 * CYb) / (8 * mub * KZ2) real part

xi_dr_plus = real(lambda_dr_plus);
eta_dr_plus = imag(lambda_dr_plus);

xi_dr_min = real(lambda_dr_min);
eta_dr_min = imag(lambda_dr_min);


%Spiral

lambda_s = (2 * CL * (Clb * Cnr - Cnb * Clr)) / (Clp * (CYb * Cnr + 4 * mub * Cnb) - Cnp * (CYb * Clr + 4 * mub * Clb))

xi_s = real(lambda_s);

%% Substitution of lambda
% it was assumed that the velocity needed in the formula was the true
% airpseed in stationary flight

%Phugoid

P_p_plus = 2 * pi * c / ( eta_p_plus * V0)
P_p_min = 2 * pi * c / ( eta_p_min * V0)

zeta_p_plus = -xi_p_plus / abs(lambda_p_plus)
zeta_p_min = -xi_p_min / abs(lambda_p_min)

T_p_half_plus = log(1/2)*c / (xi_p_plus *V0)
T_p_half_min = log(1/2)*c / (xi_p_min *V0)

w_p_plus = abs(lambda_p_plus)*V0/c
w_p_min = abs(lambda_p_min)*V0/c


%Short Period

P_sp_plus = 2 * pi * c / ( eta_sp_plus * V0)
P_sp_min = 2 * pi * c / ( eta_sp_min * V0)

zeta_sp_plus = - xi_sp_plus / abs(lambda_sp_plus)
zeta_sp_min = - xi_sp_min / abs(lambda_sp_min)

T_sp_half_plus = log(1/2)*c / (xi_sp_plus *V0)
T_sp_half_min = log(1/2)*c / (xi_sp_min *V0)

w_sp_plus = abs(lambda_sp_plus)*V0/c
w_sp_min = abs(lambda_sp_min)*V0/c

%aperiodic roll
zeta_ar = -xi_ar / abs(lambda_ar)
T_half_ar = log(1/2)*b / (xi_ar *V0)
w_ar = abs(lambda_ar)*V0/b

%Dutch roll

P_dr_plus = 2 * pi * b / ( eta_dr_plus * V0)
P_dr_min = 2 * pi * b / ( eta_dr_min * V0)

zeta_dr_plus = - xi_dr_plus / abs(lambda_dr_plus)
zeta_dr_min = - xi_dr_min / abs(lambda_dr_min)

T_dr_half_plus = log(1/2)*b / (xi_dr_plus *V0)
T_dr_half_min = log(1/2)*b / (xi_dr_min *V0)

w_dr_plus = abs(lambda_dr_plus)*V0/b
w_dr_min = abs(lambda_dr_min)*V0/b

% zeta_dr = -xi_dr / abs(lambda_dr)
% T_half_dr = log(1/2)*b / (xi_dr *V0)
% w_dr = abs(lambda_dr)*V0/b

%spiral
zeta_s = -xi_s / abs(lambda_s)
T_half_s = log(1/2)*b / (xi_s *V0)
w_s = abs(lambda_s)*V0/b