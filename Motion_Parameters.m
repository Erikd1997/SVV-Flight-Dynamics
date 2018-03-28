function [P, xi, T1_2, omega_0, omega_n] = Motion_Parameters(lambda, c, V)
P = 2*pi/imag(lambda)*c/V;
xi = -real(lambda)/abs(lambda);
T1_2 = log(1/2)/real(lambda)*c/V;
omega_0 = abs(lambda)*V/c;
omega_n = omega_0*sqrt(1-xi);