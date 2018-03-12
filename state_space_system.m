function [sysS, sysA] = state_space_system(c)
%This function transforms the C1, C2 and C3 matrices of the linear model to
%a state space representation for both the symmetric and asymmetric model.

%% Symmetric model
% Defining C1, C2 and C3
C1_S = [-2*c.muc*c.c/c.V0^2,             0,               0,               0             ;
               0,           (c.CZadot-2*c.muc)*c.c/c.V0,  0,               0             ;
               0,                        0,               1,               0             ;
               0                 c.Cmadot*c.c/c.V0,       0,  -2*c.muc*c.KY2*(c.c/c.V0)^2];
           
C2_S = [c.CXu/c.V0, c.CXa,  c.CZ0,      c.CXq*c.c/c.V0     ;
        c.CZu/c.V0, c.CZa, -c.CX0, c.c/c.V0*(c.CZq+2*c.muc);
             0,       0,      0,             -1            ;
        c.Cmu/c.V0, c.Cma,    0,       c.Cmq*c.c/c.V0      ];

C3_S = [c.CXde;
        c.CZde;
          0 ; 
        c.Cmde];

% Create A, B, C and D matrices
A_S = -C1_S\C2_S;
B_S = -C1_S\C3_S;
C_S = eye(4);
D_S = 0;
    
%Create state space system
sysS = ss(A_S, B_S, C_S, D_S);
sysS.OutputName = {'u', '\alpha', '\theta', 'q'};

%% Asymmetric model
% Defining C1, C2 and C3
C1_A = [(c.CYbdot-2*c.mub)*c.b/c.V0, 0,                  0,                          0             ;
                    0,               1,                  0,                          0             ;
                    0,               0,   -2*c.mub*c.KX2*(c.b/c.V0)^2,   2*c.mub*c.KXZ*(c.b/c.V0)^2;
           c.Cnbdot*c.b/c.V0,        0,    2*c.mub*c.KXZ*(c.b/c.V0)^2,  -2*c.mub*c.KZ2*(c.b/c.V0)^2];
            
C2_A = [c.CYb, c.CL, c.CYp*c.b/(2*c.V0), c.b/(2*c.V0)*(c.CYr-4*c.mub);
          0,     0,        -1,                       0               ;
        c.Clb,   0,  c.Clp*c.b/(2*c.V0),      c.Clr*c.b/(2*c.V0)     ;
        c.Cnb,   0,  c.Cnp*c.b/(2*c.V0),      c.Cnr*c.b/(2*c.V0)     ];

C3_A = [c.CYda, c.CYdr;
          0,       0  ;
        c.Clda, c.Cldr;
        c.Cnda, c.Cndr];

% Create A, B, C and D matrices
A_A = -C1_A\C2_A;
B_A = -C1_A\C3_A;
C_A =    eye(4);
D_A =      0;

%Create state space system
sysA = ss(A_A, B_A, C_A, D_A);
sysA.OutputName = {'\beta', '\phi', 'p', 'r'};
end