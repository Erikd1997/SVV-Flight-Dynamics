% Citation 550 - Linear simulation

% xcg = 0.25*c

% Stationary flight condition

datasheet = 'Post_Flight_Datasheet_Flight_2_DD_14_3_2018.xlsx';
flight = load('FTISxprt-20180314_101817.mat');
fd = flight.flightdata;

%First save necessary data
da = fd.delta_a.data;
de = fd.delta_e.data;
dr = fd.delta_r.data;
time = fd.time.data';

Roll = fd.Ahrs1_Roll.data;                              %[rad]
Pitch = fd.Ahrs1_Pitch.data;                            %[rad]
RollRate = fd.Ahrs1_bRollRate.data;                     %[rad/s]
PitchRate = fd.Ahrs1_bPitchRate.data;                   %[rad/s]
YawRate = fd.Ahrs1_bYawRate.data;                       %[rad/s]
AOA = fd.vane_AOA.data;                                 %[rad]
TAS = fd.Dadc1_tas.data*0.514444444;                    %[m/s]
hp = fd.Dadc1_alt.data*0.3048;                          %[m]
FMF_left = fd.lh_engine_FMF.data/3600;                  %[kg]
FMF_right = fd.rh_engine_FMF.data/3600;                 %[kg]
Wi = (xlsread(datasheet,'D18:D18')+9165)*0.45359237...  %[kg]
    + sum(xlsread(datasheet,'H8:H16'));

t_Spiral = (54*60+45)       *10;
Spiral_length = 100         *10;

FU = zeros(length(time),1);
for i = 2:length(time)
    FU(i)   = trapz(time(1:i),FMF_left(1:i)) + trapz(time(1:i),FMF_right(1:i)); %[kg]
end

% Stationary flight condition
FU_Spiral = FU(t_Spiral);     % [kg]
hp0    = hp(t_Spiral);               % pressure altitude in the stationary flight condition [m]
V0     = TAS(t_Spiral);              % true airspeed in the stationary flight condition [m/sec]
alpha0 = deg2rad(AOA(t_Spiral));     % angle of attack in the stationary flight condition [rad]
th0    = deg2rad(Pitch(t_Spiral));   % pitch angle in the stationary flight condition [rad]
m      = Wi - FU_Spiral;         % [kg]

% aerodynamic properties
e      = 0.7857;            % Oswald factor [ ]
CD0    = 0.022048;          % Zero lift drag coefficient [ ]
CLa    = 4.481;            % Slope of CL-alpha curve [ ]

% Longitudinal stability
Cma    = -0.5524;            % longitudinal stabilty [ ]
Cmde   = -1.2669;            % elevator effectiveness [ ]

% Aircraft geometry

S      = 30.00;	          % wing area [m^2]
Sh     = 0.2*S;           % stabiliser area [m^2]
Sh_S   = Sh/S;	          % [ ]
lh     = 0.71*5.968;      % tail length [m]
c      = 2.0569;          % mean aerodynamic cord [m]
lh_c   = lh/c;	          % [ ]
b      = 15.911;	      % wing span [m]
bh     = 5.791;	          % stabilser span [m]
A      = b^2/S;           % wing aspect ratio [ ]
Ah     = bh^2/Sh;         % stabilser aspect ratio [ ]
Vh_V   = 1;		          % [ ]
ih     = -2*pi/180;       % stabiliser angle of incidence [rad]

% Constant values concerning atmosphere and gravity

rho0   = 1.2250;          % air density at sea level [kg/m^3] 
lambda = -0.0065;         % temperature gradient in ISA [K/m]
Temp0  = 288.15;          % temperature at sea level in ISA [K]
R      = 287.05;          % specific gas constant [m^2/sec^2K]
g      = 9.81;            % [m/sec^2] (gravity constant)
p0     = 101325;          % [Pa]

rho    = rho0*((1+(lambda*hp0/Temp0)))^(-((g/(lambda*R))+1));   % [kg/m^3]  (air density)
W      = m*g;

% Constant values concerning aircraft inertia

muc    = m/(rho*S*c);
mub    = m/(rho*S*b);
KX2    = 0.019;
KZ2    = 0.042;
KXZ    = 0.002;
KY2    = 1.25*1.114;

% Aerodynamic constants

Cmac   = 0;                     % Moment coefficient about the aerodynamic centre [ ]
CNwa   = CLa;   		        % Wing normal force slope [ ]
CNha   = 2*pi*Ah/(Ah+2);        % Stabiliser normal force slope [ ]
depsda = 4/(A+2);               % Downwash gradient [ ]

% Lift and drag coefficient

CL = 2*W/(rho*V0^2*S);               % Lift coefficient [ ]
CD = CD0 + (CLa*alpha0)^2/(pi*A*e);  % Drag coefficient [ ]

% Stabiblity derivatives

CX0    = W*sin(th0)/(0.5*rho*V0^2*S);
CXu    = -0.02792;
CXa    = -0.47966;
CXadot = +0.08330;
CXq    = -0.28170;
CXde   = -0.03728;

CZ0    = -W*cos(th0)/(0.5*rho*V0^2*S);
CZu    = -0.37616;
CZa    = -5.74340;
CZadot = -0.00350;
CZq    = -5.66290;
CZde   = -0.69612;

Cmu    = +0.06990;
Cmadot = +0.17800;
Cmq    = -8.79415;

CYb    = -0.7500;
CYbdot =  0     ;
CYp    = -0.0304;
CYr    = +0.8495;
CYda   = +0.0400;
CYdr   = +0.2300;

Clb    = -0.10260;
Clp    = -0.71085;
Clr    = +0.23760;
Clda   = -0.23088;
Cldr   = +0.03440;

Cnb    = +0.1348;
Cnbdot = 0     ;
Cnp    = -0.0602;
Cnr    = -0.2061;
Cnda   = -0.0120;
Cndr   = -0.0939;