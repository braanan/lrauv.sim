% vehicle_coeffs.m
% July 14, 2014


global zg Mqq xg


% Mass Properties:
rho         = 1025;              % kg/m3   
g           = 9.80665;           % m/s2
mass        = 147.8671;          % kg Flooded Vehicle mass
% mass        = .110*rho;
volume      = 0.144260585365854; % m3 (equals 1450 N buoyancy in 1025 kg/m3 water)
% excludes buoyancy bladder at default setting
m = mass;           % kg, mass
W = m*g ;           % N, Weight                         % B = mass*g; 
B = rho*volume*g;   % N, Buoyancy

% Geometric Parameters (Used only by the simulation):
% rG - vehicle centers of gravity
% xg =  0.0;          % m
yg =  0.0;          % m ***-0.000236***
% zg =  0.0067940;    % m ***0.0067940***



% rB - vehicle centers of buoyancy
xb =  0.0;          % m ***0.1181***
yb =  0.0;          % m
zb =  0.0;          % m


Wp = W ;
Bp = B ;
Xgp = xg;
Ygp = yg;
Zgp = zg;
Xbp = xb;
Ybp = yb;
Zbp = zb;

%{
% Dropweight1 parameters
dropWt1Volume   = 0.0000881;    % m3    Volume of the drop weight #1, m3
dropWt1Mass     = 1.0;          % kg	Mass of the drop weight #1, kg
dropWt1X        = 0.1655;       % m     X location of the drop weight #1, m
dropWt1Y        = 0.0;          % m     Y location of the drop weight #1, m
dropWt1Z        = -0.20;        % m 	Z location of the drop weight #1, m

%  Zero these for now.  Fix this later in VehicleConstants.cc
Wdw1 = 0.0;
Bdw1 = 0.0;
   
Wp = W + Wdw1;
Bp = B + Bdw1;

Xgp = (xg*W + dropWt1X*Wdw1) / Wp;
Ygp = (yg*W + dropWt1Y*Wdw1) / Wp;
Zgp = (zg*W + dropWt1Z*Wdw1) / Wp;
Xbp = (xb*B + dropWt1X*Bdw1) / Bp;
Ybp = (yb*B + dropWt1Y*Bdw1) / Bp;
Zbp = (zb*B + dropWt1Z*Bdw1) / Bp;
%}

Sfin    =   1.15e-2;        % m^2     Total area of elevator = 2 x fin.
bfin    =   18.57e-2;       % m       Fin span 
zfin    =   0.152;          % m       Centerline to fin  
xfin    =  -0.633;          % m       Midpoint to elevator axle (x)
% xf      =   xfin+2.33e-2;   % m       Aft end to Fin section
% xf2     =   xfin-3.84e-2;   % m       Forward end of fin section  
% a       =   bfin+zfin; % m mean fin height above center line
% Kpdot   = -(xf-xf2)*2*rho*((a)^4)/pi;
% (3.84/8.5)*(zfin + bfin) + (4.66/8.5)*(9.3e-2 + zfin)

% Mass Properties:
Ixx =  3.000000;    % kg-m2     Diagonal inertia tensor
Iyy =  41.980233;   % kg-m2     Diagonal inertia tensor
Izz =  41.980233;   % kg-m2     Diagonal inertia tensor

% Thruster parameters:
Kpp   = -0.191601;    % kg-m2*  Rolling Resistance             *kg-m2/rad2?
Kprop =  0.23;        % N-m     Propeller Torque ***0.23***

% Added Mass:
Yvdot  = -126.324739; % kg;     // Yvdot, kg.
Zwdot  = -126.324739; % kg;     // Zwdot, kg.
Xudot  =   -4.876161; % kg;     // Xudot, kg.
Mqdot  =  -33.463086; % kg-m2;  // Mqdot, kg-m^2.
Nrdot  =  -33.463086; % kg-m2;  // Nrdot, kg-m^2.
Kpdot  =    0.000000; % kg-m2;  // Kpdot, kg-m^2.
Kvdot  =    0.000000; % kg-m;   // Kvdot, kg-m.
Mwdot  =    7.117842; % kg-m;   // Mwdot, kg-m.
Zqdot  =    7.117842; % kg-m;   // Zqdot, kg-m.
Nvdot  =   -7.117842; % kg-m;   // Nvdot, kg-m.
Yrdot  =   -7.117842; % kg-m;   // Yrdot, kg-m.
Ypdot  =    0.000000; % kg-m;   // Ypdot, kg-m.

% Stability Derivatives:
Xqq =  7.117842;    % kg-m;
Xrr =  7.117842;    % kg-m;
Xvv = -54.370919;   % kg/m;   // Xvv   , kg/m
Xww = -54.370919;   % kg/m;   // Xww   , kg/m
Yvv = -601.274653;  % kg/m      Cross-flow Drag (Yv|v|)
Yrr =  0.000000;    % n/a*      Cross-flow Drag (Yr|r|)         *kg-m/rad2?
Zww = -601.274653;  % kg/m      Cross-flow Drag
Zqq =  0.000000;    % n/a*      Cross-flow Drag (Zq|q|)
Mww = -58.113144;   % kg        Cross-flow drag (-Nv|v|)
% Mqq = -632.698957;  % kg-m2*    Cross-flow drag (Mq|q|)        *kg-m2/rad2?
Nvv =  58.113144;   % kg     	Cross-flow drag (Nv|v|)
Nrr = -632.698957;  % kg-m2* 	Cross-flow drag (Nr|r|)        *kg-m2/rad2?

Yuv = -23.954759;   % kg/m      Body Lift Force and Fin Lift
Zuw = -23.954759;   % kg/m      Body Lift Force and Fin Lift
Nuv = -105.660262;  % kg        Body and Fin Lift and Munk Moment
Muw =  105.660262;  % kg        Body and Fin Lift and Munk Moment

Xwq = -126.324739;  % kg;
Xvr =  126.324739;  % kg;
Yur =  8.719853;    % kg*       Added Mass Cross-term and Fin Lift *kg/rad?
Zuq = -8.719853;    % kg*       Added Mass Cross-term and Fin Lift *kg/rad?
Nur = -61.182063;   % kg-m*     Added Mass Cross-term and Fin Lift *kg-m/rad?
Muq = -61.182063;   % kg-m*     Added Mass Cross-term and Fin Lift *kg-m/rad?

Ypq = -7.117842;    % kg-m      Added Mass Cross-term (-Zqdot)
Ywp =  126.324739;  % kg-m*     Added Mass Cross-term              *kg/rad?
Zvp = -126.324739;  % kg*       Added Mass Cross-term              *kg/rad?
Zrp = -7.117842;    % kg-m*     Added Mass Cross-term (Yrdot)      *kg/rad?
Mpr =  33.463086;   % kg-m2;  // Mpr   , kg-m^2
Mrp =  33.463086;   % kg-m2*    Added Mass Cross-term          *kg-m2/rad2?
Mvp =  7.117842;    % kg-m*     Added Mass Cross-term (-Yrdot)   *kg-m/rad?
Npq = -33.463086;   % kg-m2*    Added Mass Cross-term          *kg-m2/rad2?
Nwp =  7.117842;    % kg-m*   	Added Mass Cross-term (Zqdot)    *kg-m/rad?

%{
% Fin Lift:
dCL     =   4.130000; % n/a     Coef. of Lift Slope
% Fin lift and moment equation coefficients:
Yuudr   =  rho*dCL*Sfin;        % kg/m*    Fin lift force     *kg/m-rad? 
Zuuds   = -rho*dCL*Sfin;        % kg/m*    Fin lift force     *kg/m-rad?
Muuds   =  rho*dCL*Sfin*xfin;   % kg*      Fin lift moment      *kg/rad?
Nuudr   =  rho*dCL*Sfin*xfin;   % kg*      Fin lift moment      *kg/rad?

 [ Yuudr, Zuuds, Muuds, Nuudr ] = LiftLRAUV( rho ); % needs more work

 , Yuv, Yur, Zuq, Zuw, Muw, Muq, Nuv, Nur
 
% Fin Lift:
Sfin    =  115*2/1e4; % m^2     Total area of elevator = 2 x fin.
xfin    =  -0.633;    % m       Midpoint to elevator axle (x)
dCL     =   4.130000; % n/a     Coef. of Lift Slope

% Fin lift coefficients
Yuudr   =   0.5*rho*dCL*Sfin;      % kg/m*    Fin lift force     *kg/m-rad?
Zuuds   =  -0.5*rho*dCL*Sfin;      % kg/m*    Fin lift force     *kg/m-rad?
% Fin moment coefficients:
Muuds   =   0.5*rho*dCL*Sfin*xfin; % kg*      Fin lift moment      *kg/rad?
Nuudr   =   0.5*rho*dCL*Sfin*xfin; % kg*      Fin lift moment      *kg/rad?
%}