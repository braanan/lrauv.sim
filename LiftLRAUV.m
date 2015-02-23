function [ Yuudr, Zuuds, Muuds, Nuudr, Yuv, Yur, Zuq, Zuw, Muw, Muq, Nuv, Nur ] = LiftLRAUV( rho )

% % Body Lift Forces and Moments:
% d    = 0.3048;       % m, Diameter
% l    = 2;          % m, Length 
% xzero = l/2 + 0.1181;      % m
% 
% % Body lift force coefficient (Lbody = -0.5 *rho *(d^2) *Cydbeta *u *w):
% Cybeta  = 0.003;                  % Hoerner constant fot 6.7<=(l/d)<=10
% Cydbeta = (l/d)*Cybeta*(180/pi);  % Hoerner lift slope coefficient
% 
% Yuvl    = -0.5 *rho *(d^2) *Cydbeta;
% Zuwl    = -0.5 *rho *(d^2) *Cydbeta;
% 
% % Body lift moment coefficients:
% xcp   = -0.65*l + xzero;
% Muwl  = -0.5 *rho *(d^2) *Cydbeta *xcp;
% Nuvl  = -Muwl;

% Fin Lift:
Sfin    =   2*1.15e-2; % m^2     Area of elevator
xfin    =  -0.633;      % m       Midpoint to elevator axle (x)
dCL     =   4.130000;   % n/a     Coef. of Lift Slope

Zwdot   =  -126.324739;     % kg          
Zqdot   =   7.117842;       % kg-m/rad    
Xudot   =  -4.876161;       % kg           
Yvdot   =   Zwdot;          % kg          Terms should be the same - needs checking..
Yrdot   =  -Zqdot;          % kg-m        Terms should be the same - needs checking..

% Fin lift and moment equation coefficients:
Yuudr   =  0.5*rho*dCL*Sfin;        % kg/m*    Fin lift force     *kg/m-rad? 
Zuuds   = -0.5*rho*dCL*Sfin;        % kg/m*    Fin lift force     *kg/m-rad?
Muuds   =  0.5*rho*dCL*Sfin*xfin;   % kg*      Fin lift moment      *kg/rad?
Nuudr   =  0.5*rho*dCL*Sfin*xfin;   % kg*      Fin lift moment      *kg/rad?

% Fin lift and moment coefficients:
Yuvf    =  -Yuudr;
Yurf    =   Zuuds *xfin;
Zuwf    =   Zuuds;
Zuqf    =  -Yurf; 

Muwf    =   Muuds;
Muqf    =  -Muuds *xfin; 
Nuvf    =  -Nuudr;
Nurf    =   Muqf;

% Combined terms:
Yura    =   Xudot;           % Kg   Added mass cross-term
Zuqa    =  -Xudot;           % Kg   Added mass cross-term
Muwa    =  -(Zwdot - Xudot); % Kg   Added mass cross-term
Muqa    =  -Zqdot;           % kg-m/rad   Added mass cross-term
Nuva    =  -(Xudot - Yvdot); % Kg   Added mass cross-term
Nura    =   Yrdot;           % kg-m Added mass cross-term  

Yuv     =   Yuvf;% + Yuvl;
Yur     =   Yura + Yurf;
Zuw     =   Zuwf; % + Zuwl;
Zuq     =   Zuqa + Zuqf;
Muw     =   Muwa + Muwf; % + Muwl;
Muq     =   Muqa + Muqf;
Nuv     =   Nuva + Nuvf; % + Nuvl;
Nur     =   Nura + Nurf;

end
