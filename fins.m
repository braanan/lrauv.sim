function [ Yr, Zs, Ms, Nr ] = fins( rho, u, v, w, q, r, delta_r, delta_s )

% Last modified July 16, 2014

% Fin Lift:
Sfin    =   2*1.15e-2;      % m^2     Total area of elevator = 2 x fin.
% Sfin    =   1.15e-2;      % m^2     Area of elevator
bfin    =   18.57e-2;       % m       Fin span 
zfin    =   0.152;          % m       Centerline to fin  
xfin    =  -0.633;          % m       Midpoint to elevator axle (x)
dCL     =   4.130000;       % n/a     Coef. of Lift Slope

% ahat    =   0.9;                        % Empirical factor (Hoerner)
% ARe     =   2*((bfin^2)/Sfin);          % Effective fin aspect retio
% cLa = 1/((1/2*ahat*pi) + (1/(pi*ARe))); % cL_alpha (Hoerner [pg. 3-2])

cLa     =   dCL;

% Forces:
Yr =  0.5*rho*cLa*Sfin*((u^2)*delta_r - u*v - xfin*u*r);
Zs = -0.5*rho*cLa*Sfin*((u^2)*delta_s - u*w - xfin*u*q);

% Moments:
Ms =  0.5*rho*cLa*Sfin*xfin*((u^2)*delta_s - u*w - xfin*u*q);
Nr =  0.5*rho*cLa*Sfin*xfin*((u^2)*delta_r - u*v - xfin*u*r);

% Fin lift and moment equation coefficients:
% Yuudr   =  rho*cLa*Sfin;        % kg/m*    Fin lift force     *kg/m-rad? 
% Zuuds   = -rho*cLa*Sfin;        % kg/m*    Fin lift force     *kg/m-rad?
% Muuds   =  rho*cLa*Sfin*xfin;   % kg*      Fin lift moment      *kg/rad?
% Nuudr   =  rho*cLa*Sfin*xfin;   % kg*      Fin lift moment      *kg/rad?
end
