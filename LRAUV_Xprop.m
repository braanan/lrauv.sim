function [ Xprop, Xuabu ] = LRAUV_Xprop( velocity, rho )

% This is a simple model for the LRAUV propulsion system, here we assume 
% that the propeller thrust matches the vehicle axial drag. Function 
% caluclates thrust and drag.

 u      =   velocity;       % m/s
 dia    =   0.3048;         % m     Hull frontal diameter
 Af     =   pi*(dia/2)^2;   % m2    Hull frontal area
 AreaE  =   115*2/1e4;      % m2    Total area of elevator = 2 x fin.
 Cd     =   0.185;          % n/a   Veh drag in x, excluding control surfaces
 Cd0    =   0.03;           % n/a   Control surface drag
 
 Xuabu   =   -.5*rho*Cd*(Af + AreaE*Cd0);
 Xprop   =   -Xuabu.*(u.^2);
 
end