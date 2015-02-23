function [bL, bM] = bodyLiftMoment(u, w)

% Calc Body Lift Forces and Moments:
rho     =   1025;           % kg/m3
d       = 0.3048;           % m, Diameter
l       = 2.5;              % m, Length 
xzero   = l/2 + 0.1181;      % m
xcp   = -0.65*l + xzero; 

% Body lift force coefficient (Lbody = -0.5 *rho *(d^2) *Cydbeta *u *w):
Cybeta  = 0.003;                  % Hoerner constant fot 6.7<=(l/d)<=10
Cydbeta = (l/d)*Cybeta*(180/pi);  % Hoerner lift slope coefficient

% Body lift force:
bL = -0.5 * rho * (d^2) * Cydbeta * u * w;
% Body lift moment:
bM  =  bL * xcp;
 
end
