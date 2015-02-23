function [ M, Minv ] = inv_massmatrix( coeffs ) 
% coeffs - string of script finlename containing vehicle coeffs 
%          (just name - leave out .m) 
% vehicle_coeffs ;

run(coeffs); % Load script containing vehicle coeffs

% Form the generalized mass matrix and invert it:
M = zeros(6,6);
m1 = [ m-Xudot,          0,           0,         0,        m*zg,      -m*yg ];
m2 = [       0,    m-Yvdot,           0,     -m*zg,           0, m*xg-Yrdot ];
m3 = [       0,          0,     m-Zwdot,      m*yg, -m*xg-Zqdot,          0 ];
m4 = [       0,      -m*zg,        m*yg, Ixx-Kpdot,           0,          0 ]; 
m5 = [    m*zg,          0, -m*xg-Mwdot,         0,   Iyy-Mqdot,          0 ];
m6 = [   -m*yg, m*xg-Nvdot,           0,         0,           0,  Izz-Nrdot ];

M(:,1:6) = [m1; m2; m3; m4; m5; m6 ];
Minv = inv(M);

end