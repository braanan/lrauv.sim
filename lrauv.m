% lrauv.m Vehicle Simulator Testground
% Returns the time derivative of the state vector 
% Last modified July 17, 2014

function [ACCELERATIONS,FORCES] = lrauv(x,ui)

% TERMS
% --------------------------------------------------------------------- 
% STATE VECTOR:
% x = [u v w p q r xpos ypos zpos phi theta psi]'
%  Body-referenced Coordinates
%  u            = Surge velocity            [m/sec]
%  v            = Sway velocity             [m/sec]
%  w            = Heave velocity            [m/sec]    
%  p            = Roll rate                 [rad/sec]
%  q            = Pitch rate                [rad/sec]
%  r            = Yaw rate                  [rad/sec]
%  Earth-fixed coordinates
%  xpos         = Position in x-direction   [m]
%  ypos         = Position in y-direction   [m]
%  zpos         = Position in z-direction   [m]
%  phi          = Roll angle                [rad]
%  theta        = Pitch angle               [rad]
%  psi          = Yaw angle                 [rad]
%
% INPUT VECTOR
% ui = [delta_s delta_r]'
%  Control Fin Angles
%  delta_s  = angle of stern planes         [rad]
%  delta_r  = angle of rudder planes        [rad]

% Initialize global variables 
%--------------------------------------------------------------------- 
% load vdata          ;  % W and B, CG and CB coords
% [ ~, Minv ] = inv_massmatrix( 'vehicle_coeffs' );   % Minv matrix
vehicle_coeffs ;             % non-zero vehicle coefficients only

%
% Form the generalized mass matrix and invert it:
M = zeros(6,6);
M(1,:) = [ m-Xudot,          0,           0,         0,        m*zg,      -m*yg ];
M(2,:) = [       0,    m-Yvdot,           0,     -m*zg,           0, m*xg-Yrdot ];
M(3,:) = [       0,          0,     m-Zwdot,      m*yg, -m*xg-Zqdot,          0 ];
M(4,:) = [       0,      -m*zg,        m*yg, Ixx-Kpdot,           0,          0 ]; 
M(5,:) = [    m*zg,          0, -m*xg-Mwdot,         0,   Iyy-Mqdot,          0 ];
M(6,:) = [   -m*yg, m*xg-Nvdot,           0,         0,           0,  Izz-Nrdot ];
Minv = inv(M);
%}

% Output flags
% show_forces = 10 ;

% Get and check state variables and control inputs
%---------------------------------------------------------------------
% Get state variables
u   = x(1) ; v  = x(2) ; w  = x(3) ; p  = x(4) ; q  = x(5) ; r  = x(6);
phi = x(10) ; theta  = x(11) ; psi  = x(12) ;

% Get control inputs
% delta_s = ui(1) ; delta_r = ui(2) ; 
Xprop = ui(3) ; Xuu = ui(4) ; 

% Check control inputs (useful later)
% delta_max = 15*pi/180;
% if abs(delta_s) > delta_max
%     delta_s = sign(delta_s)*delta_max ;
% end;
% if abs(delta_r) > delta_max
%     delta_r = sign(delta_r)*delta_max ;
% end

% Initialize elements of coordinate system transform matrix
%---------------------------------------------------------------------
c1 = cos(phi); c2 = cos(theta); c3 = cos(psi); 
s1 = sin(phi); s2 = sin(theta); s3 = sin(psi); 
t2 = tan(theta);

% Get fin forces and moments
%---------------------------------------------------------------------
[ F1, F2, F3, F4, M1, M2, M3, M4 ] = robsFins( ui, x );

% [bodyLift, bodyMoment] = bodyLiftMoment(u, w);
% [ Yr, Zs, Ms, Nr ] = fins( rho, u, v, w, q, r, delta_r, delta_s );
% Set total forces from equations of motion
%---------------------------------------------------------------------

% X = 0;
% X = -(W-B)*sin(theta) + Xuu*u*abs(u) + (Xwq-m)*w*q + (Xqq + m*xg)*q^2 ...
%     + (Xvr+m)*v*r + (Xrr + m*xg)*r^2 - m*yg*p*q - m*zg*p*r ... 
%     + Xprop + F1(1) + F2(1) + F3(1) + F4(1);

X = Xprop - (Wp-Bp)*s2 + Xuu*u*abs(u)...
    + m*(v*r - w*q + Xgp*(q*q + r*r) - Ygp*p*q - Zgp*p*r) ...
    + F1(1) + F2(1) + F3(1) + F4(1)...
    + Xvv*v*v + Xww*w*w + Xvr*v*r + Xwq*w*q + Xrr*r*r + Xqq*q*q;

% Y = 0;
% Y = (W-B)*cos(theta)*sin(phi) + Yvv*v*abs(v) + Yrr*r*abs(r) + Yuv*u*v ...
%   + (Ywp+m)*w*p + (Yur-m)*u*r + (m*zg)*q*r + (Ypq - m*xg)*p*q ... 
%   + F1(2) + F2(2) + F3(2) + F4(2) ; % Yr; % + F1(2) + F2(2) + F3(2) + F4(2); % + Yuudr*u^2*delta_r ;

Y = (Wp-Bp)*c2*s3 + Yvv*v*abs(v)...
    + Yuv*u*v + Yur*u*r + Yrr*r*abs(r) + Ywp*w*p...
    + m*(w*p - u*r + Ygp*(r*r+p*p) -Zgp*q*r - Xgp*p*q)...
    + F1(2) + F2(2) + F3(2) + F4(2);

% % Z = 0;
% Z = (W-B)*cos(theta)*cos(phi) + Zww*w*abs(w) + Zqq*q*abs(q) + Zuw*u*w ...
%     + (Zuq+m)*u*q + (Zvp-m)*v*p + (m*zg)*p^2 + (m*zg)*q^2 ... 
%     + (Zrp - m*xg)*r*p + F1(3) + F2(3) + F3(3) + F4(3) ; % Zs; % F1(3) + F2(3) + F3(3) + F4(3); %+ Zuuds*u^2*delta_s ; 

Z = (Wp-Bp)*c2*c3 + Zww*w*abs(w) + Zqq*q*abs(q)...
    + Zuq*u*q + Zuw*u*w + Zvp*v*p...
    + m*(u*q - v*p + Zgp*(p*p + q*q) - Xgp*p*r - yg*q*r)...
    + F1(3) + F2(3) + F3(3) + F4(3) ; % + bodyLift ;

%  K  = 0;
K = -(yg*W-yb*B)*cos(theta)*cos(phi) - (zg*W-zb*B)*cos(theta)*sin(phi) ...
    + Kpp*p*abs(p) - (Izz-Iyy)*q*r - (m*zg)*w*p + (m*zg)*u*r...
    + Kprop + M1(1) + M2(1) + M3(1) + M4(1);

% K = -(Zgp*Wp - Zbp*Bp)*c2*s1 + (Ygp*Wp - Ybp*Bp)*c2*c1...    % ROLL MOMENTS */
%     + Kpp*p*abs(p) - (Izz - Iyy)*q*r...
%     - m*(Ygp*(v*p - u*q) + Zgp*(w*p - u*r))...
%     + M1(1) + M2(1) + M3(1) + M4(1)...
%     + Kprop ;

% M = 0;
% M = -(zg*W-zb*B)*sin(theta) - (xg*W-xb*B)*cos(theta)*cos(phi) + Mww*w*abs(w) ...
%     + Mqq*q*abs(q) + (Mrp - (Ixx-Izz))*r*p + (m*zg)*v*r - (m*zg)*w*q ... 
%     + (Muq - m*xg)*u*q + Muw*u*w + (Mvp + m*xg)*v*p ...
%     + M1(2) + M2(2) + M3(2) + M4(2) ; %Ms; % M1(2) + M2(2) + M3(2) + M4(2) ;%     + Muuds*u^2*delta_s ;

M = -(Zgp*Wp - Zbp*Bp)*s2 - (Xgp*Wp - Xbp*Bp)*c2*c1...       % PITCH MOMENTS */
    + Mww*w*abs(w) + Mqq*q*abs(q) ...
    + Muw*u*w + Muq*u*q + Mpr*p*r...
    + (Izz - Ixx)*p*r...
    - m*(Zgp*(w*q - v*r) + Xgp*(u*q - v*p))...
    + M1(2) + M2(2) + M3(2) + M4(2) ; % + bodyMoment

% N = 0;
% N = -(xg*W-xb*B)*cos(theta)*sin(phi) - (yg*W-yb*B)*sin(theta) ...
%     + Nvv*v*abs(v) + Nrr*r*abs(r) + Nuv*u*v ...
%     + (Npq - (Iyy-Ixx))*p*q + (Nwp - m*xg)*w*p + (Nur + m*xg)*u*r ... 
%     + M1(3) + M2(3) + M3(3) + M4(3) ; %Nr; % M1(3) + M2(3) + M3(3) + M4(3);% + Nuudr*u^2*delta_r ;

N = (Ygp*Wp - Ybp*Bp)*s2 + (Xgp*Wp - Xbp*Bp)*c2*s1...         % YAW MOMENTS */
    + Nvv*v*abs(v) + Nrr*r*abs(r) + Nuv*u*v ...
    + Nur*u*r + Npq*p*q...
    + (Ixx - Iyy)*p*q...
    - m*Xgp*(u*r - w*p) + m*Ygp*(w*q - v*r)...
    + M1(3) + M2(3) + M3(3) + M4(3);

FORCES = [X Y Z K M N]' ;



ACCELERATIONS = ...
  [Minv(1,1)*X+Minv(1,2)*Y+Minv(1,3)*Z+Minv(1,4)*K+Minv(1,5)*M+Minv(1,6)*N
   Minv(2,1)*X+Minv(2,2)*Y+Minv(2,3)*Z+Minv(2,4)*K+Minv(2,5)*M+Minv(2,6)*N
   Minv(3,1)*X+Minv(3,2)*Y+Minv(3,3)*Z+Minv(3,4)*K+Minv(3,5)*M+Minv(3,6)*N
   Minv(4,1)*X+Minv(4,2)*Y+Minv(4,3)*Z+Minv(4,4)*K+Minv(4,5)*M+Minv(4,6)*N
   Minv(5,1)*X+Minv(5,2)*Y+Minv(5,3)*Z+Minv(5,4)*K+Minv(5,5)*M+Minv(5,6)*N
   Minv(6,1)*X+Minv(6,2)*Y+Minv(6,3)*Z+Minv(6,4)*K+Minv(6,5)*M+Minv(6,6)*N
   c3*c2*u + (c3*s2*s1-s3*c1)*v + (s3*s1+c3*c1*s2)*w
   s3*c2*u + (c1*c3+s1*s2*s3)*v + (c1*s2*s3-c3*s1)*w
     -s2*u +            c2*s1*v +            c1*c2*w 
         p +            s1*t2*q +            c1*t2*r
                           c1*q -               s1*r
                        s1/c2*q +            c1/c2*r] ;
                    
%                     pause
end