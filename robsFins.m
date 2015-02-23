function [ F1, F2, F3, F4, M1, M2, M3, M4 ] = robsFins( ui_in, x )

% robsFins: Simulate lift/drag forces and moments applied to each of LRAUV 
%           fins.

% Initialize global variables
%---------------------------------------------------------------------
% vehicle_coeffs ;

global Sfin ARe dCL CDc

rho     =   1025;         % kg/m3
% ARe     =   6.500000;     % n/a     Fin aspect ratio
% dCL     =   4.130000;     % n/a     Coef. of Lift Slope
% CDc     =   0.600000;     % n/a     Crossflow Coef. of Drag !!try 0.6!!
Cd0     =   0.030000;     % n/a     Min reaction drag coeff
ec      =   0.9;
% Sfin    =   1.15e-2;          % m^2   Fin area
% bfin    =   18.57e-2;         % m     Fin span
% ARe     = 2*((bfin^2)/Sfin);  % n/a   Effective aspect ratio

stall_angle    =   30.000000*pi/180; % arcdeg

% Position of fins from vehicle center:
% colName = {'lowerRud','portElev','upperRud','stbdElev'};
% rowName = {'X','Y','Z'};
xi = 1 ; yi = 2 ; zi = 3 ; % dimention indecies
finPos  = [ -0.633, -0.633, -0.633, -0.633 ;
             0.012, -0.152,  0.012,  0.152 ;
            -0.152,  0.000,  0.152,  0.000 ] ;  % m

% Get and check state variables and control inputs
%---------------------------------------------------------------------
% Get state variables
u   = x(1) ; v  = x(2) ; w  = x(3) ; p  = x(4) ; q  = x(5) ; r  = x(6);

% Get control inputs
elev_ang = ui_in(1) ; rud_ang = ui_in(2) ;

% Initialize elements of coordinate system transform matrix
%---------------------------------------------------------------------
s1 = sin(rud_ang); s2 = sin(elev_ang);
c1 = cos(rud_ang); c2 = cos(elev_ang);


% Calculate fin velocities in body coords
%---------------------------------------------------------------------
v1  = [ u + q*finPos(zi,1) - r*finPos(yi,1) ;
        v - p*finPos(zi,1) + r*finPos(xi,1) ;
        w + p*finPos(yi,1) - q*finPos(xi,1) ] ;

v2  = [ u + q*finPos(zi,2) - r*finPos(yi,2) ;
        v - p*finPos(zi,2) + r*finPos(xi,2) ;
        w + p*finPos(yi,2) - q*finPos(xi,2) ] ;

v3 	= [ u + q*finPos(zi,3) - r*finPos(yi,3) ;
        v + r*finPos(xi,3) - p*finPos(zi,3) ;
        w + p*finPos(yi,3) - q*finPos(xi,3) ] ;

v4  = [ u + q*finPos(zi,4) - r*finPos(yi,4) ;
        v - p*finPos(zi,4) + r*finPos(xi,4) ;
        w + p*finPos(yi,4) - q*finPos(xi,4) ] ;


% Now get angle of attack for each fin
%---------------------------------------------------------------------
norm_v1 = sqrt( v1(1)*v1(1) + v1(2)*v1(2) + v1(3)*v1(3) );
if (norm_v1 < 0.001)
    alpha1 = 0.0;
else
    alpha1 =  rud_ang - v1(2)/v1(1); % asin(-v1(2)/norm_v1) + rud_ang;
end

norm_v2 = sqrt( v2(1)*v2(1) + v2(2)*v2(2) + v2(3)*v2(3) );
if (norm_v2 < 0.001)
    alpha2 = 0.0;
else
    alpha2 =  elev_ang + v2(3)/v2(1); % atan2(v2(3),v2(1)) + elev_ang; 
end

norm_v3 = sqrt( v3(1)*v3(1) + v3(2)*v3(2) + v3(3)*v3(3) );
if (norm_v3 < 0.001)
    alpha3=0.0;
else
    alpha3 = rud_ang - v3(2)/v3(1); %asin(-v3(2)/norm_v3) + rud_ang; % 
end

norm_v4 = sqrt( v4(1)*v4(1) + v4(2)*v4(2) + v4(3)*v4(3) );
if (norm_v4 < 0.001)
    alpha4=0.0;
else
    alpha4 = elev_ang + v4(3)/v4(1); % atan2(v4(3),v4(1)) + elev_ang;
end

% lift and drag coefficients */
CDC = CDc/ARe;

% Note that if stall angle is exceeded: NO LIFT */
if (abs(alpha1) < stall_angle)
    CL1 = dCL*alpha1 + CDC*alpha1*abs(alpha1);
else CL1 = 0. ;
end

if (abs(alpha2) < stall_angle)
    CL2 = dCL*alpha2 + CDC*alpha2*abs(alpha2);
else CL2 = 0. ;
end

if (abs(alpha3) < stall_angle)
    CL3 = dCL*alpha3 + CDC*alpha3*abs(alpha3);
else CL3 = 0. ;
end

if (abs(alpha4) < stall_angle)
    CL4 = dCL*alpha4 + CDC*alpha4*abs(alpha4);
else CL4 = 0. ;
end

aa = 1.0/(pi*ARe*ec);

CD1 = Cd0 + aa*CL1*CL1;
CD2 = Cd0 + aa*CL2*CL2;
CD3 = Cd0 + aa*CL3*CL3;
CD4 = Cd0 + aa*CL4*CL4;

% lift and drag forces, in flow coords... */
%---------------------------------------------------------------------
cons = (rho*Sfin)/2.0;

LW1 = cons*norm_v1*norm_v1*CL1;      % positive when the lift vector is close to normal vector */
LW2 = cons*norm_v2*norm_v2*CL2;
LW3 = cons*norm_v3*norm_v3*CL3;
LW4 = cons*norm_v4*norm_v4*CL4;

DW1 = cons*norm_v1*norm_v1*CD1;      % always positive */
DW2 = cons*norm_v2*norm_v2*CD2;
DW3 = cons*norm_v3*norm_v3*CD3;
DW4 = cons*norm_v4*norm_v4*CD4;

LF1 = LW1*cos(alpha1) + DW1*sin(alpha1);        % force in the fin normal direction */
LF2 = LW2*cos(alpha2) + DW2*sin(alpha2);
LF3 = LW3*cos(alpha3) + DW3*sin(alpha3);
LF4 = LW4*cos(alpha4) + DW4*sin(alpha4);

DF1 = -LW1*sin(alpha1) + DW1*cos(alpha1);       % force in the fin-aft direction */
DF2 = -LW2*sin(alpha2) + DW2*cos(alpha2);
DF3 = -LW3*sin(alpha3) + DW3*cos(alpha3);
DF4 = -LW4*sin(alpha4) + DW4*cos(alpha4);

% Finally, transform into the body frame */
%---------------------------------------------------------------------
F1  = [ -LF1*s1 + (-DF1)*c1 ;        
         LF1*c1 + (-DF1)*s1 ;
         0.0                ] ;

F2  = [ -LF2*s2 + (-DF2)*c2 ;
         0.0                ;
        -LF2*c2 - (-DF2)*s2 ] ;

F3	= [ -LF3*s1 + (-DF3)*c1 ;
         LF3*c1 + (-DF3)*s1 ;
         0.0                ] ;

F4  = [ -LF4*s2 + (-DF4)*c2 ;
         0.0                ;
        -LF4*c2 - (-DF4)*s2 ] ;

% moments induced by these forces */  
M1  = [ finPos(yi,1)*F1(3) - finPos(zi,1)*F1(2) ;   
       -finPos(xi,1)*F1(3) + finPos(zi,1)*F1(1) ;
        finPos(xi,1)*F1(2) - finPos(yi,1)*F1(1) ] ;

M2  = [  finPos(yi,2)*F2(3) - finPos(zi,2)*F2(2) ;
        -finPos(xi,2)*F2(3) + finPos(zi,2)*F2(1) ;
         finPos(xi,2)*F2(2) - finPos(yi,2)*F2(1) ] ;

M3  = [  finPos(yi,3)*F3(3) - finPos(zi,3)*F3(2) ;
        -finPos(xi,3)*F3(3) + finPos(zi,3)*F3(1) ;
         finPos(xi,3)*F3(2) - finPos(yi,3)*F3(1) ] ;

M4  = [  finPos(yi,4)*F4(3) - finPos(zi,4)*F4(2) ;
        -finPos(xi,4)*F4(3) + finPos(zi,4)*F4(1) ;
         finPos(xi,4)*F4(2) - finPos(yi,4)*F4(1) ] ;
     end