function [ time, time_step, xstruct, names, ui ] = initialize_LRAUV_SIM( filename )
% ---------------------------------------------------------------------
% initialize_LRAUV_SIM: Arrange LRAUV data passed from
%                       processLargeMAT_CRITICAL.m in state vecotr format.
%
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

% Load file
%---------------------------------------------------------------------
load( filename ) % {:}
list = fieldnames( interpVars );


% Compute first derivitves:
%---------------------------------------------------------------------
time_step   = median(diff(time))*(24*3600); % time step = 1sec
time_step   = ceil(time_step*1000)/1000;    % round

if any(~cellfun('isempty',(regexpi(list,'roll'))))
    drolldt     = diff(interpVars.platform_roll_angle)./time_step;
end
if any(~cellfun('isempty',(regexpi(list,'pitch_angle'))))
    dpitchdt    = diff(interpVars.platform_pitch_angle)./time_step;
end
if any(~cellfun('isempty',(regexpi(list,'orientation'))))
    dyawdt      = diff(interpVars.platform_orientation)./time_step;
end
if any(~cellfun('isempty',(regexpi(list,'depth'))))
    dz          = diff(interpVars.depth)./time_step;
end



% Fin motion derivatives
%---------------------------------------------------------------------
% dele = diff(ui(n(1):n(end)+1,1)*180/pi);
% drud = diff(ui(n(1):n(end)+1,2)*180/pi);



% State vector:
%---------------------------------------------------------------------
try
    xstruct.u       = interpVars.platform_speed_wrt_sea_water;
catch
    try
        xstruct.u       = interpVars.platform_x_velocity_wrt_sea_water;
    catch
        xstruct.u       = interpVars.Cmd.speedCmd;
    end
end

try
    xstruct.v       = interpVars.platform_y_velocity_wrt_sea_water;
catch
    xstruct.v       = zeros(size(time));
end

try
    xstruct.w       = interpVars.platform_w_velocity_wrt_sea_water;
catch
    xstruct.w       = zeros(size(time));
end

xstruct.p       = [  drolldt(1),  drolldt ]; %[       dz(1),       dz ];
xstruct.q       = [ dpitchdt(1), dpitchdt ];
xstruct.r       = [   dyawdt(1),   dyawdt ];
xstruct.x       = zeros(size(time));
xstruct.y       = zeros(size(time));
xstruct.z       = interpVars.depth;
xstruct.phi     = interpVars.platform_roll_angle;
xstruct.theta   = interpVars.platform_pitch_angle;
xstruct.psi     = interpVars.platform_orientation;
xstruct.z_rate_vh  = interpVars.depth_rate;
xstruct.z_rate  = [      dz(1),       dz ];


try
    xstruct.mass_p  = interpVars.platform_mass_position;
catch
    xstruct.mass_p  = interpVars.Cmd.massPositionCmd;
end

try
    xstruct.buoyancy_p = interpVars.platform_buoyancy_position;
catch
    try
        xstruct.buoyancy_p = interpVars.Cmd.buoyancyCmd;
    catch
        xstruct.buoyancy_p = NaN(size(time));
    end
end

% var names
xname = {'u','v','w','p','q','r','x','y','z','phi','theta','psi'};

% Commanded variables
%---------------------------------------------------------------------
xstruct.Cmd = interpVars.Cmd;

% Temp and Sal
try
    xstruct.Temp = interpVars.sea_water_temperature;
catch
    xstruct.Temp = NaN(size(time));
end

try
    xstruct.Sal = interpVars.sea_water_salinity;
catch
    xstruct.Sal = NaN(size(time));
end
% fill missing Cmd vars with NaN array
%---------------------------------------------------------------
cmdf = {'elevatorAngleCmd','massPositionCmd','depthCmd','depthRateCmd',...
    'pitchCmd','pitchRateCmd','buoyancyCmd','speedCmd','verticalMode'};
listCmd = fieldnames(interpVars.Cmd);
for k=1:numel(cmdf)
    checkCmd = strcmp(listCmd,cmdf{k});
    if ~any(checkCmd)
        xstruct.Cmd.(cmdf{k})=NaN(size(time));
    end
end



% Compute Xprop:
%---------------------------------------------------------------------
rho          = 1025; % kg/m3
velocity     = xstruct.u;
[ ustruct.Xprop, ustruct.Xuabu] = LRAUV_Xprop( velocity, rho );



% Controls vector:
%---------------------------------------------------------------------
ustruct.delta_s      = interpVars.platform_elevator_angle;
ustruct.delta_r      = interpVars.platform_rudder_angle;


uiname = {'delta_s','delta_r','Xprop','Xuabu'};
ui=zeros(length(ustruct.delta_s),4);
for c=1:4
    ui(:,c) = ustruct.(uiname{c})';
end

names = [ xname, uiname ];
end