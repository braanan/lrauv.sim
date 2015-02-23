% minimizeError_LRAUV_SIM.M    -   For LRAUV Vehicle Simulator

% ele_offsetLRAUV_SIM : Minimizes error intruduced by off-set in elevator
%                       angle data.
%
% Last modified July 21, 2014
% Ben Raanan


clear 
close all
clc ;
startTime = datestr(clock)
tic
% h = waitbar(0,'Initializing error minimization...');


fpath = '~/Documents/MATLAB/MBARI/mat/workver/'; % shark/
% filename = [fpath 'LRAUV_SIM_201309282307_201309301141.mat']; % shark data
filename = [fpath 'LRAUV_SIM_201309121813_201309140344.mat'];  % bottoming


% STATE AND INPUT VECTORS:
% x = [u v w p q r xpos ypos zpos phi theta psi]'
% ui = [ delta_s delta_r Xprop Kprop ]'
[ time, time_step, xstruct, names, controls ] = initialize_LRAUV_SIM( filename );


timeIn  = datenum(2013,09,12,21,00,00);
% timeOut = datenum(2013,09,30,14,28,24);


timeIni = closest(timeIn,time); 


% initiate first step and set runtime
startPoint = timeIni;
timeEval    = 240; % sec, evaluation run time
n_steps     = fix(timeEval/time_step); %size(ui,1); 
n = startPoint:startPoint+n_steps;

% Define global vars
global xg zg Sfin Mqq ARe dCL CDc

zg      =   0.0067940;          % m         Center of gravity             
Sfin    =   1.15e-2;            % m^2       Fin area
Mqq     =   0.35*-632.698957;   % kg-m2     Cross-flow drag (Mq|q|) 
ARe     =   6.500000;           % n/a     Fin aspect ratio
dCL     =   1.5*4.130000;       % n/a     Coef. of Lift Slope
CDc     =   0.030000;           % n/a     Crossflow Coef. of Drag !!try 0.6!!

mass        =   147.8671;         % kg Flooded Vehicle total mass
movableMass =   26;               % kg Battary movable mass
dropWtMass  =   1.0;              % kg Mass of the drop weight #1, kg

dropWtX     =  -0.1330;           % m  X location of the drop weight #1, m
Xmass = (movableMass.*xstruct.mass_p + dropWtMass*dropWtX)./mass;


% Define error range and resolution
ele_offset  = 0*pi/180; 
rud_offset  = 0*pi/180; 

tryVal1 = (-0.46:0.02:-0.16);
tryVal2 = (-1:0.1:0.8);

nTry=length(tryVal1)*length(tryVal2);


% Hold space
err   = NaN(length(tryVal2),length(tryVal1));
theta_m = NaN(n_steps,length(tryVal2));
% waitbar(0.01,h,'Runing error minimization loops...');

cloop=1; 
for b = 1:length(tryVal1)
   
    
    dropWtX = tryVal1(b);
    Xmass   = (movableMass.*xstruct.mass_p + dropWtMass*dropWtX)./mass;

    
    for k = 1:length(tryVal2)
        
        
        % Reset control vars
        ui = controls(n,:);
        
        % Correct for offsets in data:
        ui(:,1) = ui(:,1) + ele_offset;
        ui(:,2) = ui(:,2) + rud_offset;
        
        % Correct for hysteresis and backlash offsets
        %
        ui(ui(:,1)<0*180/pi,1) = ui(ui(:,1)<0*180/pi,1) - tryVal2(k)*pi/180;
        ui(ui(:,1)>0*180/pi,1) = ui(ui(:,1)>0*180/pi,1) + 0*pi/180;
        %}
        
        % Unpack state vector
        x = zeros(1,12);
        for c=[1:6,9,10:12];
            x(c) = xstruct.(names{c})(startPoint);
        end; 
        
        
        simlog = zeros(n_steps,16);
        % Run simulation
        for i = 1:n_steps
            
            % Account for movable mass shift
            xg = Xmass(i);
            
            % Set some vars constant
            x(1) = xstruct.u(n(i));
            %     x(2) = xstruct.v(n(i));
            %     x(3) = xstruct.w(n(i));
            x(4)  = xstruct.p(n(i));
            %     x(6)  = xstruct.r(n(i));
            x(10) = xstruct.phi(n(i));
            %     x(12) = xstruct.psi(n(i));
            
            ui_in = ui(i,:);
            
            % Calc next step
            [xdot,~] = lrauv(x,ui_in); % LRAUV % lrauvTest
            
            % Log step data
            simlog(i,:) = [x ui_in];
            
            % RUNGE-KUTTA APPROXIMATION to calculate new states
            % NOTE: ideally, should be approximating ui values for k2,k3
            k1_vec = xdot;
            k2_vec = lrauv(x+(0.5.*time_step.*k1_vec)', ((ui(i,:)+ui(i+1,:))./2)) ;
            k3_vec = lrauv(x+(0.5.*time_step.*k2_vec)', ((ui(i,:)+ui(i+1,:))./2)) ;
            k4_vec = lrauv(x+(time_step.*k3_vec)', ui(i+1,:)) ;
            x = x + time_step/6.*(k1_vec +2.*k2_vec +2.*k3_vec +k4_vec)';
            
        end; 
        
        
        % Keep log of modeled theta vectors
        theta_m(:,k) = simlog(:,11);
        
%             cloop=cloop+1;
%             prg=cloop/nTry;
%             waitbar(prg,h,['Runing error minimization loop ' num2str(b) ' ['...
%             num2str(100*prg,2) '%]'] );
        
    end
    
    % Compute squered sum of residuales and find minimizing offset value
    res = (theta_m - repmat(xstruct.theta(n(1):n(end-1)),k,1)').^2;
    err(:,b) = sum(res);
    
    
    display([ datestr(clock) ': Iteration cycle ' num2str(b) ' of ' num2str(length(tryVal1)) ' complete.']);
end

%{
% Compute squered sum of residuales and find minimizing offset value
    res2 = (theta_m2 - repmat(xstruct.theta(n),length(trySfin),1)).^2;
    err2 = sum(res2,2);
    err2(:,2) = trySfin;
    [minError2, minError2i] = min(err2(:,1));

    % Log offset results
    offset2(a,:)  = [offset(minError2i,:), minError2, zg];
    theta_m3(a,:) = theta_m2(minError2i,:);
    display([ datestr(clock) ': Zg iteration number ' num2str(a) ' complete.']);
    pause(0.1)
end;
%}

toc
% waitbar(1,h,['Error minimization complete [' num2str(100) '%]'] );

% plot results
%{
figure;
surf(tryVal1,tryVal2,err)
title('LRAUV simulator error minimization function',...
    'fontweight','bold','fontsize',22)

xlabel('Deviation from static X center of gravity (cm)',...
    'fontweight','bold','fontsize',16)
set(gca,'xticklabel',100*dropWtMass*tryVal1/mass)
set(get(gca,'xlabel'),'rotation',-1);
ylabel('Elevator offset correction (deg)',...
    'fontweight','bold','fontsize',16)
set(get(gca,'ylabel'),'rotation',74);
zlabel('\Sigma Error',...
    'fontweight','bold','fontsize',16)
box off; grid on;
%}


%{
figure;
plot(xstruct.theta(n)*(180/pi),'ro');
hold on;
plot(theta_m' * (180/pi));
%}