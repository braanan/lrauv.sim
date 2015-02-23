% ele_offsetLRAUV_SIM.M    -   For LRAUV Vehicle Simulator

% ele_offsetLRAUV_SIM : Minimizes error intruduced by off-set in elevator 
%                       angle data.
% 
% Last modified July 21, 2014
% Ben Raanan


clear 
% close all

startTime = datestr(clock)
tic
h = waitbar(0,'Initializing error minimization...');

fpath = '~/Documents/MATLAB/MBARI/mat/shark/workver/'; % shark/
% filename = [fpath 'LRAUV_SIM_201309282307_201309301141.mat']; % shark data
filename = [fpath 'LRAUV_SIM_201309301141_201310070703.mat']; % shark attack
% filename = [fpath 'LRAUV_SIM_201309121813_201309140344.mat'];  % bottoming

%--------------------------------------------------------------------------
% STATE AND INPUT VECTORS:
% x = [u v w p q r xpos ypos zpos phi theta psi]'
% ui = [ delta_s delta_r Xprop Kprop ]'

[ time, time_step, xstruct, names, controls ] = initialize_LRAUV_SIM( filename );

% time of bottoming 
% (timei>=datenum(2013,09,12,22,13,00) & u >=.8); 

% time of shark attack
% (timei>=datenum(2013,09,30,14,16,44) & timei<=datenum(2013,09,30,14,16,54));

% shark attack dive profile
% (timei>=datenum(2013,09,30,13,04,10) & timei<=datenum(2013,09,30,14,28,24));

timeIn  = datenum(2013,09,30,13,04,10);
timeOut = datenum(2013,09,30,14,28,24);

[~,timeIni] = min(abs(time - timeIn)); 

%--------------------------------------------------------------------------
% initiate first step and set runtime
startPoint = timeIni+240*20;
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

lag = -10:10;


tryVal = (-5:0.025:-4.75);


% hold space                                                
theta_m   = zeros(length(tryVal),n_steps+1);
psi_m     = zeros(length(tryVal),n_steps+1);

waitbar(0.01,h,'Runing error minimization loops...');
for k = 1:length(tryVal)
    
    % Define error range and resolution
    ele_offset  = tryVal(k)*pi/180;         
    rud_offset  = 0*pi/180;         
    
    % Reset control vars   
    ui = controls;
    
    % Correct for offsets in data:
    ui(n,1) = ui(n,1) + ele_offset; % zeros(size(ui(n,1)));
    ui(n,2) = ui(n,2) + rud_offset;
    
    % Correct for Hysteresis offset
    %{
    for q = 1:length(ui)
        if ui(q,1)<0*180/pi
            ui(q,1) = ui(q,1) - db(k)*pi/180;
        elseif ui(q,1)>0*180/pi
            ui(q,1) = ui(q,1) + 0*pi/180;
        end
    end; clear q
    %}
    
    % Unpack state vector
    x = zeros(1,12);
    %
    for c=[1:6,9,10:12];
        x(c) = xstruct.(names{c})(startPoint);
    end; clear c
    %}
    
    % Run simulation
    for i = startPoint:startPoint+n_steps
        
        % Account for movable mass shift 
        xg = Xmass(i);
        
        % Set some vars constant
             x(1) = xstruct.u(i);
        %     x(2) = xstruct.v(i);
        %     x(3) = xstruct.w(i);
             x(4)  = xstruct.p(i);
        %     x(6)  = xstruct.r(i);
             x(10) = xstruct.phi(i);
        %     x(12) = xstruct.psi(i);
        
        
        ui_in = ui(i,:);
       
        [xdot,forces] = lrauv(x,ui_in); % LRAUV % lrauvTest
        
        
        % Log outputed forces
        simlog(i,:) = [x ui_in];
        % f(:,i) = forces;
        
        % EULER INTEGRATION to calculate new states x(n+1)
        % x(i+1) = x(i) + dx/dt*delta_t
        % NOTE: overwriting old states with new states, saving back at the top of the loop
        % x = x + (xdot .* time_step)';
        
        % RUNGE-KUTTA APPROXIMATION to calculate new states
        % NOTE: ideally, should be approximating ui values for k2,k3
        k1_vec = xdot;
        k2_vec = lrauv(x+(0.5.*time_step.*k1_vec)', ((ui(i,:)+ui(i+1,:))./2)) ;
        k3_vec = lrauv(x+(0.5.*time_step.*k2_vec)', ((ui(i,:)+ui(i+1,:))./2)) ;
        k4_vec = lrauv(x+(time_step.*k3_vec)', ui(i+1,:)) ;
        x = x + time_step/6.*(k1_vec +2.*k2_vec +2.*k3_vec +k4_vec)';
        
    end
    
    theta_m(k,:) = simlog(n,11)';
    % psi_m(k,:)   = simlog(n,12)';
    
    for c=1:length(lag)
        lagErrs(c,k) = sum(( simlog(n,11)' - xstruct.theta(n+lag(c)) ).^2);
        [~,minlagi(k)] = min(lagErrs(:,k));
        lagErr(k) = lag(minlagi(k));
    end; clear c lagErrs 
    
    waitbar(k/length(tryVal)-0.01,h,['Runing error minimization loops... ['...
        num2str(100*k/length(tryVal),2) '%]'] );
    pause(0.1);
end

toc
waitbar(1,h,['Error minimization complete [' num2str(100) '%]'] );

% xstruct.theta = zeros(size(xstruct.theta));
% plot
figure;
plot(xstruct.theta(n+lag(median(minlagi)))*(180/pi),'r',...
    'linewidth',2);
hold on;
plot(theta_m' * (180/pi));
xlabel('Time (sec)');
set(gca,'xticklabel',(get(gca,'xtick')./(1/time_step)));

% plot
% figure;
% plot(xstruct.psi(n)*(180/pi),'ro');
% hold on;
% plot(psi_m' * (180/pi));

% Compute squered sum of residuales and find minimizing offset value
res = (theta_m - repmat(xstruct.theta(n+lag(median(minlagi))),k,1)).^2;
err = [ sum(res,2), tryVal' ];
[minError, minErrori] = min(err(:,1));

close(h)

% res = (psi_m - repmat(xstruct.psi(n),k,1)).^2;
% err(:,2) = rud_offset*180/pi;
% err(:,2) = ele_offset*180/pi;

% Rudder_offset   = rud_offset(minErrori)*180/pi;
Elevator_offset = err(minErrori,2)
