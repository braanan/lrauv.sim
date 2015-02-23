% LRAUV_SIM.M
% LRAUV_SIM: Main script for runing vehicle simulation.
% Last modified Aug 1, 2014
% Ben Raanan


clear
close all

h = waitbar(0,'Initializing LRAUV Vehicle Simulator...');

fpath = '~/Documents/MATLAB/MBARI/mat/workver/'; % shark/
% filename = [fpath 'LRAUV_SIM_201309282307_201309301141.mat']; % shark data
filename = [fpath 'LRAUV_SIM_hf_201309121813_201309140344.mat'];  % bottoming

%--------------------------------------------------------------------------
% STATE AND INPUT VECTORS:
% x = [u v w p q r xpos ypos zpos phi theta psi]'
% ui = [ delta_s delta_r Xprop Kprop ]'

[ time, time_step, xstruct, names, controls ] = initialize_LRAUV_SIM( filename );
%{
% time of bottoming
% time >= datenum(2013,09,12,22,13,00)

% time of shark attack
% timei > =datenum(2013,09,30,14,16,44) & time <= datenum(2013,09,30,14,16,54)

% shark attack dive profile
% time >= datenum(2013,09,30,13,04,10) & timei <= datenum(2013,09,30,14,28,24)
%}
timeIn  = datenum(2013,09,12,20,21,12);
timeOut = datenum(2013,09,30,14,28,24);

[~,timeIni] = min(abs(time - timeIn));



% Initiate first step and set runtime
%--------------------------------------------------------------------------
startPoint = timeIni+240*20;
timeEval    = 240; % sec, evaluation run time
n_steps     = fix(timeEval/time_step);
n = startPoint:startPoint+n_steps;
n_ind = 1:length(n);



% Define global vars
%--------------------------------------------------------------------------
global xg zg Sfin Mqq ARe dCL CDc

zg      =   0.0067940;          % m       Center of gravity
Sfin    =   1.15e-2;            % m^2     Fin area
Mqq     =   0.35*-632.698957;   % kg-m2   Scaled Cross-flow drag (Mq|q|)
ARe     =   6.500000;           % n/a     Fin aspect ratio
dCL     =   1.5*4.130000;       % n/a     Scaled Coef. of Lift Slope
CDc     =   0.030000;           % n/a     Crossflow Coef. of Drag


mass        =   147.8671;         % kg Flooded Vehicle total mass
movableMass =   26;               % kg Battary movable mass
dropWtMass  =   1.0;              % kg Mass of the drop weight #1, kg
dropWtX     =  -0.1330;           % m  X location of the drop weight #1, m

% Account for movable mass shift (x center of gravity)
Xmass = (movableMass.*xstruct.mass_p + dropWtMass*dropWtX)./mass;

% Get control data
ui = controls;

%--------------------------------------------------------------------------
% Correct for offsets in data:
ele_offset =  -1*pi/180;            % found in ele_offsetLRAUV_SIM.m
rud_offset =   0*pi/180;              % found in ele_offsetLRAUV_SIM.m

ui(:,1) = ui(:,1) + ele_offset; % zeros(size(ui(:,1)));
ui(:,2) = ui(:,2) + rud_offset;

% Correct for hysteresis and backlash offsets
%{
for k = 1:length(ui)
    
    if ui(k,1)<0
        ui(k,1) = ui(k,1) -0.65*pi/180;
    elseif ui(k,1)>4*pi/180
        ui(k,1) = ui(k,1) + 1*pi/180;
    end
    
end
%}

%--------------------------------------------------------------------------
% Unpack state vector
x = zeros(1,12);
%
for c=[1:6,9,10:12];
    x(c) = xstruct.(names{c})(startPoint);
end; clear c
%}

%--------------------------------------------------------------------------
% Run
waitbar(0.01,h,'Running Vehicle Simulation...');
for i = startPoint:startPoint+n_steps
    
    % Account for movable mass shift
    xg = Xmass(i) ;
    %{
    % model broken mass shifter: free play in mass position - moves
                                 forwared (back) when pitched down (up).
    if x(11)<-2*pi/180
        xg = Xmass(i) + movableMass*0.000005 ;
    elseif x(11)>2*pi/180
        xg = Xmass(i) - movableMass*0.000005 ;
    else
        xg = Xmass(i) ;
    end
    %}
    
    % Set some vars constant
    x(1) = xstruct.u(i);
    %     x(2) = xstruct.v(i);
    %     x(3) = xstruct.w(i);
    x(4)  = xstruct.p(i);
    %     x(6)  = xstruct.r(i);
    x(10) = xstruct.phi(i);
    %     x(12) = xstruct.psi(i);
    
    ui_in = ui(i,:); % | ui_in = [0, ui(i,2:4)] use to kill elevator input
    
    % Compute fin forces and moments:
    %{
    [ F1, F2, F3, F4, M1, M2, M3, M4 ] = robsFins(ui_in , x );
    fin.X(:,i) = F1(1)+F2(1)+F3(1)+F4(1);
    fin.Y(:,i) = F1(2)+F2(2)+F3(2)+F4(2);
    fin.Z(:,i) = F1(3)+F2(3)+F3(3)+F4(3);
    fin.K(:,i) = M1(1)+M2(1)+M3(1)+M4(1);
    fin.M(:,i) = M1(2)+M2(2)+M3(2)+M4(2);
    fin.N(:,i) = M1(3)+M2(3)+M3(3)+M4(3);
    %}
    % Log step data:
    simlog(i,:) = [x ui_in];
    
    % Calc next step
    [xdot,forces] = lrauv(x,ui_in); % main simulation function
    
    
    % Log outputed forces
    f(:,i) = forces;
    
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
    
    waitbar((find(n==i)/length(n)),h,['Running Vehicle Simulation... ['...
        num2str(100*(find(n==i)/length(n)),2) '%]'] );
    pause(0.01)
end
waitbar(1,h,['Vehicle Simulation Complete [' num2str(100) '%]'] );

lag = -10:10;
for c=1:length(lag)
    lagErr(c,:) = [sum(( xstruct.theta(n+lag(c))-simlog(n,11)' ).^2),...
        lag(c)];
    
end; clear c
[minlag,minlagi] = min(lagErr(:,1));


% PLOT
%--------------------------------------------------------------------------
%
for c=[9,11]
    figure;
    set(gcf,'Units','normalized','Position',[0.05 0.1 0.9 0.8],...  % [left bottom width height]
    'PaperPositionMode','auto');

    s1=subplot(2,1,1);
    p1=plot(xstruct.(names{c})(n+lag(minlagi)),'linewidth',2); hold on;
    p2=plot(simlog(n,c),'linewidth',2);
    title((names{c}),'fontweight','bold','fontsize',22);
    lg=legend('Observed','Modeled','location','nw');
    set(lg,'fontweight','bold','fontsize',16)
    xlabel('Time (sec)','fontweight','bold','fontsize',16)
    set(gca,'xticklabel',(get(gca,'xtick')./(1/time_step)))
    set(gca,'layer','top','FontWeight','Bold','fontSize',14)
    grid on; hold off;
    
    if sum(c == [4:6,10:12])
        set(p1,'ydata',get(p1,'ydata')*180/pi);
        set(p2,'ydata',get(p2,'ydata')*180/pi);
        ylabel('deg','fontweight','bold','fontsize',16)
    elseif c==9
        set(gca,'YDir','reverse');
    end
    
    
    s2=subplot(2,1,2);
    %     bar((simlog(n,11)' - xstruct.theta(n+lag(minlagi)))*180/pi)
    %
    plot(ui(n,1:2)*180/pi,'linewidth',2)
    grid on; hold on;
    plot(zeros(size(ui(n,1))),'k--');
    hold off;
    xlabel('Time (sec)','fontweight','bold','fontsize',16);
    ylabel('deg','fontweight','bold','fontsize',16);
    set(gca,'xticklabel',(get(gca,'xtick')./(1/time_step)))
    lg=legend('Elev ang','Rud ang','location','nw');
    set(lg,'fontweight','bold','fontsize',16)
    set(gca,'layer','top','FontWeight','Bold','fontSize',14)
    %     plot(Xmass(n))
    %}
    
    linkaxes([s1,s2],'x')
end; clear c
%}

% plot forces
%{

forceNames = {'X','Y','Z','K','M','N'};
for c=[3,5]
figure;
subplot(2,1,1)
plot(f(c,n))
title((forceNames{c}),'fontweight','bold','fontsize',16);

subplot(2,1,2)
plot(fin.(forceNames{c})(n))
title(['Fin ' (forceNames{c})],'fontweight','bold','fontsize',16);
end
%}


% Compute error stats: 1) sum(Error)   2) max(Error)
error = [ sum((xstruct.theta(n+lag(minlagi))-simlog(n,11)').^2),...
    max(abs(xstruct.theta(n+lag(minlagi))-simlog(n,11)')*180/pi) ]

% close waitbar
close(h)



% some control plots I was playing around with
%{

figure; 
hist(xstruct.theta(n)-simlog(n,11)',100)

figure;
subplot(2,1,1)
plot(gradient(ui(n,1)*180/pi))
hold on
plot(xstruct.theta(n)*180/pi,'o-');
plot(simlog(n,11)*180/pi,'r.-');
plot(ui(n,1)*180/pi,'g')

subplot(2,1,2)
plot(ui(n,1:2)*180/pi)
xlabel('Time (sec)'); ylabel('deg');
set(gca,'xticklabel',(get(gca,'xtick')./5))
legend('Elev ang','Rud ang');

figure;
plot(controls(:,1))

%}

