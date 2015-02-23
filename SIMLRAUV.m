
function [simlog, f] = SIMLRAUV(time_step, x, xstruct, controls, startPoint,timeEval, Xmass)  
% STATE AND INPUT VECTORS:
% x = [u v w p q r xpos ypos zpos phi theta psi]'
% ui = [ delta_s delta_r Xprop Kprop ]'


% h = waitbar(0,'Initializing LRAUV Vehicle Simulator...');
global xg

n_steps     = fix(timeEval/time_step); 
n = startPoint:startPoint+n_steps;

% Get control data
ui = controls;


% Run
% waitbar(0.01,h,'Runing Vehicle Simulation...');
for i = startPoint:startPoint+n_steps
    
    % Account for movable mass shift
    xg = Xmass(i) ;
    %{
    % simulate broken mass shifter: allow free play of mass (back when 
    %                               pitched up, forward when pitched down)
    if x(11)<-10*pi/180
        xg = Xmass(i) + movableMass*0.000005 ;
    elseif x(11)>10*pi/180
        xg = Xmass(i) - movableMass*0.000005 ;
    else
        xg = Xmass(i) ;
    end
    %}
    
% Overwrite sim outputs with actual historical data
    x(1) = xstruct.u(i);
    % x(2) = xstruct.v(i);
    % x(3) = xstruct.w(i);
    x(4)  = xstruct.p(i);
    % x(6)  = xstruct.r(i);  
    x(10) = xstruct.phi(i);
    % x(12) = xstruct.psi(i);
    
    ui_in = [0, ui(i,2:4)];

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
    [xdot,forces] = lrauv(x,ui_in); % LRAUV % lrauvTest
    
    
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
    
%     waitbar((find(n==i)/length(n)),h,['Runing Vehicle Simulation... ['...
%         num2str(100*(find(n==i)/length(n)),2) '%]'] );
    
end
% waitbar(1,h,['Vehicle Simulation Complete [' num2str(100) '%]'] );

simlog = simlog(n,:);
% close(h)
