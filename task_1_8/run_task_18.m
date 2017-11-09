close all;
clear all;

%%
tstart=0;           % Sim start time
tstop=10000;        % Sim stop time
tsamp=10;           % Sampling time for how often states are stored. (NOT ODE solver time step)

% Helping Functions
rad2grad = 180/pi;
grad2rad = pi/180;

% Initial Conditions
p0=zeros(2,1);      % Initial position (NED)
v0=[0.01 0]';       % Initial velocity (body)
psi0=0;             % Inital yaw angle
r0=0;               % Inital yaw rate
c=1;                % Current on (1)/off (0)

nc_max = 85 * 2 * pi / 60; % rad/s
dc_lim = 25 * grad2rad; 
u_d = 7;

% Velocity Control Parameters
Kp_u = 100;
Ki_u = 1;

% Heading Control Parameters
Kp_psi = 10;
Ki_psi = 0;         %1/1000*Kp_psi;
Kd_r = 1000;  


sim MSFartoystyring % The measurements from the simulink model are automatically written to the workspace.

%pathplotter(p(:,1),p(:,2),psi(:),tsamp,1,tstart,tstop,0,WP);

%Plotting

plot_time = 2000;

fig1 = figure(1);
set(fig1, 'Position', [100 50 700 400])
plot(t,v(:,1),t,u_desired,'--',t,u_e,'linewidth',1.5);
xlabel('time');
ylabel('m/s');
xlim([0,plot_time]);
legend('u','u desired','u error');
grid on

fig2 = figure(2);
set(fig2, 'Position', [100 400 700 400])
plot(t,r*rad2grad,t,psi*rad2grad,'linewidth',1.5);
xlabel('time');
ylabel('degree, degree/s');
xlim([0,plot_time]);
legend('r','\psi');
grid on

fig3 = figure(3);
set(fig3,'Position', [800 400 700 400])
plot(t,nc_out,t,nc_max*ones(1,length(t)),t,-nc_max*ones(1,length(t)),'linewidth',1.5);
xlabel('time');
ylabel('rad/s');
xlim([0,plot_time]);
legend('n_c','upper limit','lower limit');
grid on