%% Information 
% This file is only an example of how you can start the simulation. The
% sampling time decides how often you store states. The execution  time
% will increase if you reduce the sampling time.

% Please note that the file "pathplotter.m" (only used in the second part
% of the assignment) shows the ship path during the path following and
% target tracking part of the assignment. It can be clever to adjust the sampling
% time when you use that file because it draws a sketch of the ship in the
% North-East plane at each time instant. Having a small sampling time will
% lead to a lot over overlap in the ship drawing. 

% You should base all of your simulink models on the MSFartoystyring model
% and extend that as you solve the assignment. For your own sake, it is
% wise to create a new model and run file for each task. That is
% especially important in the problems you need to hand in since the files
% you deliver only should create the desired result in that task.

% The msfartoystyring.m file includes the ship model. You are not allowed
% to change anything within that file. You need to include that file in
% every folder where you have a simulink model based on
% "MSFartoystyring.slx". 

% WP.mat is a set of six waypoints that you need to use in the second part of
% the assignment. The north position is given in the first row and the east
% position in the second row. 

close all;
clear all;

%%
tstart=0;           % Sim start time
tstop=10000;        % Sim stop time
tsamp=100;           % Sampling time for how often states are stored. (NOT ODE solver time step)

% Helping Functions
rad2grad = 180/pi;
grad2rad = pi/180;

% Initial Conditions
p0=[1500;500];      % Initial position (NED)
v0=[6.63 0]';       % Initial velocity (body)
psi0=deg2rad(50);             % Inital yaw angle
r0=0;               % Inital yaw rate
c=1;                % Current on (1)/off (0)

% Velocity Dynamics 
num = [0 0.00142452745428828 8.50704549113513e-06]; %from system identification toolbox
den = [1 0.00589033596766649 8.66512373642295e-06]; %from system identification toolbox
vd = tf(num,den);

nc_max = 85 * 2 * pi / 60; % rad/s
dc_lim = 25 * grad2rad; 

% Velocity Control Parameters
Kp_u = 300;
Ki_u = 0.025;

% Heading Control Parameters
Kp_psi = 1;
Ki_psi = 0;         %1/1000*Kp_psi;
Kd_r = 100;  

% Load Waypoints
%load('WP.mat');
look_ahead_distance = 100; % Eirik: Ifølge example 10.1 er denne vanligvis
                           % mellom 1,5 og 2,5 ganger lengden på skipet (som er 304,8m).
                           % Burde endre verdien til 500 eller noe. Kan
                           % forkortes Lpp.

% 2_7 script
U_aMax = 20; %Eirik: Dette tallet tok jeg vilkårlig 
delta_pTilde = 300; %Eirik: Sier noe om når båten kjører fort. Et lavt tall fører
                    % til at den gasser på mye når interceptor nærmer seg
                    % target. Vilkårlig valgt.
nu_t = [-2.4412;
        -1.7437];
p0_t = [-3500;
        -2500];

sim MSFartoystyring % The measurements from the simulink model are automatically written to the workspace.

%Plotting

plot_time = 2000;

WP = [0 0;-3500 -2500]';
fig5 = figure(5);
set(fig5, 'Position', [100 50 700 400])
pathplotter(p(:,1),p(:,2),psi(:),tsamp,1,tstart,tstop,1,WP);

fig3 = figure(3);
set(fig3, 'Position', [100 50 700 400])
plot(t,v(:,1),t,u_d,'--',t,u_e,'linewidth',1.5);
xlabel('time');
ylabel('m/s');
xlim([0,plot_time]);
legend('u','u desired','u error');
grid on

fig4 = figure(4);
set(fig4, 'Position', [100 400 700 400])
plot(t,r*rad2grad,t,psi*rad2grad,'linewidth',1.5);
xlabel('time');
ylabel('degree, degree/s');
xlim([0,plot_time]);
legend('r','\psi');
grid on

fig5 = figure(5);
set(fig5,'Position', [800 400 700 400])
plot(t,nc_out,t,nc_max*ones(1,length(t)),t,-nc_max*ones(1,length(t)),'linewidth',1.5);
xlabel('time');
ylabel('rad/s');
xlim([0,plot_time]);
legend('n_c','upper limit','lower limit');
grid on