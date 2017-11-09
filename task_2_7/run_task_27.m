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

%% Simulation parameters

tstart=0;           % Sim start time
tstop=5000;        % Sim stop time
tsamp=100;           % Sampling time for how often states are stored. (NOT ODE solver time step)

%% Helping Functions
rad2grad = 180/pi;
grad2rad = pi/180;

%% Ship
% Initial Conditions
p0=[1500;500];      % Initial position (NED)
v0=[6.63 0]';       % Initial velocity (body)
psi0=deg2rad(50);             % Inital yaw angle
r0=0;               % Inital yaw rate
c=1;                % Current on (1)/off (0)

nc_max = 85 * 2 * pi / 60; % rad/s
dc_lim = 25 * grad2rad; 

%% Control Parameters
% Velocity Control Parameters
Kp_u = 2.5;
Ki_u = 0.25;
windup_gain = 1;

% Heading Control Parameters
Kp_psi = 10;
Kd_r = 1500;

% Course control Parameters
Kp_chi = 0;
Ki_chi = 0.00;
T_psi = 50;

%% target tracking
U_aMax = 5; %Eirik: Dette tallet tok jeg vilkårlig 
delta_pTilde = 3000; %Eirik: Sier noe om når båten kjører fort. Et lavt tall fører
                    % til at den gasser på mye når interceptor nærmer seg
                    % target. Vilkårlig valgt.
nu_t = [-2.4412;
        -1.7437];
p0_t = [-3500;
        -2500];

sim MSFartoystyring % The measurements from the simulink model are automatically written to the workspace.

%% Plotting

plot_time = 5000;

WP = [0 0;-3500 -2500]';
pathplotter(p(:,1),p(:,2),psi(:),tsamp,1,tstart,tstop,1,WP);

figure(4)
plot(t,rad2deg(chi_d),t,rad2deg(normalize(chi(:))),t,rad2deg(chi_e),'linewidth',1.5);
title('Ship Course')
xlabel('time');
ylabel('degrees');
xlim([0 plot_time]);
legend('\chi_d','\chi','\chi_e');
grid on
%{
figure(5);
plot(t,v(:,1),t,u_d,'--',t,u_e,'linewidth',1.5);
xlabel('time');
ylabel('m/s');
xlim([0,plot_time]);
legend('u','u desired','u error');
grid on

figure(6);
plot(t,r*rad2grad,t,rad2deg(mod(psi_d,2*pi)),t,rad2deg(mod(psi,2*pi)),t,rad2deg(psi_e),'linewidth',1.5);
xlabel('time');
ylabel('degree, degree/s');
xlim([0,plot_time]);
legend('r','\psi_d','\psi','\psi_e');
grid on

figure(7);
plot(t,nc_out,t,nc_max*ones(1,length(t)),t,-nc_max*ones(1,length(t)),'linewidth',1.5);
xlabel('time');
ylabel('rad/s');
xlim([0,plot_time]);
legend('n_c','upper limit','lower limit');
grid on
%}