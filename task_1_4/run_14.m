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
tsamp=10;           % Sampling time for how often states are stored. (NOT ODE solver time step)
                
p0=zeros(2,1);      % Initial position (NED)
v0=[6.63 0]';       % Initial velocity (body)
psi0=0;             % Inital yaw angle
r0=0;               % Inital yaw rate
c=1;                % Current on (1)/off (0)

nc_max = 85 * 2 * pi / 60; % rad/s
dc_lim = 25 * pi/180; 


% Heading Control Parameters
Kp_psi = 100; 
Kd_r = 1000;  

sim MSFartoystyring % The measurements from the simulink model are automatically written to the workspace.

%pathplotter(p(:,1),p(:,2),psi(:),tsamp,1,tstart,tstop,0,WP);
fig1 = figure(1);
set(fig1, 'Position', [100 300 700 400])
subplot(1,2,1);
plot(t,psi*180/pi,t,psi_d*180/pi,'--',t,psi_e*180/pi,'linewidth',1.5);
xlabel('time');
ylabel('degrees');
legend('\psi','\psi desired','\psi error');
xlim([0,1000]);
grid on
subplot(1,2,2);
plot(t,r*180/pi,t,r_d*180/pi,'--',t,r_e*180/pi,'linewidth',1.5);
xlabel('time');
ylabel('degree/s');
xlim([0,1000]);
legend('r','r desired','r error');
grid on

fig2 = figure(2);
set(fig2, 'Position', [800 300 700 400])
plot(t,dc_in*180/pi,t,dc_lim*ones(1,length(t))*180/pi,t,-dc_lim*ones(1,length(t))*180/pi,'linewidth',1.5);
xlabel('time');
ylabel('rad');
xlim([0,1000]);
legend('\delta_c','upper limit','lower limit');
grid on