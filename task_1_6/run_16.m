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
v0=[0.001 0]';       % Initial velocity (body)
psi0=0;             % Inital yaw angle
r0=0;               % Inital yaw rate
c=0;                % Current on (1)/off (0)

% Velocity Dynamics 
num = [0 0.00142452745428828 8.50704549113513e-06]; %from system identification toolbox
den = [1 0.00589033596766649 8.66512373642295e-06]; %from system identification toolbox
vd = tf(num,den);

nc_max = 85 * 2 * pi / 60; % rad/s
dc_lim = 25 * pi/180; 
nc_d = 7.3;


% Heading Control Parameters
Kp_psi = 100; 
Kd_r = 1000;  


sim MSFartoystyring % The measurements from the simulink model are automatically written to the workspace.

%pathplotter(p(:,1),p(:,2),psi(:),tsamp,1,tstart,tstop,0,WP);
fig1 = figure(1);
set(fig1, 'Position', [100 300 700 400])
plot(t,v(:,1),'linewidth',1.5);
xlabel('time');
ylabel('m/s');
xlim([0,1000]);
legend('velocity u');
grid on

fig2 = figure(2);
set(fig2,'Position', [800 300 700 400])
plot(t,nc_in,t,nc_max*ones(1,length(t)),'linewidth',1.5);
xlabel('time');
ylabel('rad/s');
xlim([0,1000]);
legend('n_c','upper limit','lower limit');
grid on