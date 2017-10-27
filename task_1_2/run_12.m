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
clear all
close all

%%
tstart=0;           % Sim start time
tsamp=10;           % Sampling time for how often states are stored. (NOT ODE solver time step)
                
p0=zeros(2,1);      % Initial position (NED)
v0=[6.63 0]';       % Initial velocity (body)
psi0=0;             % Inital yaw angle
r0=0;               % Inital yaw rate
c=0;                % Current on (1)/off (0)

K = zeros(3,1);
T = zeros(3,1);
ratio = zeros(10,1);
timesteps = [1000,10000,100000]

i = 1;

for n = timesteps
    tstop=n;        % Sim stop time
    dc = 1*pi/180;
    
sim MSFartoystyring % The measurements from the simulink model are automatically written to the workspace.

%pathplotter(p(:,1),p(:,2),psi(:),tsamp,1,tstart,tstop,0,WP);

% r_deriv = r(2)/t(2);
% tNy = 0:1:200;
% T_vector = tNy*r_deriv;
% steadyState = r(length(r));
% 
% Tc(delta_c) = steadyState/r_deriv;
% K(delta_c) = steadyState * Tc(delta_c) / dc;
% ratio(delta_c) = K(delta_c)/Tc(delta_c);

% Derive K and T for dc:
tdata = t;
rdata = r*180/pi;

% x(1) = 1/T; x(2) = K;
x0 = [0.01, 0.1]';
F = inline('exp(-tdata*x(1))*0 + x(2)*(1-exp(-tdata*x(1)))*5','x','tdata');
x = lsqcurvefit(F,x0,tdata,rdata);

T(i) = 1/x(1);
K(i) = x(2);
figure(i)
plot(t,r,t,psi);
legend('r','psi','dc');

% figure(delta_c)
% plot(t,r);
% hold on
% plot(tNy,T_vector);
% hold on
% plot(ones(1,200)*steadyState);
% legend('yaw rate','timeconstantvector','steadystate');
 i = i + 1;
end
disp('25 grader');
disp('K:');
disp(K);
disp('T:')
disp(T);