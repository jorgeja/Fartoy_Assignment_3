function course = desiredCourse(u)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

p = u(1:2);
lad = u(3);

persistent k;
persistent wp;
persistent first;
k = 1;

if first == 0
    wp = load('WP.mat');
    first = 1;
end

% Angle of the path
x_0 = wp(1,k);
y_0 = wp(2,k);
x_1 = wp(1,k+1);
y_1 = wp(2,k+1);
a_k = atan2(y_1 - y_0, x_1 - x_0);

Rot = [cos(a_k) -sin(a_k);
       sin(a_k)  cos(a_k)];

% Error in boat position
epsilon = Rot' * (p - wp(:,k));
s_error = epsilon(1);
e_error = epsilon(2);

% Calculation of desired course
course = a_k + atan(-e_error/lad);
end

