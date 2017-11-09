function course = desiredCourse(u)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

p = u(1:2);
lookahead_distance = u(3);

persistent k;
persistent wp;
<<<<<<< HEAD:task_2_7/desiredCourse.m
k = 1;


=======

if isempty(wp)
    wp = load('WP.mat');
    k = 1;
end
>>>>>>> 4abd72938ff1ed5c8d61a2dcf99f7683b06b2681:task_2_4/desiredCourse.m

if sqrt((p(1)-wp.WP(1,k+1))^2+(p(2)-wp.WP(2,k+1))^2) < lookahead_distance && k < 5
    k = k+1;
end    

% Angle of the path
x_0 = wp.WP(1,k);
y_0 = wp.WP(2,k);
x_1 = wp.WP(1,k+1);
y_1 = wp.WP(2,k+1);
a_k = atan2(y_1 - y_0, x_1 - x_0);

Rot = [cos(a_k) -sin(a_k);
       sin(a_k)  cos(a_k)];

% Error in boat position
epsilon = Rot' * (p - wp.WP(:,k));
s_error = epsilon(1);
e_error = epsilon(2);

% Calculation of desired course
course = a_k + atan2(-e_error,lookahead_distance);
%course = normalize(course);
end

