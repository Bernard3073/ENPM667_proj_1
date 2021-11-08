% heading angle
phi = 90 * pi / 180;
% wheel radius
r = 0.17;
% distance from the wheel to the center of mass
b = 0.3;
% constant
c = r / 2*b;
% 
d = 0.05;
% look ahead distance
L_a = 0.1;
% wheel angular velocity
theta_dot_r = 0;
theta_dot_l = 0;
% vector of wheel velocities
niu = [theta_dot_r; theta_dot_l];  
% desired velocity
v_d  = [0, 1.414];
mat_A = [-sin(phi) cos(phi) -d 0 0;
            -cos(phi) -sin(phi) -b r 0;
            -cos(phi) -sin(phi) b 0 r;];
mat_S = [c*(b*cos(phi) - d*sin(phi)) c*(b*cos(phi) + d*sin(phi));
        c*(b*sin(phi) + d*cos(phi)) c*(b*sin(phi) - d*cos(phi));
        c -c;
        1 0;
        0 1;];
% sampling time
h = 0.02;
%********straight line path Jacobian ****
A = -1;
B = 1;
C = 0;
J_h_1 = 1/sqrt(A^2 + B^2) * [A, B, B*L_a*cos(phi) - A*L_a*sin(phi), 0, 0];
J_h_2 = [r/2, r/2];
% 2*2 matrix
phi_func = [J_h_1 * mat_S; J_h_2];
phi_func_diff = diff(phi_func);
det_phi_func = r^2*(d + L_a)*(B*cos(phi) - A*sin(phi)) / 2*b*sqrt(A^2 + B^2); 
% error between the desired velocity and the current velocity
v = [v_d(1) - theta_dot_r; v_d(2) - theta_dot_l];
phi_func_inv = inv(phi_func);
u = phi_func_inv * (v - phi_func_diff * niu);
x_dot = [mat_S*niu; 0; 0] + [0 0;0 0;0 0;0 0;0 0;1 0;0 1]*niu;
