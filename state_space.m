function [q_dot, pos_x_hat, pos_y,  vel_x_hat, vel_y] = state_space(t, q, pos_x_hat, pos_y,  vel_x_hat, vel_y)
x_c = q(1);
y_c = q(2);
% heading angle
phi = q(3);
% wheel angular velocity
theta_dot_r = q(6);
theta_dot_l = q(7);

% wheel radius
r = 0.17;
% distance from the wheel to the center of mass
b = 0.30;
% constant
c = r / 2*b;
% Distance between Po and Pc
d = 0.05;
% look ahead distance
L_a = 0.1;

% vector of wheel velocities
niu = [theta_dot_r; theta_dot_l];  
% desired velocity
v_d  = [0, 1.414];
% Mass of chassis
m_c = 60;
% Mass of each wheel + rotor of motor 
m_w = 1;
% Total mass of the robot
m_t = m_c + 2*m_w;
% Inertia of chassis
I_c = 15.625;
% Inertia of Wheels
I_w = 0.005;
% Inertia about a defined axis in the plane of wheel
I_m = 0.0025;
% I in matrix
I = ((m_c * d^2) + 2*m_w*(b^2 + d^2) + I_c + 2*I_m);

% error between the desired velocity and the current velocity
% v = [v_d(1) - theta_dot_r; v_d(2) - theta_dot_l];

% ________Path related variables______%
% Radius of circular path
R = 7.50;


% Steady State Kalman Gains
K_k = [0.0683; 0.0965; 0.0965];


mat_A = [-sin(phi) cos(phi) -d 0 0;
            -cos(phi) -sin(phi) -b r 0;
            -cos(phi) -sin(phi) b 0 r];
mat_S = [c*(b*cos(phi) - d*sin(phi)) c*(b*cos(phi) + d*sin(phi));
        c*(b*sin(phi) + d*cos(phi)) c*(b*sin(phi) - d*cos(phi));
        c -c;
        1 0;
        0 1];
% sampling time
h = 0.02;
%********straight line path Jacobian ****
A = -1;
B = 1;
C = 0;

J_h_1 = (1/sqrt(A^2 + B^2)).*[A, B, (B*L_a*cos(phi) - A*L_a*sin(phi)), 0, 0];
J_h_2 = [r/2, r/2];
% 2*2 matrix
phi_func = [J_h_1 * mat_S; J_h_2];
phi_func_diff = diff(phi_func);
det_phi_func =((r.^2)*(d + L_a)*(B*cos(phi) - A*sin(phi)))/(2*b*sqrt(A^2 + B^2));
det_straightpath = det(phi_func);

phi_func_inv = inv(phi_func);


% the coordinates of the point P_l 
x_l = x_c + L_a * cos(phi);
y_l = y_c + L_a * sin(phi);


% discretization of position model
% matrix for the discrete time state space
A_pos =[0, 1; 0, 0];
B_pos = [0;1];
C_pos = [1, 0];
sys_pos = ss(A_pos, B_pos, C_pos, 0);
sys_pos_discrete= c2d(sys_pos, h);
phi_r_pos = sys_pos_discrete.A;
gamma_r_pos = sys_pos_discrete.B;
C_r_pos = sys_pos_discrete.C;
% discretization of velocity model
A_vel = [1 h; 0 0];
B_vel = [0; 1];
C_vel = [1 0];
sys_vel = ss(A_vel, B_vel, C_vel, 0);
sys_vel_discrete = c2d(sys_vel, h);
phi_r_vel = sys_vel_discrete.A;
gamma_r_vel = sys_vel_discrete.B;
C_r_vel = sys_vel_discrete;

% C_a = [C_r, 0]
% y_k_hat = C_a * [phi_func_r - gamma_r*L, 0; 0, 1]*
Q_x_r = 10^(-5)*eye(4)
Q_p_k = 10^(-4)
R_k = 10^(-3)
%P_0 = [Q_x_r  0; ]
% Kalman gain for position model
K_k_pos = [0.1407, 0.3018, 0.2931, 0.2931, 0.2931]'
% Kalman gain for velocity model
K_k_vel = [0.0683 0.0965 0.0965]'
phi_1 = [1, h; 0, 1]
gamma_1 = [h^2/2; h]
% pole of the position model
P = pole(sys_pos)
L = acker(phi_r_pos, gamma_r_pos, P) 

% position update
% v_1_update =  K_k*q_0


pos_AOB_A = zeros(5, 5);
pos_AOB_A(1:2, 1:2) = phi_r_pos - gamma_r_pos*L;
pos_AOB_A(5, 5) = 1;
vel_AOB_A = zeros(3, 3);
vel_AOB_A(1:2, 1:2) = phi_r_pos - gamma_r_pos*L;
vel_AOB_A(3, 3) = 1;
% AOB_B = [gamma_r; 0]
C_a_pos = [1 0 0 0 0];
C_a_vel = [1 0 0];

pos_x_hat(:, t+1) = pos_AOB_A*pos_x_hat(:, t) + K_k_pos*(pos_y(t) - C_a_pos*pos_x_hat(:, t));
vel_x_hat(:, t+1) = vel_AOB_A*vel_x_hat(:, t) + K_k_vel*(vel_y(t) - C_a_vel*vel_x_hat(:, t));
v = [pos_x_hat(1, t+1); vel_x_hat(1, t+1);];
h1= (-x_l + y_l)/2;
h2 = r*(v(1) + v(2))/2;
y_k = [h1, h2];
u = phi_func_inv * (v - phi_func_diff * niu);
x_dot = [mat_S*niu; 0; 0] + [0 0;0 0;0 0;0 0;0 0;1 0;0 1]*u;
q_dot = x_dot';