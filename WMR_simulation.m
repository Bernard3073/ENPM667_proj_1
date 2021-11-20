format default;
% heading angle
phi = 90 * pi / 180;
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
% wheel angular velocity
theta_dot_r = 0;
theta_dot_l = 0;
% vector of wheel velocities
niu = [theta_dot_r; theta_dot_l];  
% Initial velocity
v_0  = [0, 0];
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
x_c = 0;
y_c = 0;
J_h_1 = (1/sqrt(A^2 + B^2)).*[A, B, (B*L_a*cos(phi) - A*L_a*sin(phi)), 0, 0];
J_h_2 = [r/2, r/2];
% 2*2 matrix
phi_func = [J_h_1 * mat_S; J_h_2];
phi_func_diff = diff(phi_func);
det_phi_func =((r.^2)*(d + L_a)*(B*cos(phi) - A*sin(phi)))/(2*b*sqrt(A^2 + B^2));
det_straightpath = det(phi_func);
% error between the desired velocity and the current velocity
v = [v_d(1) - theta_dot_r; v_d(2) - theta_dot_l];
phi_func_inv = inv(phi_func);
u = phi_func_inv * (v - phi_func_diff * niu);
x_dot = [mat_S*niu; 0; 0] + [0 0;0 0;0 0;0 0;0 0;1 0;0 1]*niu;

% the coordinates of the point P_l 
x_l = x_c + L_a * cos(phi);
y_l = y_c + L_a * sin(phi);
h1= (-x_l + y_l)/2;
h2 = r*(v(1) + v(2))/2;

% position model
y_k = [h1, h2];
C_r = [1,  0, 0, 0];
C_a = [C_r, 0];
y_k_hat = C_a * 
Q_x_r = 10^(-5)*eye(4);
Q_p_k = 10^(-4);
R_k = 10^(-3);
P_0 = [Q_x_r  0; ];
% Kalman gain
K_k = [0.1407, 0.3018, 0.2931, 0.2931, 0.2931].T;
phi_1 = [1, h; 0, 1];
gamma_1 = [h^2/2; h];
% plot a straight line
    b0 = 3;
    b1 = 4;
    f = @(x) b0+b1*x;
    ezplot( f, 0, 5 );
