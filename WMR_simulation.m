format default;
clear all
% initial state
q_0 = [0, 0, pi/2, 0, 0, 0, 0];

t_0 = 0;
t = 20;
timestep = 0.5;
pos_x_hat = zeros(5, t/timestep+1);
vel_x_hat = zeros(3, t/timestep+1);
pos_x = zeros(5, t/timestep+1);
vel_x = zeros(3, t/timestep+1);
pos_y = zeros(5, t/timestep+1);
vel_y = zeros(3, t/timestep+1);

% pos_x_hat = zeros(2, t/timestep+1);
% vel_x_hat = zeros(2, t/timestep+1);
% pos_x = zeros(2, t/timestep+1);
% vel_x = zeros(2, t/timestep+1);
% pos_y = zeros(2, t/timestep+1);
% vel_y = zeros(2, t/timestep+1);
i = 1;
x_c_list = [];
y_c_list = [];
while t_0 < t
    [q_dot, pos_x_hat, pos_x, pos_y,  vel_x_hat, vel_x, vel_y] = state_space(i, q_0, pos_x_hat, pos_x, pos_y,  vel_x_hat, vel_x, vel_y);
    q_0 = q_0 + q_dot*t_0;
    i = i + 1;
    t_0 = t_0 + timestep;
    x_c_list = [x_c_list q_0(1)];
    y_c_list = [y_c_list q_0(2)];
end
plot(x_c_list, y_c_list);

