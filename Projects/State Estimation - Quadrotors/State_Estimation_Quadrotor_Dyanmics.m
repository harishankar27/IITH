%% Quadrotor Nonlinear Simulation and State Estimation
%
% This script simulates the full nonlinear dynamics of a quadrotor UAV
% using specified physical parameters and control inputs. It numerically
% integrates the equations of motion over time with `ode45` to generate the
% vehicle's trajectory (positions, velocities, Euler angles, and angular
% rates).

%% System and Inputs
system_params = struct();

system_params.b_torque = 2.92 * 1e-9;
system_params.L = 0.222;
system_params.k_motor = 1.49 * 1e-7;
system_params.m = 1.023;  
system_params.kd = 0.1;


system_params.g = 9.81;   
system_params.Ixx = 9.5 * 1e-3;
system_params.Iyy = 9.5 * 1e-3;
system_params.Izz = 18.6 * 1e-3;
system_params.I = diag([system_params.Ixx, system_params.Iyy, system_params.Izz]);

t_start = 0;
t_stop = 12;
t_step = 0.005;
t_span = t_start:t_step:t_stop;
t1 = t_start:t_step:t_stop/8 ;

omega = 2*pi*1;
operating_point_input = ((system_params.m*system_params.g)/(4*system_params.k_motor));
input1 = zeros(1,length(t_span));
input3 = zeros(1,length(t_span));

part21 = 0.0025*operating_point_input(1)*(sin(omega*t1).^2).*exp(-10*t1);
part22 = -0.0025*operating_point_input(1)*(sin(omega*t1).^2).*exp(-10*t1);
part23 = -0.0025*operating_point_input(1)*(sin(omega*t1).^2).*exp(-10*t1);
part24 =  0.0025*operating_point_input(1)*(sin(omega*t1).^2).*exp(-10*t1);


part41 = -0.0025*operating_point_input(1)*(sin(omega*t1).^2).*exp(-10*t1);
part42 = 0.0025*operating_point_input(1)*(sin(omega*t1).^2).*exp(-10*t1);
part43 = 0.0025*operating_point_input(1)*(sin(omega*t1).^2).*exp(-10*t1);
part44 =  -0.0025*operating_point_input(1)*(sin(omega*t1).^2).*exp(-10*t1);

input2 = horzcat( part21, part22(1:end-1),part23(1:end-1),part24(1:end-1),-part21(1:end-1),-part22(1:end-1),-part23(1:end-1),-part24(1:end-1) );
input4 = horzcat(part41, part42(1:end-1),part43(1:end-1),part44(1:end-1),-part41(1:end-1),-part42(1:end-1),-part43(1:end-1),-part44(1:end-1));

delta_state = zeros(12,1);
operating_point_states = zeros(12,1);
operating_point_states(1) = 1;
operating_point_states(2) = 1;
operating_point_states(3) = 1;

x0 = operating_point_states + delta_state;
delta_input = transpose( vertcat( input1,input2,input3,input4 ) );
u_input  =  operating_point_input +  delta_input;


state_noise = 0.1*diag([5e-4 5e-4 5e-4 5e-4 5e-4 5e-4 5e-4 5e-4 5e-4 5e-4 5e-4 5e-4]);       %1e-5* eye(12);

BdQBdT = state_noise.^2;


figure;
subplot(2, 2, 1);
plot(t_span, u_input(:,1), 'r', 'LineWidth', 1.5); 
title('Motor 1  ');
xlabel('Time (s)');
ylabel('$\omega^2$', 'Interpreter', 'latex');
grid on;


subplot(2, 2, 2);
plot(t_span, u_input(:,2), 'r', 'LineWidth', 1.5); 
title('Motor 2  ');
xlabel('Time (s)');
ylabel('$\omega^2$', 'Interpreter', 'latex');
grid on;

subplot(2, 2, 3);
plot(t_span,u_input(:,3), 'r', 'LineWidth', 1.5); 
title('Motor 3 ');
xlabel('Time (s)');
ylabel('$\omega^2$', 'Interpreter', 'latex');
grid on;


subplot(2, 2, 4);
plot(t_span, u_input(:,4), 'r', 'LineWidth', 1.5); 
title('Motor 4');
xlabel('Time (s)');
ylabel('$\omega^2$', 'Interpreter', 'latex');
grid on;


sgtitle('Quadrotor: Input to the system');
[lc,ld,state_noise_covar] = linearize_quadrotor(system_params,BdQBdT,t_step);

%% Model
[t, x_nc] =  ode45(@(t, x) quadrotor_nonlinear(t, x, u_input,t_span,BdQBdT,system_params), t_span, x0);

position = x_nc(:, 1:3); 
velocity = x_nc(:, 4:6); 
orientation = x_nc(:, 7:9);  
angular_velocity = x_nc(:, 10:12);  

figure;
subplot(2, 2, 1);
plot(t, position(:, 1), 'r', 'LineWidth', 1.5); 
hold on;
plot(t, position(:, 2), 'g', 'LineWidth', 1.5);
plot(t, position(:, 3), 'b', 'LineWidth', 1.5);
title('Position  ');
xlabel('Time (s)');
ylabel('Position (m)');
legend('x', 'y', 'z');

subplot(2, 2, 2);
plot(t, velocity(:, 1), 'r', 'LineWidth', 1.5); 
hold on;
plot(t, velocity(:, 2), 'g', 'LineWidth', 1.5);
plot(t, velocity(:, 3), 'b', 'LineWidth', 1.5);
title('Velocity  ');
xlabel('Time (s)');
ylabel('Velocity (m/s)');
legend('vx', 'vy', 'vz');

subplot(2, 2, 3);
plot(t, orientation(:, 1), 'r', 'LineWidth', 1.5); 
hold on;
plot(t, orientation(:, 2), 'g', 'LineWidth', 1.5);
plot(t, orientation(:, 3), 'b', 'LineWidth', 1.5);
title('Orientation (Roll, Pitch, Yaw)  ');
xlabel('Time (s)');
ylabel('Angle (rad)');
legend('Roll', 'Pitch', 'Yaw');

subplot(2, 2, 4);
plot(t, angular_velocity(:, 1), 'r', 'LineWidth', 1.5); 
hold on;
plot(t, angular_velocity(:, 2), 'g', 'LineWidth', 1.5);
plot(t, angular_velocity(:, 3), 'b', 'LineWidth', 1.5);
title('Angular Velocity (p, q, r)  ');
xlabel('Time (s)');
ylabel('Angular Velocity (rad/s)');
legend('p', 'q', 'r');

sgtitle('Quadrotor Dynamics: Non linear Continous Model');

x_nd = zeros(length(t_span),12);
x_nd(1,:) =  operating_point_states+ delta_state;

for i=2:length(x_nd)

    input = transpose(u_input(i-1,:));
    noise = transpose(mvnrnd(zeros(12,1), state_noise_covar));
    
    x_nd(i,:) =transpose(quadrotor_update(t_step,transpose(x_nd(i-1,:)),input,system_params) +  noise);
end


position = x_nd(:, 1:3);  
velocity = x_nd(:, 4:6);  
orientation = x_nd(:, 7:9);  
angular_velocity = x_nd(:, 10:12);  


figure;
subplot(2, 2, 1);
plot(t, position(:, 1), 'r', 'LineWidth', 1.5); 
hold on;
plot(t, position(:, 2), 'g', 'LineWidth', 1.5);
plot(t, position(:, 3), 'b', 'LineWidth', 1.5);
title('Position  ');
xlabel('Time (s)');
ylabel('Position (m)');
legend('x', 'y', 'z');


subplot(2, 2, 2);
plot(t, velocity(:, 1), 'r', 'LineWidth', 1.5); 
hold on;
plot(t, velocity(:, 2), 'g', 'LineWidth', 1.5);
plot(t, velocity(:, 3), 'b', 'LineWidth', 1.5);
title('Velocity  ');
xlabel('Time (s)');
ylabel('Velocity (m/s)');
legend('vx', 'vy', 'vz');

subplot(2, 2, 3);
plot(t, orientation(:, 1), 'r', 'LineWidth', 1.5); 
hold on;
plot(t, orientation(:, 2), 'g', 'LineWidth', 1.5);
plot(t, orientation(:, 3), 'b', 'LineWidth', 1.5);
title('Orientation (Roll, Pitch, Yaw)  ');
xlabel('Time (s)');
ylabel('Angle (rad)');
legend('Roll', 'Pitch', 'Yaw');

subplot(2, 2, 4);
plot(t, angular_velocity(:, 1), 'r', 'LineWidth', 1.5); 
hold on;
plot(t, angular_velocity(:, 2), 'g', 'LineWidth', 1.5);
plot(t, angular_velocity(:, 3), 'b', 'LineWidth', 1.5);
title('Angular Velocity (p, q, r)  ');
xlabel('Time (s)');
ylabel('Angular Velocity (rad/s)');
legend('p', 'q', 'r');

sgtitle('Quadrotor Dynamics: Non linear Discrete Model');

x_ld = zeros(length(t_span),12);
x_ld(1,:) = delta_state;
for i=2:length(x_ld)

    input = transpose(delta_input(i-1,:));

    in = ld.B*input;

    noise  = transpose(mvnrnd(zeros(12,1), state_noise_covar));

    x_ld(i,:) = transpose(  ld.A*transpose(x_ld(i-1,:)) + in + noise);
end


position = 1 + x_ld(:, 1:3);  
velocity = x_ld(:, 4:6);  
orientation = x_ld(:, 7:9);  
angular_velocity = x_ld(:, 10:12);  

figure;
subplot(2, 2, 1);
plot(t, position(:, 1), 'r', 'LineWidth', 1.5); 
hold on;
plot(t, position(:, 2), 'g', 'LineWidth', 1.5);
plot(t, position(:, 3), 'b', 'LineWidth', 1.5);
title('Position  ');
xlabel('Time (s)');
ylabel('Position (m)');
legend('x', 'y', 'z');


subplot(2, 2, 2);
plot(t, velocity(:, 1), 'r', 'LineWidth', 1.5); 
hold on;
plot(t, velocity(:, 2), 'g', 'LineWidth', 1.5);
plot(t, velocity(:, 3), 'b', 'LineWidth', 1.5);
title('Velocity  ');
xlabel('Time (s)');
ylabel('Velocity (m/s)');
legend('vx', 'vy', 'vz');


subplot(2, 2, 3);
plot(t, orientation(:, 1), 'r', 'LineWidth', 1.5); 
hold on;
plot(t, orientation(:, 2), 'g', 'LineWidth', 1.5);
plot(t, orientation(:, 3), 'b', 'LineWidth', 1.5);
title('Orientation (Roll, Pitch, Yaw)  ');
xlabel('Time (s)');
ylabel('Angle (rad)');
legend('Roll', 'Pitch', 'Yaw');


subplot(2, 2, 4);
plot(t, angular_velocity(:, 1), 'r', 'LineWidth', 1.5); 
hold on;
plot(t, angular_velocity(:, 2), 'g', 'LineWidth', 1.5);
plot(t, angular_velocity(:, 3), 'b', 'LineWidth', 1.5);
title('Angular Velocity (p, q, r)  ');
xlabel('Time (s)');
ylabel('Angular Velocity (rad/s)');
legend('p', 'q', 'r');


sgtitle('Quadrotor Dynamics: Linear Discrete Model');

C1 = zeros(6,12);
C1(1,1) = 1;
C1(2,2) = 1;
C1(3,3) = 1;

C1(4,7) = 1;
C1(5,8) = 1;
C1(6,9) = 1;

y_nc = zeros(length(t_span),6);
y_nd = zeros(length(t_span),6);
y_ld = zeros(length(t_span),6);

sensor_noise = 1e-4*eye(6);

R = sensor_noise.^2   ;

for i=1:length(t_span)

    noise  = transpose(mvnrnd(zeros(6,1), R));

    y_nc(i,:) = transpose(C1*transpose(x_nc(i,:)) + noise);

    noise = transpose(mvnrnd(zeros(6,1), R));

    y_nd(i,:) = transpose(C1*transpose(x_nd(i,:)) + noise);

    noise = transpose(mvnrnd(zeros(6,1), R));

    y_ld(i,:) = transpose(C1*transpose(x_ld(i,:)) + noise);

end

%% Kalman Filter
x_initial_estimate = zeros(12,1);
P_initial = 10*eye(12);
end_sim = length(t_span);

x_kalman_estimate = zeros(end_sim,12);
tic;
x_kalman_estimate(1, :) = x_initial_estimate;
P = P_initial;
for index = 1:end_sim-1
    u =  transpose(delta_input(i,:));
    [x_est, P_update] = Kalman(transpose(y_ld(index + 1, :)), R, ld.A, ld.B, ld.C, P, transpose(x_kalman_estimate(index, :)), u, state_noise_covar);
    x_kalman_estimate(index + 1, :) = x_est;
    P = P_update;
end
elapsed_time = toc;
% Display result (optional)
fprintf('Time Required for Kalman Filter: %.4f seconds\n', elapsed_time)
position = x_kalman_estimate(:, 1:3);  
velocity = x_kalman_estimate(:, 4:6);  
orientation = x_kalman_estimate(:, 7:9);  
angular_velocity = x_kalman_estimate(:, 10:12); 

position_actual = x_ld(:, 1:3);  
velocity_actual = x_ld(:, 4:6);  
orientation_actual = x_ld(:, 7:9);  
angular_velocity_actual = x_ld(:, 10:12); 
figure;


subplot(2, 2, 1);
plot(t, position(:, 1), 'r', 'LineWidth', 1.5); hold on;
plot(t, position(:, 2), 'g', 'LineWidth', 1.5);
plot(t, position(:, 3), 'b', 'LineWidth', 1.5);
title('Position (x, y, z)');
xlabel('Time (s)');
ylabel('Position (m)');
legend('x', 'y', 'z');


subplot(2, 2, 2);
plot(t, velocity(:, 1), 'r', 'LineWidth', 1.5); hold on;
plot(t, velocity(:, 2), 'g', 'LineWidth', 1.5);
plot(t, velocity(:, 3), 'b', 'LineWidth', 1.5);
title('Velocity (vx, vy, vz)');
xlabel('Time (s)');
ylabel('Velocity (m/s)');
legend('vx', 'vy', 'vz');

subplot(2, 2, 3);
plot(t, orientation(:, 1), 'r', 'LineWidth', 1.5); hold on;
plot(t, orientation(:, 2), 'g', 'LineWidth', 1.5);
plot(t, orientation(:, 3), 'b', 'LineWidth', 1.5);
title('Orientation (Roll, Pitch, Yaw)');
xlabel('Time (s)');
ylabel('Angle (rad)');
legend('Roll', 'Pitch', 'Yaw');


subplot(2, 2, 4);
plot(t, angular_velocity(:, 1), 'r', 'LineWidth', 1.5); hold on;
plot(t, angular_velocity(:, 2), 'g', 'LineWidth', 1.5);
plot(t, angular_velocity(:, 3), 'b', 'LineWidth', 1.5);
title('Angular Velocity (p, q, r)');
xlabel('Time (s)');
ylabel('Angular Velocity (rad/s)');
legend('p', 'q', 'r');

sgtitle('Estimated Quadrotor Dynamics: Kalman Filter');

figure;


subplot(6, 2, 1);
plot(t, position_actual(:, 1) - position(:, 1), 'r', 'LineWidth', 1.5);
title('Position Error: X');

ylabel('Error (m)');

subplot(6, 2, 2);
plot(t, position_actual(:, 2) - position(:, 2), 'r', 'LineWidth', 1.5);
title('Position Error: Y');

ylabel('Error (m)');

subplot(6, 2, 3);
plot(t, position_actual(:, 3) - position(:, 3), 'r', 'LineWidth', 1.5);
title('Position Error: Z');
xlabel('Time (s)');
ylabel('Error (m)');


subplot(6, 2, 4);
plot(t, velocity_actual(:, 1) - velocity(:, 1), 'r', 'LineWidth', 1.5);
title('Velocity Error: X');

ylabel('Error (m/s)');

subplot(6, 2, 5);
plot(t, velocity_actual(:, 2) - velocity(:, 2), 'r', 'LineWidth', 1.5);
title('Velocity Error: Y');

ylabel('Error (m/s)');

subplot(6, 2, 6);
plot(t, velocity_actual(:, 3) - velocity(:, 3), 'r', 'LineWidth', 1.5);
title('Velocity Error: Z');

ylabel('Error (m/s)');

subplot(6, 2, 7);
plot(t, orientation_actual(:, 1) - orientation(:, 1), 'r', 'LineWidth', 1.5);
title('Orientation Error: Roll');

ylabel('Error (rad)');

subplot(6, 2, 8);
plot(t, orientation_actual(:, 2) - orientation(:, 2), 'r', 'LineWidth', 1.5);
title('Orientation Error: Pitch');
ylabel('Error (rad)');

subplot(6, 2, 9);
plot(t, orientation_actual(:, 3) - orientation(:, 3), 'r', 'LineWidth', 1.5);
title('Orientation Error: Yaw');

ylabel('Error (rad)');


subplot(6, 2, 10);
plot(t, angular_velocity_actual(:, 1) - angular_velocity(:, 1), 'r', 'LineWidth', 1.5);
title('Angular Velocity Error: p');

ylabel('Error (rad/s)');

subplot(6, 2, 11);
plot(t, angular_velocity_actual(:, 2) - angular_velocity(:, 2), 'r', 'LineWidth', 1.5);
title('Angular Velocity Error: q');
xlabel('Time (s)');
ylabel('Error (rad/s)');

subplot(6, 2, 12);
plot(t, angular_velocity_actual(:, 3) - angular_velocity(:, 3), 'r', 'LineWidth', 1.5);
title('Angular Velocity Error: r');
xlabel('Time (s)');
ylabel('Error (rad/s)');

sgtitle('State Errors: Kalman Filter vs Actual');

%% Moving Horizon Estimator
x_initial_estimate = zeros(12,1);

end_sim = length(t_span);

x_mhe_estimator = zeros(end_sim,12);
tic;
x_mhe_estimator(1,:) = transpose(x_initial_estimate);
N = 20;
for k=2:end_sim

    if k<=N
        y_measure = y_ld(1:k,:);
        u_estimator = delta_input(1:k,:);
    else
        y_measure = y_ld(k-N+1:k,:);
        u_estimator = delta_input(k-N+1:k,:);

    end

    [omega,psi,z,W_vn,W_dn] = matrix_obt(y_measure,u_estimator,ld.C,ld.A,ld.B,R,state_noise_covar);

    sol = LSE(psi,z,omega,W_vn,W_dn);

    state = sol(1:12);

    l = length(sol)/12;

    for i = 2:l

        state = linear_model(ld.A,ld.B,sol( (i-1)*12+1:i*12 ),state,transpose(u_estimator(i-1,:)) );

    end
    x_mhe_estimator(k,:) = transpose(state);
end

elapsed_time = toc;

fprintf('Time Required for Moving Horizon Estimator: %.4f seconds\n', elapsed_time)
position = x_mhe_estimator(:, 1:3);  
velocity = x_mhe_estimator(:, 4:6); 
orientation = x_mhe_estimator(:, 7:9);  
angular_velocity = x_mhe_estimator(:, 10:12);  

position_actual = x_ld(:, 1:3);  
velocity_actual = x_ld(:, 4:6); 
orientation_actual = x_ld(:, 7:9);  
angular_velocity_actual = x_ld(:, 10:12);  
figure;


subplot(2, 2, 1);
plot(t, position(:, 1), 'r', 'LineWidth', 1.5); hold on;
plot(t, position(:, 2), 'g', 'LineWidth', 1.5);
plot(t, position(:, 3), 'b', 'LineWidth', 1.5);
title('Position (x, y, z)');
xlabel('Time (s)');
ylabel('Position (m)');
legend('x', 'y', 'z');


subplot(2, 2, 2);
plot(t, velocity(:, 1), 'r', 'LineWidth', 1.5); hold on;
plot(t, velocity(:, 2), 'g', 'LineWidth', 1.5);
plot(t, velocity(:, 3), 'b', 'LineWidth', 1.5);
title('Velocity (vx, vy, vz)');
xlabel('Time (s)');
ylabel('Velocity (m/s)');
legend('vx', 'vy', 'vz');


subplot(2, 2, 3);
plot(t, orientation(:, 1), 'r', 'LineWidth', 1.5); hold on;
plot(t, orientation(:, 2), 'g', 'LineWidth', 1.5);
plot(t, orientation(:, 3), 'b', 'LineWidth', 1.5);
title('Orientation (Roll, Pitch, Yaw)');
xlabel('Time (s)');
ylabel('Angle (rad)');
legend('Roll', 'Pitch', 'Yaw');


subplot(2, 2, 4);
plot(t, angular_velocity(:, 1), 'r', 'LineWidth', 1.5); hold on;
plot(t, angular_velocity(:, 2), 'g', 'LineWidth', 1.5);
plot(t, angular_velocity(:, 3), 'b', 'LineWidth', 1.5);
title('Angular Velocity (p, q, r)');
xlabel('Time (s)');
ylabel('Angular Velocity (rad/s)');
legend('p', 'q', 'r');

sgtitle('Estimated Quadrotor Dynamics: MHE Filter');

figure;

subplot(6, 2, 1);
plot(t, position_actual(:, 1) - position(:, 1), 'r', 'LineWidth', 1.5);
title('Position Error: X');

ylabel('Error (m)');

subplot(6, 2, 2);
plot(t, position_actual(:, 2) - position(:, 2), 'r', 'LineWidth', 1.5);
title('Position Error: Y');

ylabel('Error (m)');

subplot(6, 2, 3);
plot(t, position_actual(:, 3) - position(:, 3), 'r', 'LineWidth', 1.5);
title('Position Error: Z');
xlabel('Time (s)');
ylabel('Error (m)');


subplot(6, 2, 4);
plot(t, velocity_actual(:, 1) - velocity(:, 1), 'r', 'LineWidth', 1.5);
title('Velocity Error: X');

ylabel('Error (m/s)');

subplot(6, 2, 5);
plot(t, velocity_actual(:, 2) - velocity(:, 2), 'r', 'LineWidth', 1.5);
title('Velocity Error: Y');

ylabel('Error (m/s)');

subplot(6, 2, 6);
plot(t, velocity_actual(:, 3) - velocity(:, 3), 'r', 'LineWidth', 1.5);
title('Velocity Error: Z');

ylabel('Error (m/s)');


subplot(6, 2, 7);
plot(t, orientation_actual(:, 1) - orientation(:, 1), 'r', 'LineWidth', 1.5);
title('Orientation Error: Roll');

ylabel('Error (rad)');

subplot(6, 2, 8);
plot(t, orientation_actual(:, 2) - orientation(:, 2), 'r', 'LineWidth', 1.5);
title('Orientation Error: Pitch');
ylabel('Error (rad)');

subplot(6, 2, 9);
plot(t, orientation_actual(:, 3) - orientation(:, 3), 'r', 'LineWidth', 1.5);
title('Orientation Error: Yaw');

ylabel('Error (rad)');

subplot(6, 2, 10);
plot(t, angular_velocity_actual(:, 1) - angular_velocity(:, 1), 'r', 'LineWidth', 1.5);
title('Angular Velocity Error: p');

ylabel('Error (rad/s)');

subplot(6, 2, 11);
plot(t, angular_velocity_actual(:, 2) - angular_velocity(:, 2), 'r', 'LineWidth', 1.5);
title('Angular Velocity Error: q');
xlabel('Time (s)');
ylabel('Error (rad/s)');

subplot(6, 2, 12);
plot(t, angular_velocity_actual(:, 3) - angular_velocity(:, 3), 'r', 'LineWidth', 1.5);
title('Angular Velocity Error: r');
xlabel('Time (s)');
ylabel('Error (rad/s)');

sgtitle('State Errors: MHE Filter vs Actual');

%% Unscented_Kalman Filter
x_initial_estimate = zeros(12,1);
x_initial_estimate(1:3,1) = 1;
end_sim  = length(t_span);
P_initial = 1000*eye(12);
x_ukalman_estimate = zeros(end_sim,12);
x_ukalman_estimate(1, :) = x_initial_estimate;
P = P_initial;

tic;
for index = 1:end_sim-1
    u =  transpose(u_input(i,:));
    [x_est, P_update] = Unscented_Kalman( transpose(y_nd(index + 1, :)), R, ld.C, P, transpose(x_ukalman_estimate(index, :)), u, state_noise_covar,system_params,t_step);
    x_ukalman_estimate(index + 1, :) = x_est;
    P = P_update;
end

elapsed_time = toc;

fprintf('Time Required for Unscented Kalman Filter: %.4f seconds\n', elapsed_time)

position = x_ukalman_estimate(:, 1:3); 
velocity = x_ukalman_estimate(:, 4:6); 
orientation = x_ukalman_estimate(:, 7:9);  
angular_velocity = x_ukalman_estimate(:, 10:12);  

position_actual = x_nd(:, 1:3);  
velocity_actual = x_nd(:, 4:6); 
orientation_actual = x_nd(:, 7:9);  
angular_velocity_actual = x_nd(:, 10:12);  
figure;


subplot(2, 2, 1);
plot(t, position(:, 1), 'r', 'LineWidth', 1.5); hold on;
plot(t, position(:, 2), 'g', 'LineWidth', 1.5);
plot(t, position(:, 3), 'b', 'LineWidth', 1.5);
title('Position (x, y, z)');
xlabel('Time (s)');
ylabel('Position (m)');
legend('x', 'y', 'z');


subplot(2, 2, 2);
plot(t, velocity(:, 1), 'r', 'LineWidth', 1.5); hold on;
plot(t, velocity(:, 2), 'g', 'LineWidth', 1.5);
plot(t, velocity(:, 3), 'b', 'LineWidth', 1.5);
title('Velocity (vx, vy, vz)');
xlabel('Time (s)');
ylabel('Velocity (m/s)');
legend('vx', 'vy', 'vz');

subplot(2, 2, 3);
plot(t, orientation(:, 1), 'r', 'LineWidth', 1.5); hold on;
plot(t, orientation(:, 2), 'g', 'LineWidth', 1.5);
plot(t, orientation(:, 3), 'b', 'LineWidth', 1.5);
title('Orientation (Roll, Pitch, Yaw)');
xlabel('Time (s)');
ylabel('Angle (rad)');
legend('Roll', 'Pitch', 'Yaw');

subplot(2, 2, 4);
plot(t, angular_velocity(:, 1), 'r', 'LineWidth', 1.5); hold on;
plot(t, angular_velocity(:, 2), 'g', 'LineWidth', 1.5);
plot(t, angular_velocity(:, 3), 'b', 'LineWidth', 1.5);
title('Angular Velocity (p, q, r)');
xlabel('Time (s)');
ylabel('Angular Velocity (rad/s)');
legend('p', 'q', 'r');

sgtitle('Estimated Quadrotor Dynamics: Unscented Kalman Filter');

figure;

subplot(6, 2, 1);
plot(t, position_actual(:, 1) - position(:, 1), 'r', 'LineWidth', 1.5);
title('Position Error: X');

ylabel('Error (m)');

subplot(6, 2, 2);
plot(t, position_actual(:, 2) - position(:, 2), 'r', 'LineWidth', 1.5);
title('Position Error: Y');

ylabel('Error (m)');

subplot(6, 2, 3);
plot(t, position_actual(:, 3) - position(:, 3), 'r', 'LineWidth', 1.5);
title('Position Error: Z');
xlabel('Time (s)');
ylabel('Error (m)');

subplot(6, 2, 4);
plot(t, velocity_actual(:, 1) - velocity(:, 1), 'r', 'LineWidth', 1.5);
title('Velocity Error: X');

ylabel('Error (m/s)');

subplot(6, 2, 5);
plot(t, velocity_actual(:, 2) - velocity(:, 2), 'r', 'LineWidth', 1.5);
title('Velocity Error: Y');

ylabel('Error (m/s)');

subplot(6, 2, 6);
plot(t, velocity_actual(:, 3) - velocity(:, 3), 'r', 'LineWidth', 1.5);
title('Velocity Error: Z');

ylabel('Error (m/s)');


subplot(6, 2, 7);
plot(t, orientation_actual(:, 1) - orientation(:, 1), 'r', 'LineWidth', 1.5);
title('Orientation Error: Roll');

ylabel('Error (rad)');

subplot(6, 2, 8);
plot(t, orientation_actual(:, 2) - orientation(:, 2), 'r', 'LineWidth', 1.5);
title('Orientation Error: Pitch');
ylabel('Error (rad)');

subplot(6, 2, 9);
plot(t, orientation_actual(:, 3) - orientation(:, 3), 'r', 'LineWidth', 1.5);
title('Orientation Error: Yaw');

ylabel('Error (rad)');


subplot(6, 2, 10);
plot(t, angular_velocity_actual(:, 1) - angular_velocity(:, 1), 'r', 'LineWidth', 1.5);
title('Angular Velocity Error: p');

ylabel('Error (rad/s)');

subplot(6, 2, 11);
plot(t, angular_velocity_actual(:, 2) - angular_velocity(:, 2), 'r', 'LineWidth', 1.5);
title('Angular Velocity Error: q');
xlabel('Time (s)');
ylabel('Error (rad/s)');

subplot(6, 2, 12);
plot(t, angular_velocity_actual(:, 3) - angular_velocity(:, 3), 'r', 'LineWidth', 1.5);
title('Angular Velocity Error: r');
xlabel('Time (s)');
ylabel('Error (rad/s)');

sgtitle('State Errors: Unscented Kalman Filter vs Actual');






function dx_dt = quadrotor_nonlinear(t,x, u_in , t_input,state_noise,params)

    b_torque = params.b_torque;
    L = params.L;
    k_motor = params.k_motor;
    m = params.m;
    g = params.g;
    I = params.I;
    kd = params.kd;

    u = interp1(t_input,u_in,t,'linear', 'extrap');

    tau_phi = k_motor * L * (u(2) - u(4));
    tau_theta = k_motor * L * (u(1) - u(3));
    tau_psi = b_torque*(u(1) - u(2)  + u(3) - u(4));
    T_thrust = k_motor * (u(1) + u(2) + u(3) + u(4));
 



    x1 =(x(1:3));    
    x2 = (x(4:6));    
    x3 = (x(7:9));   
    x4 = (x(10:12));  


    phi_roll = x3(1); theta_pitch = x3(2); psi_yaw = x3(3);
   
    Rb = [1 sin(phi_roll)*tan(theta_pitch) cos(phi_roll)*tan(theta_pitch); 0 cos(phi_roll) -sin(phi_roll); 0 sin(phi_roll)/cos(theta_pitch) cos(phi_roll)*cos(theta_pitch)];
    R = rotation_matrix(phi_roll, theta_pitch, psi_yaw);

    
    gravity = [0; 0; -g];
    thrust_body_frame = [0; 0; T_thrust];
    thrust_inertial_frame = (1/m) * R * thrust_body_frame;


    tau = [tau_phi; tau_theta; tau_psi];
    omega_dot = inv(I) * (tau - cross(x4, I * x4));

    %Fd = -kd * x2;
    noise  = transpose(mvnrnd(zeros(12,1), state_noise));
 
    dx_dt1 = zeros(12, 1);
    dx_dt1(1:3) = x2;
    dx_dt1(4:6) = gravity + thrust_inertial_frame; %+ Fd;
    dx_dt1(7:9) = inv(Rb) * x4;
    dx_dt1(10:12) = omega_dot;

    dx_dt = dx_dt1 + noise;

end

function R = rotation_matrix(phi_roll, theta_pitch, psi_yaw)

    R_x = [1, 0, 0;
           0, cos(phi_roll), -sin(phi_roll);
           0, sin(phi_roll), cos(phi_roll)];
       
    R_y = [cos(theta_pitch), 0, sin(theta_pitch);
           0, 1, 0;
           -sin(theta_pitch), 0, cos(theta_pitch)];
       
    R_z = [cos(psi_yaw), -sin(psi_yaw), 0;
           sin(psi_yaw), cos(psi_yaw), 0;
           0, 0, 1];
       
   
    R = R_z * R_y * R_x;
end



function [system_lc,system_ld, Noise] = linearize_quadrotor(params,noise_covar, sampling_time)

   b_torque = params.b_torque;
    L = params.L;
    k_motor = params.k_motor;
    m = params.m;
    g = params.g;
    I = params.I;
    kd = params.kd;

    syms  x y z v_x v_y v_z phi_roll theta_pitch psi_yaw omega_x omega_y omega_z u1 u2 u3 u4

R_x = [1, 0, 0;
           0, cos(phi_roll), -sin(phi_roll);
           0, sin(phi_roll), cos(phi_roll)];
R_y = [cos(theta_pitch), 0, sin(theta_pitch);
           0, 1, 0;
           -sin(theta_pitch), 0, cos(theta_pitch)];
       
R_z = [cos(psi_yaw), -sin(psi_yaw), 0;
           sin(psi_yaw), cos(psi_yaw), 0;
           0, 0, 1];

R = R_z*R_y*R_x;

thrust = k_motor*(u1 + u2 + u3 + u4);

tau_phi = L*k_motor*(u2 - u4);

tau_theta = L*k_motor*(u1- u3);

tau_psi = b_torque*(u1 - u2  + u3 - u4);

force = -[0;0; m*g] + R*[0;0;thrust] ;


Rb = [1 sin(phi_roll)*tan(theta_pitch) cos(phi_roll)*tan(theta_pitch); 0 cos(phi_roll) -sin(phi_roll); 0 sin(phi_roll)/cos(theta_pitch) cos(phi_roll)*cos(theta_pitch)];
eta_dot  = Rb*[omega_x; omega_y; omega_z];
omega = [omega_x; omega_y; omega_z];

omega_dot = inv(I)*([tau_phi;tau_theta;tau_psi] - cross(omega,I*omega));

x_dot = v_x;
y_dot = v_y;
z_dot = v_z;

a_x = force(1)/m;
a_y = force(2)/m;
a_z = force(3)/m;


phi_dot = eta_dot(1);
theta_dot = eta_dot(2);
psi_dot = eta_dot(3);


omega_x_dot = omega_dot(1);
omega_y_dot = omega_dot(2);
omega_z_dot = omega_dot(3);

jacobian_matrix = jacobian([x_dot,y_dot,z_dot,a_x,a_y, a_z,phi_dot, theta_dot, psi_dot, omega_x_dot, omega_y_dot, omega_z_dot],[x y z v_x v_y v_z phi_roll theta_pitch psi_yaw omega_x omega_y omega_z u1 u2 u3 u4]);

operating_point_input = ((m*g)/(4*k_motor))*ones(4,1);
operating_point_states = zeros(12,1);

operating_point_states(1) = 1;

operating_point_states(2) = 1;

operating_point_states(3) = 1;

x = operating_point_states(1);
y = operating_point_states(2);
z = operating_point_states(3);


v_x = operating_point_states(4);
v_y = operating_point_states(5);
v_z = operating_point_states(6);

phi_roll = operating_point_states(7);
theta_pitch = operating_point_states(8);
psi_yaw = operating_point_states(9);

omega_x = operating_point_states(10);
omega_y = operating_point_states(11);
omega_z = operating_point_states(12);

u1 = operating_point_input(1);
u2 = operating_point_input(2);
u3 = operating_point_input(3);
u4 = operating_point_input(4);

jacobian_subs = double(subs(jacobian_matrix));

A_matrix = jacobian_subs(:,1:12);
B_matrix = jacobian_subs(:,13:16);


C1 = zeros(6,12);
C1(1,1) = 1;
C1(2,2) = 1;
C1(3,3) = 1;

C1(4,7) = 1;
C1(5,8) = 1;
C1(6,9) = 1;
system_lc  = ss(A_matrix,B_matrix,C1,zeros(6,4));

system_ld = c2d(system_lc,sampling_time);

Ad = system_ld.A;
Bd = system_ld.B;

controlability_matrix = ctrb(Ad,Bd);

disp("The rank of controllability matrix is")
disp(rank(controlability_matrix))

if rank(controlability_matrix)==12
  disp("The system is fully controllable")

else
    disp("The system is not controlable")
end

obsv_mat = obsv(Ad,C1);

disp("The rank of observability matrix is")
disp(rank(controlability_matrix))

if rank(obsv_mat)==12
  disp("The system is fully observable")

else
    disp("The system is not observable")
end

num_steps = 1000;  
dt = sampling_time / num_steps;  

Noise = zeros(size(Ad)); 

for tau = 0:dt:sampling_time
    exp_A_tau = expm(A_matrix * tau); 

    Noise = Noise + exp_A_tau * noise_covar * exp_A_tau' * dt;  
end

end




function dx_dt = quadrotor_nd(x, u,params)

    b_torque = params.b_torque;
    L = params.L;
    k_motor = params.k_motor;
    m = params.m;
    g = params.g;
    I = params.I;
    kd = params.kd;

    tau_phi = k_motor * L * (u(2) - u(4));
    tau_theta = k_motor * L * (u(1) - u(3));
    tau_psi = b_torque*(u(1) - u(2)  + u(3) - u(4));
    T_thrust = k_motor * (u(1) + u(2) + u(3) + u(4));
 


   
    x1 =(x(1:3));    
    x2 = (x(4:6));   
    x3 = (x(7:9));    
    x4 = (x(10:12));  

 
    phi_roll = x3(1); theta_pitch = x3(2); psi_yaw = x3(3);
    Rb = [1 sin(phi_roll)*tan(theta_pitch) cos(phi_roll)*tan(theta_pitch); 0 cos(phi_roll) -sin(phi_roll); 0 sin(phi_roll)/cos(theta_pitch) cos(phi_roll)*cos(theta_pitch)];
    R = rotation_matrix(phi_roll, theta_pitch, psi_yaw);

    
    gravity = [0; 0; -g];
    thrust_body_frame = [0; 0; T_thrust];
    thrust_inertial_frame = (1/m) * R * thrust_body_frame;

    tau = [tau_phi; tau_theta; tau_psi];

    omega_dot = inv(I) * (tau - cross(x4, I * x4));
    dx_dt1 = zeros(12, 1);
    dx_dt1(1:3) = x2;
    dx_dt1(4:6) = gravity + thrust_inertial_frame; %+ Fd;
    dx_dt1(7:9) = Rb * x4;
    dx_dt1(10:12) = omega_dot;

    dx_dt = dx_dt1;

end


function [x_est,P] = Kalman(y_obs,R,A,B,C,Pin,x_prev,u_in,BdQdBd)

P_k_k_prev = (A*Pin*transpose(A)) + BdQdBd;


x_tilda = (A*x_prev) + (B*u_in) ;


y_hat = C*x_tilda;

L_gain = P_k_k_prev*transpose(C)*inv((C*P_k_k_prev*transpose(C)) + R);


P = (eye(12) - L_gain*C)*P_k_k_prev;

x_est = x_tilda + L_gain*(y_obs - y_hat);
end


function [x_est,P] = Unscented_Kalman(y_obs,R,C,Pin,x_prev,u_in,BdQdBd,param,t_step)

n = 12;
alpha = 1e-3;
kappa = 3-n;
beta = 2;
lambda = ((alpha^2)*(n+kappa)) - n;

sqrt_P = sqrtm((n+lambda)*Pin);

sigma_points = zeros(2*n +1 ,n);

sigma_points(1,:) = transpose(x_prev);

for index=2:n+1
    sigma_points(index,:) = transpose(x_prev  + sqrt_P(:,index-1));
end

for index=n+2:2*n+1
    sigma_points(index,:) = transpose(x_prev  - sqrt_P(:,index-1-n));
end

estimates_model = zeros(2*n+1,n);

for index = 1:2*n+1

    estimates_model(index,:)  = quadrotor_update(t_step,transpose(sigma_points(index,:)),u_in,param);
end

weights_m = zeros(2*n + 1, 1);
weights_c = zeros(2*n + 1, 1);
weights_m(1) = lambda/(n + lambda);
weights_c(1) = lambda/(n + lambda) + (1 - alpha^2 + beta);
weights_m(2:end) = 1/(2*(n + lambda));
weights_c(2:end) = 1/(2*(n + lambda));

x_k_k_1 =transpose(estimates_model)*weights_m;

P_k_k_1 = zeros(n,n);

for index = 1:2*n+1

    P_k_k_1 = P_k_k_1 +  weights_c(index,1) * ( transpose(estimates_model(index,:)) - x_k_k_1)*transpose(transpose(estimates_model(index,:)) - x_k_k_1)  ;

end

P_k_k_1 = P_k_k_1 + BdQdBd;

sigma_points_obsv = zeros(2*n + 1,n);
sqrt_P_k = sqrtm((n+lambda)*P_k_k_1);

sigma_points_obsv(1,:) = x_k_k_1;

for index=2:n+1
    sigma_points_obsv(index,:) = transpose(x_k_k_1  + sqrt_P_k(:,index-1));
end

for index=n+2:2*n+1
    sigma_points_obsv(index,:) = transpose(x_k_k_1  - sqrt_P_k(:,index-1-n));
end

estimates_obsv = zeros(2*n +1 ,length(y_obs));


for index=1:2*n+1
    estimates_obsv(index,:) = transpose(C*transpose(sigma_points_obsv(index,:)));
end

y_k_k_1 = transpose(estimates_obsv)*weights_m;

P_k_y = zeros(length(y_obs),length(y_obs));

for index = 1:2*n+1

    P_k_y = P_k_y + weights_c(index,1) * (transpose(estimates_obsv(index,:)) - y_k_k_1)*transpose(transpose(estimates_obsv(index,:)) - y_k_k_1);

end
P_k_y = P_k_y + R;

P_k_xy = zeros(n,length(y_obs));

for index = 1:2*n+1
    weight = weights_c(index);

    matrix = (transpose(sigma_points(index,:)) - x_k_k_1)*transpose(transpose( estimates_obsv(index,:) ) - y_k_k_1);
    P_k_xy = P_k_xy + weight *matrix ;

end

L_gain = P_k_xy/(P_k_y);

x_est = x_k_k_1 + L_gain*(y_obs - y_k_k_1);
I= eye(12);

P =(I - L_gain * C) * P_k_k_1 * transpose(I - L_gain * C) + L_gain * R * transpose(L_gain);

end


function [x_next] = quadrotor_update(dt, x, u_in,params)

    dx_dt = quadrotor_nd(x, u_in, params);
    
    % Update the state using Euler's method
    x_next = x + dt * dx_dt ;
end


function [state] = linear_model(A,B,d,x,u)

state = A*x + B*u + d;

end

function[omega,psi,Z,W_vn,W_dn] = matrix_obt(y,u,C,A,B,R,Q)


N = size(y,1);
n = size(A,1);
m = size(C,1);
omega_matrix = zeros(m,n*(N-1));
Z = transpose(y(1,:));


psi = C;

inp_term = zeros(12,1);

row = omega_matrix;

sensor = inv(R);
W_vn = sensor;
system = inv(Q);

for index = 2:N

    for j =1:N-1
        if j<index
            if j==1
                col = C*(A^(index-j-1));
            else
                col = horzcat(col,C*(A^(index-j-1)));
            end
        else
            col = horzcat(col,zeros(m,n));
          
        end
    end
    row = vertcat(row,col);
end

omega = row;

for i=2:N

    psi = vertcat(psi,C*(A^(i-1)));

    W_vn = blkdiag(W_vn,sensor);
    if i==2
        W_dn = system;
    else
        W_dn = blkdiag(W_dn,system);
    end
    input = B*transpose(u(i-1,:));

    inp_term = A*inp_term  + input;

    Z = vertcat(Z, ( transpose(y(i,:)) - C*inp_term) );


end

end


function [sol] = LSE(psi,z,omega,W_vn,W_dn)

a11 = transpose(psi)*W_vn*psi;

a12 = transpose(psi)*W_vn*omega;

a21 = transpose(omega)*W_vn*psi;


a22 = (transpose(omega)*W_vn*omega) + W_dn;


b11 = transpose(psi)*W_vn;
b12 = transpose(omega)*W_vn;

B = vertcat(b11,b12) * z;

a1 = horzcat(a11,a12);
a2 = horzcat(a21,a22);

A = vertcat(a1,a2);

sol = A\B;

end

