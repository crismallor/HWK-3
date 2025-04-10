clear; clc; close all;

%% Data
mu = 398600; % Earth gravitational parameter (km^3/s^2)
h_t =  390; % Tiangong space station height (km)
RE = 6371; % Earth radius (km)
rh = h_t + RE; % Space station distance from the center of the Earth (km)

% Initial conditions
% relative position
r0 = [0; -1; 0]; % km
% relative velocity
v0 = [0; 0; 0]; % km/s

% Final condition
% relative position
rf = [0; -0.1; 0]; % km

% Time of maneuver
n = sqrt(mu/rh^3); % mean angular rate rad/s
tau = 2*pi/n; % Tiangong's orbital period (s)
t = tau/2; % time of maneuver (s)

% t = 1.49*3600;
% n = 0.0011569;

%% a.) Compute the two-impulse maneuver for the hop using the simplified CW equations
omega = n;
Deltay = rf(2)-r0(2); % km
[DV1,DV2]=hopping(omega, Deltay);

%% b.) Plotting different trajectories
% 1. Half orbital period

t_vector1 = linspace(0, tau/2, 1000); % time vector (s)
t_vector2 = linspace(0, 2*tau, 2000); % time vector 2 (s) 

for i = 1:length(t_vector1)
    t = t_vector1(i); % time (s)
    PHI=CW_TransMat(omega,t); % compute the CW matrices
    Phi_rr = PHI(1:3,1:3); % position-to-position state transition matrix
    Phi_rv = PHI(1:3,4:6); % position-to-velocity state transition matrix
    r1(:,i) = Phi_rr*r0 + Phi_rv*DV1/1000; % compute the relative position (km)
end

for i = 1:length(t_vector2)
    t = t_vector2(i); % time (s)
    PHI=CW_TransMat(omega,t); % compute the CW matrices
    Phi_rr = PHI(1:3,1:3); % position-to-position state transition matrix
    Phi_rv = PHI(1:3,4:6); % position-to-velocity state transition matrix
    r2(:,i) = Phi_rr*r0 + Phi_rv*DV1/1000; % compute the relative position (km)
end

figure;
plot(r1(1,:),r1(2,:),'b','LineWidth',2); % plot the trajectory
xlabel('x (R-bar direction) [km]');
ylabel('y (V-bar direction) [km]');
title('Relative Trajectory of the Spacecraft');
grid on;
axis equal;
set(gca, 'FontSize', 12, 'FontName', 'Arial');
legend('Trajectory', 'Location', 'best');

figure;
plot(r2(1,:),r2(2,:),'b','LineWidth',2); % plot the trajectory
xlabel('x (R-bar direction) [km]');
ylabel('y (V-bar direction) [km]');
title('Relative Trajectory of the Spacecraft for two orbital periods');
grid on;
axis equal;
set(gca, 'FontSize', 12, 'FontName', 'Arial');
legend('Trajectory', 'Location', 'best');

%% c.) Montecarlo Anaylisis: 
    
%Limits for the deviation in the ellipse: parameters sigma_x and sigmay
%define the limit of the ellipse in which random values will be taken. 
sigma_x = 10; 
sigma_y = 20; 
N = 100; %Number of samples

[x0, y0] = gaussDisp(sigma_x, sigma_y, N); %Computation of the random points, x0 and y0 are given in meters 
r0_vector = r0 + [x0; y0; zeros(1, length(x0))]/1000; %Dividing by 1000, as the vector is given in meters

%% d.) Compute and plot trajectories for all initial positions
figure;
hold on;
for j = 1:N
    r0_sample = r0_vector(:, j); % Select a random initial position
    for i = 1:length(t_vector1)
        t = t_vector1(i); % time (s)
        PHI = CW_TransMat(omega, t); % compute the CW matrices
        Phi_rr = PHI(1:3, 1:3); % position-to-position state transition matrix
        Phi_rv = PHI(1:3, 4:6); % position-to-velocity state transition matrix
        r_sample(:, i) = Phi_rr * r0_sample + Phi_rv * DV1 / 1000; % compute the relative position (km)
    end
    plot(r_sample(1, :), r_sample(2, :)); % plot the trajectory
    
end
xlabel('x (R-bar direction) [km]');
ylabel('y (V-bar direction) [km]');
title('Monte Carlo Analysis: Trajectories for Random Initial Positions');
grid on;
axis equal;
set(gca, 'FontSize', 12, 'FontName', 'Arial');
legend('Trajectories', 'Location', 'best');
hold off;


%% e.) 
% Find the trajectory that ends farthest from the origin at the final time
max_distance = 0;
farthest_trajectory = [];

for j = 1:N
    r0_sample = r0_vector(:, j); % Select a random initial position
    for i = 1:length(t_vector1)
        t = t_vector1(i); % time (s)
        PHI = CW_TransMat(omega, t); % compute the CW matrices
        Phi_rr = PHI(1:3, 1:3); % position-to-position state transition matrix
        Phi_rv = PHI(1:3, 4:6); % position-to-velocity state transition matrix
        r_sample(:, i) = Phi_rr * r0_sample + Phi_rv * DV1 / 1000; % compute the relative position (km)
    end
    % Compute the distance from the final point rf at the final time
    final_distance = norm(r_sample(:, end) - rf);

    if final_distance > max_distance
        max_distance = final_distance;
        r_far = r_sample; % Update the farthest trajectory
        r0_far = r0_sample; % Save the initial conditions for the farthest trajectory
    end

end



% Plot the farthest trajectory
figure;
plot(r_far(1, :), r_far(2, :), 'r', 'LineWidth', 2);
xlabel('x (R-bar direction) [km]');
ylabel('y (V-bar direction) [km]');
title('Farthest Trajectory at Final Time');
grid on;
axis equal;
set(gca, 'FontSize', 12, 'FontName', 'Arial');
legend('Farthest Trajectory', 'Location', 'best');


%Computation of the position at the midspan of the period of the orbit: 
t_midspan = tau/4; % time at midspan (s)
PHI_midspan = CW_TransMat(omega, t_midspan); % compute the CW matrices

Phi_rr_midspan = PHI_midspan(1:3, 1:3); % position-to-position state transition matrix
Phi_rv_midspan = PHI_midspan(1:3, 4:6); % position-to-velocity state transition matrix
Phi_vr_midspan = PHI_midspan(4:6, 1:3); % velocity-to-position state transition matrix
Phi_vv_midspan = PHI_midspan(4:6, 4:6); % velocity-to-velocity state transition matrix

r0_midspan = Phi_rr_midspan * r0_far + Phi_rv_midspan * DV1 / 1000; % compute the relative position (km)
v0_midspan = Phi_vr_midspan * r0_far + Phi_vv_midspan * DV1 / 1000; % compute the relative velocity (km/s)


%Initial Conditions at a quarter of the orbit: 
r0 = r0_midspan; 
v0 = v0_midspan; 

%Final Conditions: 
rf = [0; -0.1; 0]; % km

%Creation of the state vector: 
state0MC = [r0_midspan, v0_midspan]; % state vector (km, km/s)

[DV1_midspan, DV2_midspan] = mid_course(omega, state0MC, rf); 


% Plot the trajectory from 0 to tau/4
t_vector_mid1 = linspace(0, tau/4, 500); % time vector (s)
trajectory_mid1 = zeros(3, length(t_vector_mid1));

for i = 1:length(t_vector_mid1)
    t = t_vector_mid1(i); % time (s)
    PHI = CW_TransMat(omega, t); % compute the CW matrices
    Phi_rr = PHI(1:3, 1:3); % position-to-position state transition matrix
    Phi_rv = PHI(1:3, 4:6); % position-to-velocity state transition matrix
    trajectory_mid1(:, i) = Phi_rr * r0_far + Phi_rv * DV1 / 1000; % compute the relative position (km)
end

% Plot the trajectory from tau/4 to tau/2
t_vector_mid2 = linspace(0, tau/4, 500); % time vector (s)
trajectory_mid2 = zeros(3, length(t_vector_mid2));

for i = 1:length(t_vector_mid2)
    t = t_vector_mid2(i); % time (s)
    PHI = CW_TransMat(omega, t); % compute the CW matrices
    Phi_rr = PHI(1:3, 1:3); % position-to-position state transition matrix
    Phi_rv = PHI(1:3, 4:6); % position-to-velocity state transition matrix
    trajectory_mid2(:, i) = Phi_rr * r0_midspan + Phi_rv * ((DV1_midspan)/1000 + v0_midspan); % compute the relative position (km)
end

% Plot both segments
figure;
hold on;
plot(trajectory_mid1(1, :), trajectory_mid1(2, :), 'b', 'LineWidth', 2); % 0 to tau/4
plot(trajectory_mid2(1, :), trajectory_mid2(2, :), 'r', 'LineWidth', 2); % tau/4 to tau/2
xlabel('x (R-bar direction) [km]');
ylabel('y (V-bar direction) [km]');
title('Trajectory Segments: 0 to tau/4 and tau/4 to tau/2');
grid on;
axis equal;
set(gca, 'FontSize', 12, 'FontName', 'Arial');
legend('0 to tau/4', 'tau/4 to tau/2', 'Location', 'best');
hold off;



