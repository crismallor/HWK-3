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

%% Compute the two-impulse maneuver for the hop using the simplified CW equations
omega = n;
Deltay = rf(2)-r0(2); % km
[DV1,DV2]=hopping(omega,Deltay);

%% Plotting different trajectories
% 1. Half orbital period

t_vector1 = linspace(0, tau/2, 1000); % time vector (s)

for i = 1:length(t_vector1)
    t = t_vector1(i); % time (s)
    PHI=CW_TransMat(omega,t); % compute the CW matrices
    Phi_rr = PHI(1:3,1:3); % position-to-position state transition matrix
    Phi_rv = PHI(1:3,4:6); % position-to-velocity state transition matrix
    r(:,i) = Phi_rr*r0 + Phi_rv*v0; % compute the relative position (km)
end

figure;
plot(r(1,:),r(2,:),'b','LineWidth',2); % plot the trajectory
xlabel('x (R-bar direction) [km]');
ylabel('y (V-bar direction) [km]');
title('Relative Trajectory of the Spacecraft');
grid on;
axis equal;
set(gca, 'FontSize', 12, 'FontName', 'Arial');
legend('Trajectory', 'Location', 'best');