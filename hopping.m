function [DV1,DV2]=hopping(omega,Deltay)

    % Inputs:
    % omega - orbit mean angular rate (rad/s); aka 'n'
    % Deltay - change in y position (km)
    
    % Outputs:
    % DV1 - first impulse velocity change (km/s)
    % DV2 - second impulse velocity change (km/s)

    %% Data

    % Relative motion parameters
    % Initial relative position
    r0 = [0; -1; 0]; % km
    % Final relative position
    rf = [0; -0.1; 0]; % km

    % Time of maneuver
    n = omega; % mean angular rate rad/s
    tau = 2*pi/n; % Tiangong's orbital period (s)
    t = tau/2; % time of maneuver (s)

    %t = 5364;
    %% Algorithm
    % 1. Determination of the CW matrices
    [Phi_rr, Phi_rv, Phi_vr, Phi_vv] = CW_matrices(n, t);

    % 2. Initial velocity after the first impulse
    v0_plus = inv(Phi_rv)*(rf-Phi_rr*r0);

    % 3. Final velocity before the second impulse
    vf_minus = Phi_vr*r0 + Phi_vv*v0_plus;

    % 4. Initial velocity before first impulse
    v0_minus = zeros(3,1); % km/s; zeros since both the chaser and target are initially in the same circular orbit

    % 5. Final velocity after the second impulse
    vF_plus = zeros(3,1); % km/s; zeros since both the chaser and target are at the end in the same circular orbit

    % 6. Velocity change for the first impulse
    DV1 = v0_plus - v0_minus; % km/s
    DV1 = norm(DV1); % km/s; magnitude of the velocity change

    % 7. Velocity change for the second impulse
    DV2 = vF_plus - vf_minus; % km/s
    DV2 = norm(DV2); % km/s; magnitude of the velocity change


    function [Phi_rr, Phi_rv, Phi_vr, Phi_vv] = CW_matrices(n, t)
        
        % Inputs:
        %   n - Mean motion (rad/s)
        %   t - Time (s)

        % Outputs:
        %   Phi_rr - Position-to-position state transition matrix
        %   Phi_rv - Position-to-velocity state transition matrix
        %   Phi_vr - Velocity-to-position state transition matrix
        %   Phi_vv - Velocity-to-velocity state transition matrix
    
        % Position-to-position transition matrix (Phi_rr)
        Phi_rr = [4 - 3*cos(n*t), 0, 0;
                  6*(sin(n*t) - n*t), 1, 0;
                  0, 0, cos(n*t)];
    
        % Position-to-velocity transition matrix (Phi_rv)
        Phi_rv = [1/n*sin(n*t), 2/n*(1 - cos(n*t)), 0;
                  2/n*(cos(n*t) - 1), 1/n*(4*sin(n*t) - 3*n*t), 0;
                  0, 0, 1/n*sin(n*t)];
    
        % Velocity-to-position transition matrix (Phi_vr)
        Phi_vr = [3*n*sin(n*t), 0, 0;
                  6*n*(cos(n*t) - 1), 0, 0;
                  0, 0, -n*sin(n*t)];
    
        % Velocity-to-velocity transition matrix (Phi_vv)
        Phi_vv = [cos(n*t), 2*sin(n*t), 0;
                  -2*sin(n*t), 4*cos(n*t) - 3, 0;
                  0, 0, cos(n*t)];
    end
    



end