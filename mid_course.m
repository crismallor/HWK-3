function [DVx, DVy] = mid_course(omega,state0MC,rf)

    %Inputs: 
    % omega: mean angular rate (rad/s)
    % state0MC: mid course state vector (position and velocity) of the chaser (km, km/s)
    % rf: final relative position (km)

    %Outputs: 
    % DVx: Delta velocity change in x direction (m/s)
    % DVy: Delta velocity change in y direction (m/s)

     % Time of maneuver
    n = omega; % mean angular rate rad/s
    tau = 2*pi/n; % Tiangong's orbital period (s)
    t = tau/4; % time of maneuver (s)
    
    rf = rf/1000; %Change from meters to kilometers. 

    %1.) 
    r0_midspan = state0MC(1:3, 1); % initial position (km)
    v0_midspan = state0MC(1:3, 2); % initial velocity (km/s)

    % 2. Determination of the CW matrices
    [Phi_rr, Phi_rv, Phi_vr, Phi_vv] = CW_matrices(omega, t); % compute the CW matrices %Duda, que tiempo pongo 0 o tau/4?

    % 3. Initial velocity after the first impulse
    v0_plus = inv(Phi_rv)*(rf-Phi_rr*r0_midspan);

    % 4. Final velocity before the second impulse
    vf_minus = Phi_vr*r0_midspan + Phi_vv*v0_plus;

    % 5. Initial velocity before first impulse
    v0_minus = v0_midspan; % km/s; zeros since both the chaser and target are initially in the same circular orbit

    % 7. Final velocity after the second impulse
    vF_plus = zeros(3,1); % km/s; zeros since both the chaser and target are at the end in the same circular orbit

    % 8. Velocity change for the first impulse
    DV1 = (v0_plus - v0_minus) * 1000; % m/s
    DVx = DV1(1); 
    DVy = DV1(2); 


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