function PHI=CW_TransMat(omega,t)

    % Inputs:
        %   omega - Mean motion (rad/s)
        %   t - Time (s)

        % Outputs:
        %   Phi_rr - Position-to-position state transition matrix
        %   Phi_rv - Position-to-velocity state transition matrix
        %   Phi_vr - Velocity-to-position state transition matrix
        %   Phi_vv - Velocity-to-velocity state transition matrix

        n = omega;

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

        PHI = [Phi_rr Phi_rv; Phi_vr Phi_vv];




end