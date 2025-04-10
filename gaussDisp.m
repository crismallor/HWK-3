function [x_in, y_in] = gaussDisp(sigma_x, sigma_y, N)

% Generate N random points inside the ellipse
x_in = [];
y_in = [];

while length(x_in) < N
    % Generate random points within the bounds of the ellipse
    x0 = -sigma_x + 2 * sigma_x * rand; % Uniform distribution in [-sigma_x, sigma_x]
    y0 = -sigma_y + 2 * sigma_y * rand; % Uniform distribution in [-sigma_y, sigma_y]
    
    % Check if the point lies within the ellipse
    if (x0^2 / sigma_x^2 + y0^2 / sigma_y^2) <= 1
        x_in = [x_in, x0];
        y_in = [y_in, y0];
    end
end

figure; 
% Display the points inside the ellipse
scatter(x_in, y_in, 'filled');
xlabel('x');
ylabel('y');
title('Points within the ellipse');
axis equal;
grid on;

% Plot the ellipse boundary
theta = linspace(0, 2*pi, 100);
ellipse_x = sigma_x * cos(theta);
ellipse_y = sigma_y * sin(theta);
hold on;
plot(ellipse_x, ellipse_y, 'r', 'LineWidth', 2);
legend('Points', 'Ellipse');
hold off;



end 