
function PlotTrueRange()
    % Define vertex coordinates
    I = [1, -1];
    C = [1, 1];
    A = [-1, 1];

    % Define the diameter and radius of the ball
    diameter = 0.37;
    radius = diameter / 2;

    % Define total time and number of steps for the movement
    num_steps = 100; % Number of steps per segment

    % Define the order of vertices
    vertices = [I; C; A; I]; % Closed triangle

    % Predefine trajectory arrays
    x_traj = [];
    y_traj = [];

    % Calculate the trajectory
    for i = 1:3
        % Get the start and end points of the current segment
        start_point = vertices(i, :);
        end_point = vertices(i + 1, :);
        
        % Calculate the trajectory for each segment
        x_segment = linspace(start_point(1), end_point(1), num_steps);
        y_segment = linspace(start_point(2), end_point(2), num_steps);
        
        % Add the current segment's trajectory to the total trajectory
        x_traj = [x_traj, x_segment];
        y_traj = [y_traj, y_segment];
    end

    % Calculate the boundary points of the ball at each position
    boundary_points_x = [];
    boundary_points_y = [];

    for k = 1:length(x_traj)
        theta = linspace(0, 2*pi, 100);
        x_circle = radius * cos(theta);
        y_circle = radius * sin(theta);
        boundary_points_x = [boundary_points_x, x_traj(k) + x_circle];
        boundary_points_y = [boundary_points_y, y_traj(k) + y_circle];
    end

    % Use the boundary function to calculate the outer contour of the swept area
    k_outer = boundary(boundary_points_x', boundary_points_y', 0.03);

    % Inner vertices
    I_inner = [0.815, -0.5509];
    C_inner = [0.815, 0.815];
    A_inner = [-0.5534, 0.815];

    % Plot the result
    % figure;
    hold on;

    % Outer contour
    h1 = plot(boundary_points_x(k_outer), boundary_points_y(k_outer), 'r--', 'LineWidth', 1.5, "DisplayName", "Outer Contour");

    % Inner contour
     plot([I_inner(1) C_inner(1) A_inner(1) I_inner(1)], [I_inner(2) C_inner(2) A_inner(2) I_inner(2)], 'r--', 'LineWidth', 1.5, "DisplayName", "Inner Contour");

    for k = 1:length(x_traj)
        % Draw the ball at each position
        theta = linspace(0, 2*pi, 100);
        x_circle = radius * cos(theta);
        y_circle = radius * sin(theta);
        fill(x_traj(k) + x_circle, y_traj(k) + y_circle, 'g', 'FaceAlpha', 0.05, 'EdgeColor', 'none');
    end

    % Draw the edges of the triangle
    h3 = plot([I(1) C(1) A(1) I(1)], [I(2) C(2) A(2) I(2)], 'k--', 'LineWidth', 1.5, "DisplayName", "True Path");

    % Set axis scale and labels
    axis equal;
    xlabel('X');
    ylabel('Y');
    ylim([-2.5,2.5]);
    xlim([-3,3]);
    title('Outer and Inner Contour of the Swept Area of the Ball');

    % Create legend only for the desired objects
    legend([h1, h3], 'Path Contour', 'True Path');
    % legend([ h3], 'True Path');

    % hold off;
    grid on;

end