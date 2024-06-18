% FILEPATH: /e:/DATA/TUD/Master/TUD_Master_Y1/Q4/ET4173 Introduction to UWB Technology, Systems and Applications (202223 Q4) - 2252024 - 523 PM/Project/FalseAlarmRate.m
% FalseAlarmRate calculates the false alarm rate based on the estimated trajectory.
% The function assumes that the human movement follows a triangular path, with key points I, C, and A.
% The function generates a trajectory by linearly interpolating between the key points.
% It then calculates the boundary points of a ball with a given width at each position along the trajectory.
% The outer contour of the swept area is calculated using the boundary function.
% The function checks if data points from the estimated trajectory are inside the outer movement area but outside the inner triangle.
% The false alarm rate is calculated as the ratio of the number of points inside the movement area to the total number of data points.
% 
% Inputs:
%   - Estimated_trajectory: A matrix representing the estimated trajectory of the human movement.
%                           Each row represents a data point with x and y coordinates.
%
% Output:
%   - false_alarm_rate: The calculated false alarm rate.
%
% Example usage:
%   estimated_trajectory = [0, 0; 1, 0; -1, 0; 0.5, 0.5; -0.5, -0.5; 2, 2; -1, 1.5; 1, -1.5];
%   rate = FalseAlarmRate(estimated_trajectory);
function num_points_inside = NumberInRegion(Estimated_trajectory)

    % Suppose the human width is 0.37m
    % the shape of the movement is a triangle
    % Path of the human is I - C - A - I
    % Define the key points
    % Define the key points
    % Define the key points
    I = [1, -1];
    C = [1, 1];
    A = [-1, 1];

    num_steps = 100; % Number of steps per segment

    % Human width
    width = 0.37;
    radius = width / 2;

    % Define the order of vertices
    vertices = [I; C; A; I]; % Closed triangle

    % Predefine trajectory arrays
    x_traj = [];
    y_traj = [];

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

    % Define the inner vertices
    I_inner = [0.815, -0.5509];
    C_inner = [0.815, 0.815];
    A_inner = [-0.5534, 0.815];
    inner_vertices = [I_inner; C_inner; A_inner; I_inner];

    % Plot the trajectory and the area covered
    % figure;
    % hold on;
    for k = 1:length(x_traj)
        % Draw the ball at each position
        theta = linspace(0, 2*pi, 100);
        x_circle = radius * cos(theta);
        y_circle = radius * sin(theta);
        % fill(x_traj(k) + x_circle, y_traj(k) + y_circle, 'g', 'FaceAlpha', 0.05, 'EdgeColor', 'none');
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

    % % Plot the inner triangle
    % fill(inner_vertices(:, 1), inner_vertices(:, 2), 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');

    % Define the dataset points (example data points)
    % data_points = [0, 0; 1, 0; -1, 0; 0.5, 0.5; -0.5, -0.5; 2, 2; -1, 1.5; 1, -1.5];
    data_points = Estimated_trajectory;

    % % Plot the data points
    % plot(data_points(:, 1), data_points(:, 2), 'ro');

    % Check if points are inside the outer movement area
    in_outer = inpolygon(data_points(:, 1), data_points(:, 2),boundary_points_x(k_outer), boundary_points_y(k_outer));

    % Check if points are inside the inner triangle
    in_inner = inpolygon(data_points(:, 1), data_points(:, 2), inner_vertices(:, 1), inner_vertices(:, 2));

    % Points inside the outer area but outside the inner triangle
    in_movement_area = in_outer & ~in_inner;

    % % Highlight points inside the movement area
    % plot(data_points(in_movement_area, 1), data_points(in_movement_area, 2), 'bo');

    % Display the number of points inside the movement area
    num_points_inside = sum(in_movement_area);
    % false_alarm_rate =1- num_points_inside / length(data_points);
    % disp(['False Alarm Rate: ', num2str(false_alarm_rate)]);
    % disp(['Number of points inside the movement area: ', num2str(num_points_inside)]);
    % hold off;
    disp(['Number of points inside the movement area: ', num2str(num_points_inside)]);

end