
function [false_alarm_rate,MeanDistance] = FalseAlarmRate(Estimated_trajectory,ts,JJ)

    % Suppose the human width is 0.37m
    % the shape of the movement is a triangle
    % Path of the human is I - C - A - I
    % Define the key points
    % Define the key points
    % Define the key points
    I = [1, -1];
    C = [1, 1];
    A = [-1, 1];

    JJ=JJ-50;
    PathDistance=sqrt((I(1)-C(1))^2+(I(2)-C(2))^2)+sqrt((C(1)-A(1))^2+(C(2)-A(2))^2)+sqrt((A(1)-I(1))^2+(A(2)-I(2))^2);
    TotalTime=(JJ)*ts;
    Speed=PathDistance/TotalTime;
    num_steps = JJ;

    change_point_1=round(sqrt((I(1)-C(1))^2+(I(2)-C(2))^2)/Speed/ts);
    change_point_2=round(sqrt((C(1)-A(1))^2+(C(2)-A(2))^2)/Speed/ts)+change_point_1;
    change_point_3=round(sqrt((A(1)-I(1))^2+(A(2)-I(2))^2)/Speed/ts)+change_point_2;

    TrueTrajectory = zeros(num_steps,2);
    
    % for i = 1: 50
    %     TrueTrajectory(i,1) = I(1);
    %     TrueTrajectory(i,2) = I(2);  
    % end
    % TrueTrajectory(1,1) = I(1);
    % TrueTrajectory(1,2) = I(2);  

    for i = 1: num_steps+1
        if i <= change_point_1
            TrueTrajectory(i,1) = I(1)  ;
            TrueTrajectory(i,2) = I(2) +(i-1)* Speed*ts; 
            

        elseif i <= change_point_2
            TrueTrajectory(i,1) = C(1) - (i-change_point_1)* Speed*ts;
            TrueTrajectory(i,2) = C(2) ;
        else 
            TrueTrajectory(i,1) = A(1)+ (i-change_point_2)* Speed*ts/sqrt(2);
            TrueTrajectory(i,2) = A(2) - (i-change_point_2)* Speed*ts/sqrt(2);
        end
    end

    padding = repmat([1, -1], 50, 1);
    TrueTrajectory = [padding; TrueTrajectory];

    % figure;
    % plot(I(1),I(2),'o');
    % hold on;
    % plot(C(1),C(2),'o');
    % plot(A(1),A(2),'o');

    % plot(TrueTrajectory(:,1),TrueTrajectory(:,2),'r.');
    % hold off

    % size(TrueTrajectory)
    % Human width
    width = 0.37;
    radius = width / 2;
    in_movement_area = zeros(length(Estimated_trajectory),1);
    distance_list = zeros(length(Estimated_trajectory),1);
    for i = 1: length(Estimated_trajectory)
        distance=sqrt((TrueTrajectory(i,1)-Estimated_trajectory(i,1))^2+(TrueTrajectory(i,2)-Estimated_trajectory(i,2))^2);
        distance_list(i) = distance;
        if distance <= radius
            in_movement_area(i) = 1;
        else
            in_movement_area(i) = 0;
        end
    end
    MeanDistance=mean(distance_list);
    % figure;
    % plot(distance_list);
    % plot(I(1),I(2),'o');
    % hold on;
    % plot(C(1),C(2),'o');
    % plot(A(1),A(2),'o');
    % plot(TrueTrajectory(:,1),TrueTrajectory(:,2),'r.');
    % plot(Estimated_trajectory(:,1),Estimated_trajectory(:,2),'b.');
    % hold off;

    % figure;
    % plot(Estimated_trajectory(:,1));
    % hold on;
    % plot(TrueTrajectory(:,1));
    % hold off;

    % figure;
    % plot(Estimated_trajectory(:,2));
    % hold on;
    % plot(TrueTrajectory(:,2));
    % hold off;

    delay = finddelay(Estimated_trajectory(:,1),TrueTrajectory(:,1))+finddelay(Estimated_trajectory(:,2),TrueTrajectory(:,2));
    delay = delay/2
    % % Define the inner vertices
    % I_inner = [0.815, -0.5509];
    % C_inner = [0.815, 0.815];
    % A_inner = [-0.5534, 0.815];
    % inner_vertices = [I_inner; C_inner; A_inner; I_inner];

   
    num_points_inside = sum(in_movement_area);
    false_alarm_rate =1- num_points_inside / length(Estimated_trajectory);
    disp(['False Alarm Rate: ', num2str(false_alarm_rate)]);
    disp(['Number of points inside the movement area: ', num2str(num_points_inside)]);
    % hold off;

end