function Plot_Locations(Position,Trajectory,Trajectory_Flag,Title,Save_File_path,Save_Flag)


    % Radar positions
    % x1=4.38/2+1; y1=0; %m
    % x2=(4.38/2+1)/sqrt(2); y2=-(4.38/2+1)/sqrt(2); %m
    % x3=0; y3=-(4.38/2+1); %m
    % x4=-(4.38/2+1)/sqrt(2); y4=-(4.38/2+1)/sqrt(2); %m
    % x5=-(4.38/2+1); y5=0; %m
    x1=-2;y1=2;
    x2=2;y2=2;
    x3=-2;y3=-2;
    x4=2;y4=-2;

    radar_x=[x1,x2,x3,x4];
    radar_y=[y1,y2,y3,y4];

    label_point = {'A','B','C','D','E','F','G','H','I'};
    label_point_x = [-1,0,1,-1,0,1,-1,0,1];
    label_point_y = [1,1,1,0,0,0,-1,-1,-1];

    figure;
    % Set the figure width and height
    set(gcf, 'Units', 'pixels', 'Position', [0, 0, 1920/1.5, 1080/1.5]);

    hold on;
    PlotTrueRange();
    % plot(Position(:,1),Position(:,2),'.','Color','b','DisplayName','Position');
    
    if Trajectory_Flag
        plot(Trajectory(:,1),Trajectory(:,2),'-','color','r','DisplayName','Trajectory','LineWidth',1.5);
    end
    plot(radar_x,radar_y,'rs',"linewidth",3,'DisplayName','Radar');
    text(x1,y1-0.2,'uwb 101')
    text(x2,y2-0.2,'uwb 102')
    text(x3,y3-0.2,'uwb 104')
    text(x4,y4-0.2,'uwb 106')
    
    % text(x5,y5-0.2,'radar 5')

    % Plot the label points
    plot(label_point_x,label_point_y,'S',"linewidth",2,'color','black','MarkerSize',7,'DisplayName','Label Points');
    for i=1:length(label_point)
        text(label_point_x(i)-0.1,label_point_y(i)+0.1,label_point(i));
    end
    hold off
    xlabel('x (m)'); ylabel('y (m)');
    ylim([-2.5,2.5]);
    xlim([-3,3]);
    title(Title);
    legend show
    grid on
    axis equal

    if Save_Flag
        set(gcf, 'PaperPositionMode', 'auto');
        print('-dpng', fullfile(Save_File_path, Title), '-r300');
    end
