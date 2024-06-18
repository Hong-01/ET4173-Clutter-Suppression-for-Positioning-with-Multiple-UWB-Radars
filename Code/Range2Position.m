function [theta] = Range2Position(Estimated_Range)


    x1=-2;y1=2;
    x2=2;y2=2;
    x3=-2;y3=-2;
    x4=2;y4=-2;

    radar_x=[x1,x2,x3,x4];
    radar_y=[y1,y2,y3,y4];

    h=2*[x4-x1,y4-y1;
    x4-x2,y4-y2;
    x4-x3,y4-y3];

    theta=zeros(length(Estimated_Range(1,:)),2);
    for i=1:length(Estimated_Range(1,:))
        r1 = Estimated_Range(1,i);
        r2 = Estimated_Range(2,i);
        r3 = Estimated_Range(3,i);
        r4 = Estimated_Range(4,i);
        % r5 = Estimated_Range(5,i);
        x=[r1^2-r4^2-x1^2+x4^2-y1^2+y4^2;
        r2^2-r4^2-x2^2+x4^2-y2^2+y4^2;
        r3^2-r4^2-x3^2+x4^2-y3^2+y4^2];
        theta(i,:)=transpose((h'*h)\h'*x);
    end
    theta = theta(~any(isnan(theta),2),:);
end