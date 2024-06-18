function y=AdaptiveCS(Data,a_max,a_min)
    % Data: II x JJ x KK matrix
    %Model: y(t)=r(t)-c(t), y(t): Clean signal, r(t): Received signal, c(t): Clutter signal
    [II,JJ,~]=size(Data);
    c=zeros(II,JJ);
    y=zeros(II,JJ);
    r=Data(:,:,1);
    for j=1:JJ
        r_j=r(:,j);
        if sum(r_j)~=0  % if the signal is not zero column, then do the following
            y(:,j)=r(:,j)-c(:,j);
            norm_env_r_j=abs(hilbert(r_j))/max(abs(hilbert(r_j)));
            norm_env_y_j=abs(hilbert(y(:,j)))/max(abs(hilbert(y(:,j))));
            for i=1:II
                d=min(norm_env_r_j(i),norm_env_y_j(i))/norm_env_r_j(i);
                a=a_min+(a_max-a_min)*d;
                c(i,j+1)=a*c(i,j)+(1-a)*r(i,j);
            end
        else
            y(:,j)=r(:,j);
            c(:,j+1)=c(:,j);
        end
    end