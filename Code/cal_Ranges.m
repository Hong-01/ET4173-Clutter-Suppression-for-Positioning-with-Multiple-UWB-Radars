function [estimated_range,rng] = cal_Ranges(pathuwb,range_0, range_max)

    [Nrng, Nscans] = size(pathuwb); 
    rng = linspace(range_0, range_max, Nrng);   % fast time [ns]
    
    % rng = 3e8*(t-t(1))/2e9; % range [m] 

    xtalk_rng = 1;       % crosstalk max range [m] 
    threshold = 0.035;    % threshold for detection 

    estimated_range = zeros(1, Nscans);
    for i1 = 1:Nscans   % Slow time index 
        i2 = find(rng <= xtalk_rng, 1, 'last'); % Avoiding picking antenna coupling  
        while pathuwb(i2, i1) < threshold && i2 < Nrng-1 
            i2 = i2 + 1; 
        end 
        estimated_range(i1) = rng(i2); 
    end 

    % estimated_range =2+ estimated_range;
    estimated_range = medfilt1(estimated_range, 50)+2 ; 
    estimated_range = filloutliers(estimated_range,'linear','median');

end

