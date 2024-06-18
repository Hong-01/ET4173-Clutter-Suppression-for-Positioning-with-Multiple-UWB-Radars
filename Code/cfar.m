function [rangeArray,rng ]=cfar(radar_data,guardBand,trainingBand,False_alarm,range_0,range_max,ts)
    % guardBand =10;
    % trainingBand = 25;
    % False_alarm = 9e-2;
    disp('CFAR Detection Processing...')
    copy_edge=80;
    radar_data_edge = radar_data;
    radar_data_edge=[radar_data(end-copy_edge+1:end,:);radar_data_edge];
    radar_data_edge=[radar_data_edge;radar_data(1:copy_edge,:)];
    radar_data_edge=[radar_data_edge(:,end-copy_edge+1:end),radar_data_edge];
    radar_data_edge=[radar_data_edge,radar_data_edge(:,1:copy_edge)];
    radar_data=radar_data_edge;

    cfar2D = phased.CFARDetector2D('GuardBandSize', guardBand, 'TrainingBandSize',trainingBand, ...
        'ProbabilityFalseAlarm', False_alarm,"ThresholdFactor","Auto");
    JJ = size(radar_data, 2);
    II = size(radar_data, 1);
    timeGrid = transpose((0:JJ-1) * ts);
    rngGrid = linspace(range_0, range_max, II)';
    rng = linspace(range_0, range_max, II);
    % Determine valid index ranges considering guard and training cells

    totalGuardTraining = max(guardBand) + max(trainingBand);

    validRangeStart = totalGuardTraining + 1;
    validRangeEnd = II - totalGuardTraining;

    validTimeStart = totalGuardTraining + 1;
    validTimeEnd = JJ - totalGuardTraining;

    % Ensure the indices are within the valid range
    [~, rangeIndx] = min(abs(rngGrid - [0 4.8]));
    rangeIndx(1) = max(rangeIndx(1), validRangeStart);
    rangeIndx(2) = min(rangeIndx(2), validRangeEnd);

    [~, timeIndx] = min(abs(timeGrid - [0 (JJ-1) * ts]));
    timeIndx(1) = max(timeIndx(1), validTimeStart);
    timeIndx(2) = min(timeIndx(2), validTimeEnd);

    % Generate indices for CUT cells within the valid range
    [columnInds, rowInds] = meshgrid(timeIndx(1):timeIndx(2), rangeIndx(1):rangeIndx(2));

    


    % [columnInds, rowInds] = meshgrid(copy_edge+1:size(radar_data_edge,1)-copy_edge, copy_edge+1:size(radar_data_edge,2)-copy_edge);
    % columnInds = transpose(columnInds);
    % rowInds = transpose(rowInds);
    % size(columnInds)
    % size(rowInds)

    CUTIdx = [rowInds(:) columnInds(:)]';
    % CUTIdx = [rowInds(2:end-1) columnInds(2:end-1)]';


    % Run CFAR detector
    detections = cfar2D(radar_data, CUTIdx);
    % HelperDetectionsMap(radar_data,rngGrid,timeGrid,rangeIndx,timeIndx,detections)

    %save the detections


    % save('E:\DATA\TUD\Master\TUD_Master_Y1\Q4\ET4173 Introduction to UWB Technology, Systems and Applications (202223 Q4) - 2252024 - 523 PM\Project\detections.mat','detections');

    se = strel('square', 12);
    clean_detect = imopen(detections, se);
    % clean_detect = detections;
 
    % se=strel('square', 30);
    % clean_detect = imclose(clean_detect, se);
    % HelperDetectionsMap(radar_data,rngGrid,timeGrid,rangeIndx,timeIndx,clean_detect)
    detections_rs = reshape(clean_detect, rangeIndx(2) - rangeIndx(1) + 1, timeIndx(2) - timeIndx(1) + 1);
    % zero padding
    padLeft = rangeIndx(1) - 1;
    padRight = size(radar_data, 1) - rangeIndx(2);
    padTop = timeIndx(1) - 1;
    padBottom = size(radar_data, 2) - timeIndx(2);

    detections_padded = padarray(detections_rs, [padLeft, padTop], 0, 'pre');
    detections_padded = padarray(detections_padded, [padRight, padBottom], 0, 'post');
    estimated_detect = zeros(II, JJ);
    for i=1:JJ
        if sum(detections_padded(:,i))>0
            indicesOfOnes = find(detections_padded(:,i));
            middleIndex = indicesOfOnes(ceil(length(indicesOfOnes) /3));
            estimated_detect(middleIndex,i) = 1;
            % minIndex = min(indicesOfOnes);
            % estimated_detect(minIndex,i) = 1;
        end
    end

    % convert to range
    rangeArray = zeros(1, JJ);
    for i=1:JJ
        if sum(estimated_detect(:,i))>0
            indicesOfOnes = find(estimated_detect(:,i));
            rangeArray(i) = rngGrid(indicesOfOnes);
        end
    end

    %median filter
    rangeArray = medfilt1(rangeArray, 15);
    rangeArray(isnan(rangeArray)) = 0;  %replace NaN with 0
    %interpolation
    rangeArray = interp1(find(rangeArray), rangeArray(rangeArray~=0), 1:JJ, 'linear','extrap');
    rangeArray = medfilt1(rangeArray, 75);
    rangeArray = rangeArray(copy_edge+1:end-copy_edge)+2;

    rangeArray = filloutliers(rangeArray,'linear','median');

end