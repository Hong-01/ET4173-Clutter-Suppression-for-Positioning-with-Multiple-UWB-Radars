function HelperDetectionsMap(resp,rng_grid,dop_grid,rangeIndx,dopplerIndx,detections)
  % This function is only in support of CFARDetectionExample. It may be
  % removed in a future release.

  %   Copyright 2016 The MathWorks, Inc.

  figure
  detectionMap = zeros(size(resp));
  detectionMap(rangeIndx(1):rangeIndx(2),dopplerIndx(1):dopplerIndx(2)) = ...
    reshape(double(detections),rangeIndx(2)-rangeIndx(1)+1,dopplerIndx(2)-dopplerIndx(1)+1);
  h = imagesc(dop_grid,rng_grid,detectionMap);
  xlabel('Time bin'); ylabel('Range bin'); title('CFAR Detections');
  h.Parent.YDir = 'normal';
end