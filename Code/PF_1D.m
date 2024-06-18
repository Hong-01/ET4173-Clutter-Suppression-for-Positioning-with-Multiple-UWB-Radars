function FilteredPostion=PF_1D(InputPosition,numParticles,windowSize)
        x_input = InputPosition;
        pf = stateEstimatorPF;
        pf.StateEstimationMethod = 'mean';
        pf.ResamplingMethod = 'systematic';
        initialize(pf,numParticles,x_input(1),eye(1));
        % 设置观测更新窗口大小
        observationBuffer = zeros(windowSize, 1);
        % 预测和校正的结果
        PosPred = zeros(length(x_input), 1);
        PosCorrected = zeros(length(x_input), 1);
        for i = 1:length(x_input)
            % 更新观测缓冲区
            observationBuffer(mod(i-1, windowSize) + 1) = x_input(i);
            % 预测步骤
            [PosPred(i), PosCov] = predict(pf);       
            % 校正步骤
            for k = 1:windowSize
                currentObservation = observationBuffer(k);
                if all(currentObservation == 0)
                    continue; % 跳过未填充的缓冲区条目
                end
                [PosCorrected(i), PosCov] = correct(pf, currentObservation);
            end
        end
        FilteredPostion=PosCorrected;
end