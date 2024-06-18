function FilteredPostion=PF_2D(InputPosition,numParticles,windowSize)
        theta = InputPosition;
        pf = stateEstimatorPF;
        pf.StateEstimationMethod = 'mean';
        pf.ResamplingMethod = 'systematic'; % 可以尝试 'multinomial' 或 'stratified'
        % numParticles=100000;
        initialize(pf,numParticles,[theta(1,1) theta(1,2)],eye(2),'CircularVariables',[-pi,pi]);

        % for i = 1:length(theta)
        %     [PosPred(i,:),PosCov] =predict(pf);
        %     [PosCorrected(i,:),robotCov] = correct(pf,theta(i,:));
        % end

        % 设置观测更新窗口大小
        % windowSize = 100;
        observationBuffer = zeros(windowSize, 2);

        % 预测和校正的结果
        PosPred = zeros(length(theta), 2);
        PosCorrected = zeros(length(theta), 2);

        for i = 1:length(theta)
            % 更新观测缓冲区
            observationBuffer(mod(i-1, windowSize) + 1, :) = theta(i, :);
            
            % 预测步骤
            [PosPred(i,:), PosCov] = predict(pf);
            
            % 添加过程噪声
            % particles = pf.Particles;
            % for j = 1:numParticles
            %     particles(j, :) = particles(j, :) + mvnrnd([0 0], eye(2), 1); % 添加过程噪声，调整为适合你的系统的协方差矩阵
            % end
            % pf.Particles = particles;
            
            % 校正步骤
            for k = 1:windowSize
                currentObservation = observationBuffer(k, :);
                if all(currentObservation == 0)
                    continue; % 跳过未填充的缓冲区条目
                end
                [PosCorrected(i,:), PosCov] = correct(pf, currentObservation);
                
                % 添加测量噪声
                % weights = pf.ParticleWeights;
                % for j = 1:numParticles
                %     weights(j) = weights(j) * mvnpdf(currentObservation, PosPred(i,:), eye(2)); % 使用合适的测量噪声协方差矩阵
                % end
                % pf.ParticleWeights = weights / sum(weights); % 归一化权重
            end
        end
        FilteredPostion=PosCorrected;
end