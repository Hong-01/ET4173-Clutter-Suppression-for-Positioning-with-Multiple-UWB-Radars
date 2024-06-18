% Clutter Suppression for Detection and Positioning with Multiple IR-UWB Radars
%-------------------------------------------------------------------------------------------------------------------
% This code is for the project of ET4173 Introduction to UWB Technology, Systems and Applications (2023/24 Q4)

clc;
clear ;
close all;
%--------------------------------------------------------------------------------------------------------------------
%% Import the data and Define Paramter
%--------------------------------------------------------------------------------------------------------------------
uwb101=load("E:\DATA\TUD\Master\TUD_Master_Y1\Q4\ET4173 Introduction to UWB Technology, Systems and Applications (202223 Q4) - 2252024 - 523 PM\Project\UWB Project Data\1\101_P1.mat");
uwb102= load("E:\DATA\TUD\Master\TUD_Master_Y1\Q4\ET4173 Introduction to UWB Technology, Systems and Applications (202223 Q4) - 2252024 - 523 PM\Project\UWB Project Data\1\102_P1.mat");
uwb104= load("E:\DATA\TUD\Master\TUD_Master_Y1\Q4\ET4173 Introduction to UWB Technology, Systems and Applications (202223 Q4) - 2252024 - 523 PM\Project\UWB Project Data\1\104_P1.mat");
uwb106= load("E:\DATA\TUD\Master\TUD_Master_Y1\Q4\ET4173 Introduction to UWB Technology, Systems and Applications (202223 Q4) - 2252024 - 523 PM\Project\UWB Project Data\1\106_P1.mat");

Nscansuwb101 = uwb101.scn.Nscn;Datauwb101=[uwb101.scn.scn];Datauwb101=reshape(Datauwb101, Nscansuwb101,[]);
Nscansuwb102 = uwb102.scn.Nscn;Datauwb102=[uwb102.scn.scn];Datauwb102=reshape(Datauwb102, Nscansuwb102,[]);
Nscansuwb104 = uwb104.scn.Nscn;Datauwb104=[uwb104.scn.scn];Datauwb104=reshape(Datauwb104, Nscansuwb104,[]);
Nscansuwb106 = uwb106.scn.Nscn;Datauwb106=[uwb106.scn.scn];Datauwb106=reshape(Datauwb106, Nscansuwb106,[]);

UWB_Data=cat(3,Datauwb101,Datauwb102,Datauwb104,Datauwb106);    % Concatenate the data from all the nodes

radar_nodes = ["101","102","104","106"];

[II,JJ,KK] = size(UWB_Data);    % Get the size of the data
fprintf('The slowtime bins are: \t%i \nthe range bins are: \t%i\nthe radar nodes are: \t%i\n', JJ,II,KK');

%Define the parameters
t = linspace(uwb101.scn(1).Tstrt, uwb101.scn(1).Tstp, uwb101.scn(1).Nscn)/1000; 
ts= t(2)-t(1); % fast time [ns], sample time
range_scope = 3e8*(t-t(1))/2e9; % range [m]
range_0 = min(range_scope); % min range in m
range_max = max(range_scope); % max range in m

%--------------------------------------------------------------------------------------------------------------------
%% Plot the range maps
%--------------------------------------------------------------------------------------------------------------------
%Plot the corresponding range maps with units
figure(1);
set(gcf, 'Units', 'pixels', 'Position', [0, 0, 1920, 1080]);
for k = 1:KK
    subplot(KK+1,1,k);
    imagesc([0,ts*JJ],[range_0,range_max],20*log10(abs(UWB_Data(:,:,k)))); axis xy
    colormap('jet'); axis xy; colorbar('east'); 
    xlabel("slowtime (sec)"); ylabel("range (m)");
    title("Radar node: "+radar_nodes(k));
end
set(gcf, 'PaperPositionMode', 'auto');
print('-dpng', fullfile('E:\DATA\TUD\Master\TUD_Master_Y1\Q4\ET4173 Introduction to UWB Technology, Systems and Applications (202223 Q4) - 2252024 - 523 PM\Project\Saved Images', 'Radar Range Time Plot (dB)'), '-r300');

figure(2);
set(gcf, 'Units', 'pixels', 'Position', [0, 0, 1920, 1080]);
for k = 1:KK
    subplot(KK+1,1,k);
    imagesc([0,ts*JJ],[range_0,range_max],abs(UWB_Data(:,:,k))); axis xy
    colormap('jet'); axis xy; colorbar('east'); 
    xlabel("slowtime (sec)"); ylabel("range (m)");
    title("Radar node: "+radar_nodes(k));
end
set(gcf, 'PaperPositionMode', 'auto');
print('-dpng', fullfile('E:\DATA\TUD\Master\TUD_Master_Y1\Q4\ET4173 Introduction to UWB Technology, Systems and Applications (202223 Q4) - 2252024 - 523 PM\Project\Saved Images', 'Radar Range Time Plot'), '-r300');
%% Cut the range maps
 % indexes based on previous plot (manual selection)
UWB_Data = UWB_Data(:,55:760,:);
[II,JJ,KK] = size(UWB_Data);    % Get the size of the data after cropping

%--------------------------------------------------------------------------------------------------------------------
%% Clutter Suppression
%--------------------------------------------------------------------------------------------------------------------

% 1. MTI filter
% w=[1, -0.6, -0.3, -0.1];
% for k = 1:KK
%     UWB_Data(:,:,k) = filter(w,1,UWB_Data(:,:,k),[],2);
% end

% figure;
% set(gcf, 'Units', 'pixels', 'Position', [0, 0, 1920, 1080]);
% for k = 1:KK
%     subplot(KK+1,1,k);
%     imagesc([0,ts*JJ],[range_0,range_max],abs(UWB_Data(:,:,k))); axis xy
%     colormap('jet'); axis xy; colorbar('east'); 
%     xlabel("slowtime (sec)"); ylabel("range (m)");
%     title("Radar node: "+radar_nodes(k));
% end
% set(gcf, 'PaperPositionMode', 'auto');
% print('-dpng', fullfile('E:\DATA\TUD\Master\TUD_Master_Y1\Q4\ET4173 Introduction to UWB Technology, Systems and Applications (202223 Q4) - 2252024 - 523 PM\Project\Saved Images', 'MTI Filtered Range Time Plot'), '-r300');


% 2. Adaptive Clutter Suppression
for k = 1:KK
    UWB_Data(:,:,k) = AdaptiveCS(UWB_Data(:,:,k),0.6,0.2);
end

figure(3);
set(gcf, 'Units', 'pixels', 'Position', [0, 0, 1920, 1080]);
for k = 1:KK
    subplot(KK+1,1,k);
    imagesc([0,ts*JJ],[range_0,range_max],(abs(UWB_Data(:,:,k)))); axis xy
    colormap('jet'); axis xy; colorbar('east'); 
    xlabel("slowtime (sec)"); ylabel("range (m)");
    title("Radar node: "+radar_nodes(k));
end
set(gcf, 'PaperPositionMode', 'auto');
print('-dpng', fullfile('E:\DATA\TUD\Master\TUD_Master_Y1\Q4\ET4173 Introduction to UWB Technology, Systems and Applications (202223 Q4) - 2252024 - 523 PM\Project\Saved Images', 'Adaptive Clutter Suppression Filtered Range Time Plot'), '-r300');


%--------------------------------------------------------------------------------------------------------------------
%% Normalization
%--------------------------------------------------------------------------------------------------------------------
for k = 1:KK
    UWB_Data(:,:,k) = abs(hilbert(UWB_Data(:,:,k))); 
    UWB_Data(:,:,k) = UWB_Data(:,:,k)/max(UWB_Data(:,:,k),[],'all'); 
end
%Plot the Normailized range maps
figure(4);
set(gcf, 'Units', 'pixels', 'Position', [0, 0, 1920, 1080]);
for k = 1:KK
    subplot(KK+1,1,k);
    imagesc([0,ts*JJ],[range_0,range_max],UWB_Data(:,:,k)); axis xy
    colormap('jet'); axis xy; colorbar('east'); 
    xlabel("slowtime (sec)"); ylabel("range (m)");
    title("Radar node: "+radar_nodes(k));
end
set(gcf, 'PaperPositionMode', 'auto');
print('-dpng', fullfile('E:\DATA\TUD\Master\TUD_Master_Y1\Q4\ET4173 Introduction to UWB Technology, Systems and Applications (202223 Q4) - 2252024 - 523 PM\Project\Saved Images', 'Normalized Range Time Plot'), '-r300');


%--------------------------------------------------------------------------------------------------------------------
%% Detection (NP)
%--------------------------------------------------------------------------------------------------------------------
%% 1. Centralized Detection 
estimated_range_np = zeros(KK,JJ);
for k = 1:KK    %No PF
    [estimated_range_np(k,:),rng] = cal_Ranges(UWB_Data(:,:,k), range_0, range_max);
    % % convert nan into 0
    estimated_range_np(k,isnan(estimated_range_np(k,:))) = 0;

    figure(5);
    set(gcf, 'Units', 'pixels', 'Position', [0, 0, 1920, 1080]);
    subplot(KK+1,1,k);
    imagesc([0,ts*JJ],rng,UWB_Data(:,:,k));title("Radar node: "+radar_nodes(k));axis xy;
    colormap('jet'); axis xy; colorbar('east'); 
    xlabel('slow time index [s]'); ylabel('range [m]')
    hold on
    plot([0:length(estimated_range_np(k,:))-1]*ts,estimated_range_np(k,:)-2,'r','LineWidth',1)
    % plot([0:length(estimated_range_np(k,:))-1]*ts,estimated_range_np(k,:)-2,'w','LineWidth',0.75)
    hold off
end
set(gcf, 'PaperPositionMode', 'auto');
print('-dpng', fullfile('E:\DATA\TUD\Master\TUD_Master_Y1\Q4\ET4173 Introduction to UWB Technology, Systems and Applications (202223 Q4) - 2252024 - 523 PM\Project\Saved Images', 'Estimated Range by NP Detector'), '-r300');
save('E:\DATA\TUD\Master\TUD_Master_Y1\Q4\ET4173 Introduction to UWB Technology, Systems and Applications (202223 Q4) - 2252024 - 523 PM\Project\Saved Data\estimated_range_np_CS.mat','estimated_range_np');

% 2. Decentralized Detection - Particle Filter
estimated_range_np_pf = zeros(KK,JJ);
% estimated_range_np_pf = load('E:\DATA\TUD\Master\TUD_Master_Y1\Q4\ET4173 Introduction to UWB Technology, Systems and Applications (202223 Q4) - 2252024 - 523 PM\Project\estimated_range_np_pf.mat').estimated_range_np;
for k = 1:KK    %PF
    % [estimated_range_np(k,:),rng] = cal_Ranges(UWB_Data(:,:,k), range_0, range_max);
    [~,rng] = cal_Ranges(UWB_Data(:,:,k), range_0, range_max);
    % % convert nan into 0
    % estimated_range_np(k,isnan(estimated_range_np(k,:))) = 0;

    %particle filter
    estimated_range_np_pf(k,:)=PF_1D(estimated_range_np(k,:),350,35);

    figure(6);
    set(gcf, 'Units', 'pixels', 'Position', [0, 0, 1920, 1080]);
    subplot(KK+1,1,k);
    imagesc([0,ts*JJ],rng,UWB_Data(:,:,k));title("Radar node: "+radar_nodes(k));axis xy;
    colormap('jet'); axis xy; colorbar('east'); 
    xlabel('slow time index [s]'); ylabel('range [m]')
    hold on
    plot([0:length(estimated_range_np_pf(k,:))-1]*ts,estimated_range_np_pf(k,:)-2,'r','LineWidth',1)
    % plot([0:length(estimated_range_np_pf(k,:))-1]*ts,estimated_range_np_pf(k,:),'w','LineWidth',0.75)
    hold off
end
set(gcf, 'PaperPositionMode', 'auto');
print('-dpng', fullfile('E:\DATA\TUD\Master\TUD_Master_Y1\Q4\ET4173 Introduction to UWB Technology, Systems and Applications (202223 Q4) - 2252024 - 523 PM\Project\Saved Images', 'Estimated Range by NP Detector with PF'), '-r300');
save('E:\DATA\TUD\Master\TUD_Master_Y1\Q4\ET4173 Introduction to UWB Technology, Systems and Applications (202223 Q4) - 2252024 - 523 PM\Project\Saved Data\estimated_range_np_pf_CS.mat','estimated_range_np_pf');



%--------------------------------------------------------------------------------------------------------------------
% %% Detection (CFAR)
%--------------------------------------------------------------------------------------------------------------------
% 1. Centralized Detection
estimated_range_cfar = zeros(KK,JJ);
% estimated_range_cfar=load('E:\DATA\TUD\Master\TUD_Master_Y1\Q4\ET4173 Introduction to UWB Technology, Systems and Applications (202223 Q4) - 2252024 - 523 PM\Project\estimated_range_cfar.mat').estimated_range_cfar;
% estimated_range_cfar_pf = load('E:\DATA\TUD\Master\TUD_Master_Y1\Q4\ET4173 Introduction to UWB Technology, Systems and Applications (202223 Q4) - 2252024 - 523 PM\Project\estimated_range_cfar_pf.mat').estimated_range_cfar;
for k = 1:KK    %No PF

    [estimated_range_cfar(k,:),rng] = cfar(UWB_Data(:,:,k), [10,15],[36,54],0.156,range_0,range_max,ts);

    % convert nan into 0
    estimated_range_cfar(k,isnan(estimated_range_cfar(k,:))) = 0;

    figure(7);
    set(gcf, 'Units', 'pixels', 'Position', [0, 0, 1920, 1080]);
    subplot(KK+1,1,k);
    imagesc([0,ts*JJ],rng,UWB_Data(:,:,k));title("Radar node: "+radar_nodes(k));axis xy;
    colormap('jet'); axis xy; colorbar('east'); 
    xlabel('slow time index [s]'); ylabel('range [m]')
    hold on
    plot([1:length(estimated_range_cfar(k,:))]*ts,estimated_range_cfar(k,:)-2,'r','LineWidth',1)
    % plot([1:length(estimated_range_cfar(k,:))]*ts,estimated_range_cfar(k,:)-2,'w','LineWidth',0.75)
    hold off
end
set(gcf, 'PaperPositionMode', 'auto');
print('-dpng', fullfile('E:\DATA\TUD\Master\TUD_Master_Y1\Q4\ET4173 Introduction to UWB Technology, Systems and Applications (202223 Q4) - 2252024 - 523 PM\Project\Saved Images', 'Estimated Range by CFAR Detector'), '-r300');
save('E:\DATA\TUD\Master\TUD_Master_Y1\Q4\ET4173 Introduction to UWB Technology, Systems and Applications (202223 Q4) - 2252024 - 523 PM\Project\Saved Data\estimated_range_cfar_CS.mat','estimated_range_cfar');


estimated_range_cfar_pf = zeros(KK,JJ);
for k = 1:KK    %with PF
    % [estimated_range_cfar(k,:),rng] = cfar(UWB_Data(:,:,k), 10,25,9e-2,range_0,range_max,ts);

    % % convert nan into 0
    % estimated_range_cfar(k,isnan(estimated_range_cfar(k,:))) = 0;

    %particle filter
    estimated_range_cfar_pf(k,:)=PF_1D(estimated_range_cfar(k,:),350,35);

    figure(8);
    set(gcf, 'Units', 'pixels', 'Position', [0, 0, 1920, 1080]);
    subplot(KK+1,1,k);
    imagesc([0,ts*JJ],rng,UWB_Data(:,:,k));title("Radar node: "+radar_nodes(k));axis xy;
    colormap('jet'); axis xy; colorbar('east'); 
    xlabel('slow time index [s]'); ylabel('range [m]')
    hold on
    plot([1:length(estimated_range_cfar_pf(k,:))]*ts,estimated_range_cfar_pf(k,:)-2,'r','LineWidth',1)
    % plot([1:length(estimated_range_cfar_pf(k,:))]*ts,estimated_range_cfar_pf(k,:)-2,'w','LineWidth',0.75)
    
    hold off
end
set(gcf, 'PaperPositionMode', 'auto');
print('-dpng', fullfile('E:\DATA\TUD\Master\TUD_Master_Y1\Q4\ET4173 Introduction to UWB Technology, Systems and Applications (202223 Q4) - 2252024 - 523 PM\Project\Saved Images', 'Estimated Range by CFAR Detector with PF'), '-r300');
save('E:\DATA\TUD\Master\TUD_Master_Y1\Q4\ET4173 Introduction to UWB Technology, Systems and Applications (202223 Q4) - 2252024 - 523 PM\Project\Saved Data\estimated_range_cfar_pf_CS.mat','estimated_range_cfar_pf');




%--------------------------------------------------------------------------------------------------------------------
%% localization by LS and PF
%--------------------------------------------------------------------------------------------------------------------

% 0. use np detector (centralized)
theta=Range2Position(estimated_range_np);

%Plot the estimated locations
Plot_Locations(theta,nan,0,'Estimated locations by NP detector with CS (Centralized)','E:\DATA\TUD\Master\TUD_Master_Y1\Q4\ET4173 Introduction to UWB Technology, Systems and Applications (202223 Q4) - 2252024 - 523 PM\Project\Saved Images',1);

%PF filter
% FilteredPostion_np=load('E:\DATA\TUD\Master\TUD_Master_Y1\Q4\ET4173 Introduction to UWB Technology, Systems and Applications (202223 Q4) - 2252024 - 523 PM\Project\FilteredPostion_np.mat').FilteredPostion_np;
FilteredPostion_np_c=PF_2D(theta,10000,50);
save('E:\DATA\TUD\Master\TUD_Master_Y1\Q4\ET4173 Introduction to UWB Technology, Systems and Applications (202223 Q4) - 2252024 - 523 PM\Project\Saved Data\FilteredPostion_np_CS.mat','FilteredPostion_np_c');

%Plot the estimated locations
Plot_Locations(theta,FilteredPostion_np_c,1,'Estimated locations by Particle filter and NP detector (Centralized) with CS','E:\DATA\TUD\Master\TUD_Master_Y1\Q4\ET4173 Introduction to UWB Technology, Systems and Applications (202223 Q4) - 2252024 - 523 PM\Project\Saved Images',1);




% 1. use np detector (Decentralized)
theta=Range2Position(estimated_range_np_pf);

%Plot the estimated locations
Plot_Locations(theta,nan,0,'Estimated locations by NP detector with CS (Decentralized)','E:\DATA\TUD\Master\TUD_Master_Y1\Q4\ET4173 Introduction to UWB Technology, Systems and Applications (202223 Q4) - 2252024 - 523 PM\Project\Saved Images',1);


%PF filter
% FilteredPostion_np=load('E:\DATA\TUD\Master\TUD_Master_Y1\Q4\ET4173 Introduction to UWB Technology, Systems and Applications (202223 Q4) - 2252024 - 523 PM\Project\FilteredPostion_np_pf.mat').FilteredPostion_np;
FilteredPostion_np_d=PF_2D(theta,10000,50);
save('E:\DATA\TUD\Master\TUD_Master_Y1\Q4\ET4173 Introduction to UWB Technology, Systems and Applications (202223 Q4) - 2252024 - 523 PM\Project\Saved Data\FilteredPostion_np_pf_CS.mat','FilteredPostion_np_d');
%Plot the estimated locations
Plot_Locations(theta,FilteredPostion_np_d,1,'Estimated locations by Particle filter and NP detector (Decentralized) with CS','E:\DATA\TUD\Master\TUD_Master_Y1\Q4\ET4173 Introduction to UWB Technology, Systems and Applications (202223 Q4) - 2252024 - 523 PM\Project\Saved Images',1);


% 2. use CFAR detector (centralized)
theta=Range2Position(estimated_range_cfar);

%Plot the estimated locations
Plot_Locations(theta,nan,0,'Estimated locations by CFAR detector with CS (Centralized)','E:\DATA\TUD\Master\TUD_Master_Y1\Q4\ET4173 Introduction to UWB Technology, Systems and Applications (202223 Q4) - 2252024 - 523 PM\Project\Saved Images',1);


%PF filter
% FilteredPostion_cfar=load('E:\DATA\TUD\Master\TUD_Master_Y1\Q4\ET4173 Introduction to UWB Technology, Systems and Applications (202223 Q4) - 2252024 - 523 PM\Project\FilteredPostion_cfar.mat').FilteredPostion_cfar;
FilteredPostion_cfar_c=PF_2D(theta,10000,50);
save('E:\DATA\TUD\Master\TUD_Master_Y1\Q4\ET4173 Introduction to UWB Technology, Systems and Applications (202223 Q4) - 2252024 - 523 PM\Project\Saved Data\FilteredPostion_cfar_CS.mat','FilteredPostion_cfar_c');
%Plot the estimated locations
Plot_Locations(theta,FilteredPostion_cfar_c,1,'Estimated locations by Particle filter and CFAR detector (Centralized) with CS','E:\DATA\TUD\Master\TUD_Master_Y1\Q4\ET4173 Introduction to UWB Technology, Systems and Applications (202223 Q4) - 2252024 - 523 PM\Project\Saved Images',1);


% 3. use CFAR detector (Decentralized)
theta=Range2Position(estimated_range_cfar_pf);

%Plot the estimated locations
Plot_Locations(theta,nan,0,'Estimated locations by CFAR detector with CS (Decentralized)','E:\DATA\TUD\Master\TUD_Master_Y1\Q4\ET4173 Introduction to UWB Technology, Systems and Applications (202223 Q4) - 2252024 - 523 PM\Project\Saved Images',1);

%PF filter
% FilteredPostion_cfar=load('E:\DATA\TUD\Master\TUD_Master_Y1\Q4\ET4173 Introduction to UWB Technology, Systems and Applications (202223 Q4) - 2252024 - 523 PM\Project\FilteredPostion_cfar_pf.mat').FilteredPostion_cfar;
FilteredPostion_cfar_d=PF_2D(theta,10000,50);
save('E:\DATA\TUD\Master\TUD_Master_Y1\Q4\ET4173 Introduction to UWB Technology, Systems and Applications (202223 Q4) - 2252024 - 523 PM\Project\Saved Data\FilteredPostion_cfar_pf_CS.mat','FilteredPostion_cfar_d');

%Plot the estimated locations
Plot_Locations(theta,FilteredPostion_cfar_d,1,'Estimated locations by Particle filter and CFAR detector (Decentralized) with CS','E:\DATA\TUD\Master\TUD_Master_Y1\Q4\ET4173 Introduction to UWB Technology, Systems and Applications (202223 Q4) - 2252024 - 523 PM\Project\Saved Images',1);

%--------------------------------------------------------------------------------------------------------------------
% False Alarm Rate
%--------------------------------------------------------------------------------------------------------------------

% 1. False Alarm Rate for NP detector
%Centralized
[false_alarm_rate_np_c,MeanDistance_np_c] = FalseAlarmRate(FilteredPostion_np_c,ts,JJ);

%Decentralized
[false_alarm_rate_np_d,MeanDistance_np_d]= FalseAlarmRate(FilteredPostion_np_d,ts,JJ);

% 2. False Alarm Rate for CFAR detector
%Centralized
[false_alarm_rate_cfar_c,MeanDistance_cfar_c] = FalseAlarmRate(FilteredPostion_cfar_c,ts,JJ);

%Decentralized
[false_alarm_rate_cfar_d,MeanDistance_cfar_d] = FalseAlarmRate(FilteredPostion_cfar_d,ts,JJ);

FA= [false_alarm_rate_np_c,false_alarm_rate_np_d,false_alarm_rate_cfar_c,false_alarm_rate_cfar_d];
MeanDistance= [MeanDistance_np_c,MeanDistance_np_d,MeanDistance_cfar_c,MeanDistance_cfar_d];

figure;
set(gcf, 'Units', 'pixels', 'Position', [0, 0, 1920/2, 1080/2]);
bar(FA);
xticklabels({'NP Centralized','NP Decentralized','CFAR Centralized','CFAR Decentralized'});
ylabel('False Alarm Rate');
% ylim([0,0.9]);
title('False Alarm Rate for NP and CFAR detectors with CS');
grid on;
for k = 1:length(FA)
    text(k, FA(k), num2str(FA(k), '%.4f'), ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom')
end
set(gcf, 'PaperPositionMode', 'auto');
print('-dpng', fullfile('E:\DATA\TUD\Master\TUD_Master_Y1\Q4\ET4173 Introduction to UWB Technology, Systems and Applications (202223 Q4) - 2252024 - 523 PM\Project\Saved Images', 'False Alarm Rate for NP and CFAR detectors with CS'), '-r300');

figure;
set(gcf, 'Units', 'pixels', 'Position', [0, 0, 1920/2, 1080/2]);
bar(MeanDistance);
xticklabels({'NP Centralized','NP Decentralized','CFAR Centralized','CFAR Decentralized'});
ylabel('Mean Distance Error [m]');
ylim([0,0.35]);
title('Mean Distance Error for NP and CFAR detectors with CS');
grid on;
for k = 1:length(FA)
    text(k, MeanDistance(k), num2str(MeanDistance(k), '%.4f'), ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom')
end
set(gcf, 'PaperPositionMode', 'auto');
print('-dpng', fullfile('E:\DATA\TUD\Master\TUD_Master_Y1\Q4\ET4173 Introduction to UWB Technology, Systems and Applications (202223 Q4) - 2252024 - 523 PM\Project\Saved Images', 'Mean Distance error for NP and CFAR detectors with CS'), '-r300');

%--------------------------------------------------------------------------------------------------------------------
% How many points are inside the region of human movement
%--------------------------------------------------------------------------------------------------------------------

PointNumber_np_c = NumberInRegion(FilteredPostion_np_c);

%Decentralized
PointNumber_np_d = NumberInRegion(FilteredPostion_np_d);

% 2. False Alarm Rate for CFAR detector
%Centralized
PointNumber_cfar_c = NumberInRegion(FilteredPostion_cfar_c);

%Decentralized
PointNumber_cfar_d = NumberInRegion(FilteredPostion_cfar_d);

Points_list= [PointNumber_np_c,PointNumber_np_d,PointNumber_cfar_c,PointNumber_cfar_d];

figure;
set(gcf, 'Units', 'pixels', 'Position', [0, 0, 1920/2, 1080/2]);
bar(Points_list);
xticklabels({'NP Centralized','NP Decentralized','CFAR Centralized','CFAR Decentralized'});
ylabel('Number of Points ');
title('Number of Points inside the region of human movement');
grid on;
for k = 1:length(FA)
    text(k, Points_list(k), num2str(Points_list(k), '%.f'), ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom')
end
set(gcf, 'PaperPositionMode', 'auto');
print('-dpng', fullfile('E:\DATA\TUD\Master\TUD_Master_Y1\Q4\ET4173 Introduction to UWB Technology, Systems and Applications (202223 Q4) - 2252024 - 523 PM\Project\Saved Images', 'Number of the points in the region for NP and CFAR detectors without CS'), '-r300');








% Analysis of the delay between the central and decentralized estimations
delay_np = (finddelay(FilteredPostion_np_c(:,1),FilteredPostion_np_d(:,1))+finddelay(FilteredPostion_np_c(:,2),FilteredPostion_np_d(:,2)))/2;
delay_cfar = (finddelay(FilteredPostion_cfar_c(:,1),FilteredPostion_cfar_d(:,1))+finddelay(FilteredPostion_cfar_c(:,2),FilteredPostion_cfar_d(:,2)))/2;

fprintf('The delay between the central and decentralized estimations for NP detector is: %f\n',delay_np);
fprintf('The delay between the central and decentralized estimations for CFAR detector is: %f\n',delay_cfar);


%--------------------------------------------------------------------------------------------------------------------
% END
%--------------------------------------------------------------------------------------------------------------------