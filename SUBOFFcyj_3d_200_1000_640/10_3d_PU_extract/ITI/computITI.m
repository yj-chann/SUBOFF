timestr='0';
probeNum=15;
filename = sprintf('D:\\OpenFOAM\\v2312\\SUBOFF\\SUBOFFcyj_3d_200_1000_640\\09_3d_LES\\postProcessing\\probes_line\\%s\\U', timestr);
fid = fopen(filename, 'r');
formatSpec = '# Probe %*d (%f %*f %*f)';
headerData = textscan(fid, formatSpec, probeNum);
fclose(fid);
x_coords = headerData{1};

fid = fopen(filename, 'r');
formatSpec = ['%f', repmat(' (%f %f %f)', 1, probeNum)];
dataCell = textscan(fid, formatSpec, 'HeaderLines', probeNum+1);
fclose(fid);
dataMatrix = cell2mat(dataCell);

% Extract time and U_data

time = dataMatrix(:, 1);
Ux = dataMatrix(:, 2:3:end-2);
Uy = dataMatrix(:, 3:3:end-1);
Uz = dataMatrix(:, 4:3:end);
plot(time,Ux(:,14))
hold on
plot(time,Ux(:,15))
plot(time,mean(Ux(:,14))*ones(size(time)))

%% 





[rms_Ux,mean_Ux]  = std(Ux, 1, 1);
[rms_Uy,mean_Uy]  = std(Uy, 1, 1);
[rms_Uz,mean_Uz]  = std(Uz, 1, 1);




% filename = 'D:\\graduationProject\\SUBOFFnjy\\09_3d_post\\InletTurbIntensity\\Tecplot_InputFiles\\InletTurbIntensity.plt';
% data = readmatrix(filename, 'FileType', 'text', 'NumHeaderLines', 3);
% time=data(:,1);
% Ux0=data(:,2);
% [rms_Ux0,mean_Ux0]=std(Ux0,1,1);
% Uy0=data(:,3);
% [rms_Uy0,mean_Uy0]=std(Uy0,1,1);
% Uz0=data(:,4);
% [rms_Uz0,mean_Uz0]=std(Uz0,1,1);





% --- 绘制湍流强度空间分布 ---
figure('Name', '湍流强度分布', 'NumberTitle', 'off');
% plot([-0.23,x_coords]/0.8712, [(rms_Ux0+rms_Uy0+rms_Uz0)/3/1.649*100 (rms_Uz+rms_Ux+rms_Uy)/3/1.649*100], 's-', 'LineWidth', 1.5, ...
%     'MarkerSize', 8, 'MarkerFaceColor', 'r');
% plot(x_coords/0.8712, (rms_Uz+rms_Ux+rms_Uy)/3/1.649*100, 's-', 'LineWidth', 1.5, ...
%     'MarkerSize', 8, 'MarkerFaceColor', 'r');
plot(x_coords/0.8712, (rms_Ux)/1.649*100, 's-', 'LineWidth', 1.5, ...
    'MarkerSize', 8, 'MarkerFaceColor', 'r');hold on;
plot(x_coords/0.8712, (rms_Uy)/1.649*100, 's-', 'LineWidth', 1.5, ...
    'MarkerSize', 8, 'MarkerFaceColor', 'r');
plot(x_coords/0.8712, (rms_Uz)/1.649*100, 's-', 'LineWidth', 1.5, ...
    'MarkerSize', 8, 'MarkerFaceColor', 'r');
xlabel('$x/L$','Interpreter','latex');
ylabel('$(U)_{\mathrm{rms}} / U_{\mathrm{ref}} \, (\%)$','Interpreter','latex');
ax=gca;
ax.FontName='Times New Roman';
ax.FontSize=20;
title('各采样点湍流强度空间分布');
grid on;

xline(0, '--k', '驻点');
hold off;
