%%
clear all
close all
clc
%% Nice Colors for plots
%% .........................................................................
Color(:,1) = [85;170;170]/255;
Color(:,2) = [60;60;230]/255;
Color(:,3) = [170;0;170]/255;
Color(:,4) = [200;0;0]/255;
Color(:,5) = [50;80;125]/255;
Color(:,6) = [60;110;80]/255;
Color(:,7) = [220;125;0]/255;
Color(:,8) = [150;150;255]/255;
%% ........................................................................


% parameters
T = 4;
dt = 0.005;
Steps = T/dt;
gap = 100;

load('xTrue.txt')
load('yTrue.txt')
load('zTrue.txt')
t = (0:Steps-1)*dt; 


load('xData.txt')
load('yData.txt')
load('zData.txt')


t_data = [t(1:gap:end) t(end)+dt];

% figure
subplot(311),hold on, plot(t,xTrue,'Color',Color(:,1),'LineWidth',2)
hold on, plot(t_data,xData,'.','Color',Color(:,4),'MarkerSize',20)
subplot(312),hold on, plot(t,yTrue,'Color',Color(:,1),'LineWidth',2)
hold on, plot(t_data,yData,'.','Color',Color(:,4),'MarkerSize',20)
subplot(313),hold on, plot(t,zTrue,'Color',Color(:,1),'LineWidth',2)
hold on, plot(t_data,zData,'.','Color',Color(:,4),'MarkerSize',20)
set(gcf,'Color','w')

%% Prior
load('xPrior.txt')
load('yPrior.txt')
load('zPrior.txt')
subplot(311),hold on, plot(t,xPrior,'--','Color',Color(:,2),'LineWidth',2)
subplot(312),hold on, plot(t,yPrior,'--','Color',Color(:,2),'LineWidth',2)
subplot(313),hold on, plot(t,zPrior,'--','Color',Color(:,2),'LineWidth',2)

%% MAP 
load('xMAP.txt')
load('yMAP.txt')
load('zMAP.txt')
subplot(311),hold on, plot(t,xMAP,'Color',Color(:,3),'LineWidth',2)
subplot(312),hold on, plot(t,yMAP,'Color',Color(:,3),'LineWidth',2)
subplot(313),hold on, plot(t,zMAP,'Color',Color(:,3),'LineWidth',2)

figure(4)
for kk=0:99;
    xtmp = load(strcat('xSample',num2str(kk),'.txt'));
    subplot(311),hold on,plot(t,xtmp,'Color',Color(:,2),'LineWidth',1)
    
    ytmp = load(strcat('ySample',num2str(kk),'.txt'));
    subplot(312),hold on, plot(t,ytmp,'Color',Color(:,2),'LineWidth',1)
    
    ztmp = load(strcat('zSample',num2str(kk),'.txt'));
    subplot(313),hold on, plot(t,ztmp,'Color',Color(:,2),'LineWidth',1)
end

for kk=0:99
    xtmp = load(strcat('xpSample',num2str(kk),'.txt'));
    subplot(311),hold on, plot(t,xtmp,'Color',Color(:,3),'LineWidth',1)
    
    ytmp = load(strcat('ypSample',num2str(kk),'.txt'));
    subplot(312),hold on, plot(t,ytmp,'Color',Color(:,3),'LineWidth',1)
    
    ztmp = load(strcat('zpSample',num2str(kk),'.txt'));
    subplot(313),hold on, plot(t,ztmp,'Color',Color(:,3),'LineWidth',1)
end

subplot(311),hold on, plot(t,xTrue,'Color',Color(:,1),'LineWidth',2)
hold on, plot(t_data,xData,'.','Color',Color(:,4),'MarkerSize',20)
subplot(312),hold on, plot(t,yTrue,'Color',Color(:,1),'LineWidth',2)
hold on, plot(t_data,yData,'.','Color',Color(:,4),'MarkerSize',20)
subplot(313),hold on, plot(t,zTrue,'Color',Color(:,1),'LineWidth',2)
hold on, plot(t_data,zData,'.','Color',Color(:,4),'MarkerSize',20)
set(gcf,'Color','w')

load('weights.txt')
figure,plot(weights,'r.','MarkerSize',20)

%% compute mean
xMean = zeros(Steps,1);
yMean = zeros(Steps,1);
zMean = zeros(Steps,1);
for kk=0:999
    xtmp = load(strcat('xpSample',num2str(kk),'.txt'));
    xMean = xMean+weights(kk+1)*xtmp;
    ytmp = load(strcat('ypSample',num2str(kk),'.txt'));
    yMean = yMean+weights(kk+1)*ytmp;
    ztmp = load(strcat('zpSample',num2str(kk),'.txt'));
    zMean = zMean+weights(kk+1)*ztmp;
end
figure(1)
subplot(311),hold on, plot(t,xMean,'--','Color',Color(:,4),'LineWidth',2)
subplot(312),hold on, plot(t,yMean,'--','Color',Color(:,4),'LineWidth',2)
subplot(313),hold on, plot(t,zMean,'--','Color',Color(:,4),'LineWidth',2)


% H=[-188065 6.11209e+06 -1.19832e+06
%      6.12317e+06 -2.21872e+11 -2.15198e+10
%      -1.19972e+06 3.82246e+07 -7.56928e+06];
% H = 0.5*(H+H')
% [a b]=eig(H)
% b(1)
% b(2)
% b(3)
% binv =b\eye(3);binv(1,1)=0;
% a'*binv
