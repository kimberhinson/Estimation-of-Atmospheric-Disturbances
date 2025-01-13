close all; clear all

load('./Data/WindTunnel/10psf_3aoa_CG3_1p5Hz_r2');

%% Load linear model

load LARGE.mat

dt = 1/100;
G = [0 0 0 0; 1 0 0 0; 0 0 0 0; 0 1 0 0; 0 0 0 0; 0 0 1 0; 0 0 0 0; 0 0 0 1];
H = zeros(3,4);

OBS = rank(obsv(A,C));                                                                      % Check observability
sysc_noise = ss(A,[B G],C,[D H]);                                                           % Continuous state-space model 
sys = c2d(sysc_noise,dt);                                                                   % Convert to discrete time ss with 1/100 time step 
c_scale = sqrtm(diag([10, 200, 1]));                                                        % Output scaling from "Ill-Conditioned Autocovariance Least Square Noise Covariance Identification for Colored Noise Gust Estimation," 

sysc_KF = ss(A,[B G],c_scale*C,[c_scale*D H]);                                              % Define scaled continuous state space plant model with process noise
sysKF = c2d(sysc_KF,dt);                                                                    % Convert to discrete time ss with 1/100 time step 

%% Colored noise Kalman filter

Ag = sysc_noise.A;
Ag(:,9) = [sysc_noise.B(:,7)];                                                              % Appending gust input to A matrix
Ag(9,:) = [zeros(1,8) 1];                                                                   % Random bias state
Bg = [sysc_noise.B(:,[1:6 8:11]), zeros(8,1); [zeros(1,10),1]];
Cg = [sysc_noise.C, sysc_noise.D(:,7)];
Dg = [sysc_noise.D(:,[1:6 8:11]),zeros(3,1)];

syscg_KF = ss(Ag,Bg,c_scale*Cg,c_scale*Dg); 
sysgKF = c2d(syscg_KF,dt);                                                                   % Convert to discrete time ss with 1/100 time step 

%% Wind Tunnel Test Data
   
tt = data{4}.Values.B2_A.Time(900:3:end)-data{4}.Values.B2_A.Time(900);                                                
      
[b_aa,a_aa] = butter(25,50/(300/2));                                                               % Anti aliasing before down sample

y_aa = [filtfilt(b_aa,a_aa,data{4}.Values.WT_A1.Data)...
        filtfilt(b_aa,a_aa,data{4}.Values.WT_A2.Data)...
        filtfilt(b_aa,a_aa,data{4}.Values.WT_A3.Data)];
     
y_noise = c_scale*[y_aa(900:3:end,1)-mean(y_aa(900:end,1))...
                   y_aa(900:3:end,2)-mean(y_aa(900:end,2))...
                   y_aa(900:3:end,3)-mean(y_aa(900:end,3))]';

u_aa = filtfilt(b_aa,a_aa,data{6}.Values.Data(900:3:end-90));
gust_cmd = [zeros(30,1); u_aa]*pi/180;                                                 

%% Estimator with ALS Q and R for colored noise Kalman filter

load('.\Results\WindTunnel_processed\meanQR.mat');

Qgest = zeros(5,5);
for ii = 1:4
Qgest(ii,ii) = meanQb(ii);                                                                                        % arbitrary guess of inital noise covariance
end
Qgest(5,5) = 1;
Rgest = zeros(3,3);
for jj = 1:3
Rgest(jj,jj) = meanRb(jj);
end
Rgest = c_scale'*Rgest*c_scale;
ugKF = zeros(6,length(tt));

[kgest,Lg,Pg,Mg] = kalman(sysgKF,Qgest,Rgest);                                                                    % define kalman filter

%% Play Data through ALS colored noise Estimator

[yg_est,tg_est,xg_est] = lsim(kgest,[ugKF; y_noise]',tt);  

xg_hat = [];
zyg_hat = [];
inng = [];

xg_hat =  xg_est';

% measurements
zyg_hat = sysgKF.C*xg_hat+sysgKF.D(:,1:6)*ugKF;

% innovations
zy = y_noise;
inng = zy-zyg_hat;      

outputs = {'WT1','WT2','WT3'};
    
for bb = 1:3
figure(bb)
subplot(2,1,1)
plot(tt,zy(bb,:))
hold on; grid on;
plot(tt, zyg_hat(bb,:),'-.')
legend('Measurement','ALS Estimate')
ylabel(outputs{bb},'Interpreter','none')
subplot(2,1,2)
plot(tt,zy(bb,:)-zyg_hat(bb,:),'-.','Color',[0.8500 0.3250 0.0980])
hold on; grid on;
ylabel('inn')
end

for bb = 9
figure(4) 
hold on; grid on;
plot(tt,gust_cmd./(max(gust_cmd)),'Color',[0 0.4470 0.7410])
plot(tt,xg_est(:,bb)./(max(gust_cmd)),'-.','Color',[0.4940 0.1840 0.5560],'LineWidth',1.5)
end
ylabel('Normalized Gust')
xlabel('Time (sec)')
legend('Actual','Estimate')

% set(gcf,'Units','inches');
% set(gcf,'Position',[0 1 12 5 ])
% screenposition = get(gcf,'Position');
% set(findall(gcf,'-property','FontSize'),'FontSize',15)
% set(gcf,...
%     'PaperPosition',[0 0 screenposition(3:4)],...
%     'PaperSize',[screenposition(3:4)]);
% print -dpng -painters -r400 Gust_est_WT_norm