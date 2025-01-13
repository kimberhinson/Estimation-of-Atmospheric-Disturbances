close all; clear all

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


%% Simulate Plant

Q = 0.01*diag([1,2,6,3]);                                                                                  % process and measurement noise to simulate
R = c_scale'*0.00001*diag([9,6,4])*c_scale;
                                                                            
tt = 0:dt:10;
u = zeros(7,length(tt));                                                                                   % construct input
gust_in = 1/2*sin(3*2*pi*tt); 

u(7,1:100) = zeros(1,100);
u(7,101:900) = gust_in(1:800);
u(7,901:1000) = zeros(1,100);

rng(1,'twister');
w = mvnrnd([0;0;0;0],Q,length(tt));                                                                        % process noise
v = mvnrnd(zeros(3,1),R,length(tt));                                                                       % measurement noise

[y,t,x] = lsim(sys,[u' w],tt);                                                                             % simulate the plant
y_noise = c_scale*(y') + v';                                                                               % add measurement noise

%% Estimator with ALS Q and R

load('.\Results\SIM\meanQR.mat');

Qest = zeros(4,4);
for ii = 1:4
Qest(ii,ii) = meanQb(ii);                                                                                  % arbitrary guess of inital noise covariance
end
Rest = zeros(3,3);
for jj = 1:3
Rest(jj,jj) = meanRb(jj);
end
Rest = c_scale'*Rest*c_scale;

[kest,L,P,M] = kalman(sysKF,Qest,Rest);                                                                    % define kalman filter


%% Play Data through ALS Estimator

[y_est,t_est,x_est] = lsim(kest,[u; y_noise]',tt);  

x_hat = [];
zy_hat = [];
inn = [];

x_hat =  x_est';

% measurements
zy = y_noise;
zy_hat = sysKF.C*x_hat+sysKF.D(:,1:7)*u;

% innovations
inn = zy-zy_hat;

%% Estimator with Perfect Q and R

[kperf,~,Pperf,~] = kalman(sysKF,Q,R);                                                                    % define kalman filter

%% Play Data through Perfect Estimator

[y_perf,t_perf,x_perf] = lsim(kperf,[u; y_noise]',tt);  

x_hat_p = [];
zy_hat_p = [];
inn_p = [];

x_hat_p =  x_perf';

% measurements
zy_p = y_noise;
zyo_p = y' + v';
zy_hat_p = sysKF.C*x_hat_p+sysKF.D(:,1:7)*u;

% innovations
inn_p = zy-zy_hat_p; 

%% Estimator with ALS Q and R for colored noise Kalman filter

load('.\Results\SIM\meanQR.mat');

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
inng = zy-zyg_hat; 

%% Plots

outputs = {'WT1','WT2','WT3'};

S = sysKF.C*(Pperf)*sysKF.C'+R;

for bb = 1
subplot(2,1,1);
plot(tt,zy(bb,:))
hold on; grid on;
plot(tt, zy_hat(bb,:),'--')
plot(tt, zy_hat_p(bb,:),'-.')
plot(tt, zyg_hat(bb,:),':')

legend('Measurement','ALS Estimate','Ideal Estimate','Colored Noise Estimate')
ylabel(outputs{bb},'Interpreter','none')
ylim([-20 20])
subplot(2,1,2);
plot(tt,zy(bb,:)-zy_hat(bb,:),'-.','Color',[0.8500 0.3250 0.0980])
hold on; grid on;
plot(tt,zy(bb,:)-zy_hat_p(bb,:),':','Color',[0.9290 0.6940 0.1250])
plot(tt,zy(bb,:)-zyg_hat(bb,:),':','Color',[0.4940 0.1840 0.5560])
plot(tt,-3*(sqrt(S(1,1)))*ones(length(x_hat),1),'--k');
plot(tt,3*(sqrt(S(1,1)))*ones(length(x_hat),1),'--k');
ylabel('innovations')
ylim([-3 3])
end

[average_autocovariance_ALS, norm_diff_ALS] = KF_optimality_LTI(inn,sysKF.A,sysKF.C,Qest,Rest,P);
[average_autocovariance_ALSp, norm_diff_ALSp] = KF_optimality_LTI(inn_p,sysKF.A,sysKF.C,Q,R,Pperf);
[average_autocovariance_ALSg, norm_diff_ALSg] = KF_optimality_LTI(inng,sysgKF.A,sysgKF.C,Qgest,Rgest,Pg);

E_opt_plot
