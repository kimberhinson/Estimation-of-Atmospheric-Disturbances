close all; clear all

for seed = 1:50

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

%% Simulate Plant

Q = 0.01*diag([1,2,6,3]);                                                                                  % process and measurement noise to simulate
R = c_scale'*0.00001*diag([9,6,4])*c_scale;
                                                                            
tt = 0:dt:50;
u = zeros(7,length(tt));                                                                                   % zero input

rng(seed,'twister');
w = mvnrnd([0;0;0;0],Q,length(tt));                                                                        % process noise
v = mvnrnd(zeros(3,1),R,length(tt));                                                                       % measurement noise

[y,t,x] = lsim(sys,[u' w],tt);                                                                             % simulate the plant
y_noise = c_scale*(y') + v';                                                                               % add measurement noise


%% Estimator with Inital Q and R

Qest = 0.01*eye(4);                                                                                        % arbitrary guess of inital noise covariance
Rest = 0.01*eye(3);
[kest,L,P,M] = kalman(sysKF,Qest,Rest);                                                                    % define kalman filter


%% Play Data through Suboptimal Estimator

[y_est,t_est,x_est] = lsim(kest,[u;y_noise],tt);  

x_hat = [];
y_hat = [];
y_tilde = [];
inn = [];

x_hat =  x_est';


% measurements
zy = y_noise;
zyo = y' + v';
zy_hat = sysKF.C*x_hat+sysKF.D(:,1:7)*u;

% innovations
inn = zy-zy_hat;       

if seed == 1
figure 
ax1(1) = subplot(4,1,1);
plot(tt,zy(1,:))
hold on; grid on;
plot(tt, zy_hat(1,:),'-.')
legend('Measurement','Measurment Estimate')
ylabel('WT1')


ax1(2) = subplot(4,1,2);
plot(tt,zy(2,:))
hold on; grid on;
plot(tt, zy_hat(2,:),'-.')
ylabel('WT2')


ax1(3) = subplot(4,1,3);
plot(tt,zy(3,:))
hold on; grid on;
plot(tt, zy_hat(3,:),'-.')
ylabel('WT3')


ax1(4) = subplot(4,1,4);
hold on; grid on;
plot(tt,inn(1,:))
plot(tt,inn(2,:))
plot(tt,inn(3,:))
ylabel('Innovations')
xlabel('Time (sec)')


figure 
axs(1) = subplot(4,2,1);
plot(tt,x(:,1))
hold on; grid on;
plot(tt, x_hat(1,:),'-.')
ylabel('\eta_1')


axs(2) = subplot(4,2,2);
plot(tt,x(:,2))
hold on; grid on;
plot(tt, x_hat(2,:),'-.')
ylabel('$\dot{\eta_1}$','Interpreter','latex')


axs(3) = subplot(4,2,3);
plot(tt,x(:,3))
hold on; grid on;
plot(tt, x_hat(3,:),'-.')
ylabel('\eta_2')


axs(4) = subplot(4,2,4);
plot(tt,x(:,4))
hold on; grid on;
plot(tt, x_hat(4,:),'-.')
ylabel('$\dot{\eta_2}$','Interpreter','latex')


axs(5) = subplot(4,2,5);
plot(tt,x(:,5))
hold on; grid on;
plot(tt, x_hat(5,:),'-.')
ylabel('\eta_3')


axs(6) = subplot(4,2,6);
plot(tt,x(:,6))
hold on; grid on;
plot(tt, x_hat(6,:),'-.')
ylabel('$\dot{\eta_3}$','Interpreter','latex')


axs(7) = subplot(4,2,7);
plot(tt,x(:,7))
hold on; grid on;
plot(tt, x_hat(7,:),'-.')
ylabel('\eta_4')
xlabel('Time (sec)')

axs(8) = subplot(4,2,8);
plot(tt,x(:,8))
hold on; grid on;
plot(tt, x_hat(8,:),'-.')
ylabel('$\dot{\eta_4}$','Interpreter','latex')
xlabel('Time (sec)')
legend('State','State Estimate')
end

AA = sysKF.A;
BB = sysKF.B(:,1:7);
GG = sysKF.B(:,8:11);
CC = sysKF.C;
DD = sysKF.D(:,1:7);
HH = sysKF.D(:,8:11);

if ~isfolder('./Data/SIM');mkdir('./Data/SIM');end
save(['./Data/SIM/LW_sys_Noise_' num2str(seed) ],'AA','BB','CC','GG','DD','HH','x_hat','zy','c_scale','Qest','Rest');

end