close all

datafiles = {'10psf_3aoa_CG1p5_2p5Hz_r2.mat'...
    '10psf_3aoa_CG3_0p5Hz_r1.mat'...
    '10psf_3aoa_CG3_0p5Hz_r2.mat'...
    '10psf_3aoa_CG3_1p0Hz_r1.mat'...
    '10psf_3aoa_CG3_1p0Hz_r2.mat'...
    '10psf_3aoa_CG3_1p5Hz_r1.mat'...
    '10psf_3aoa_CG3_1p5Hz_r2.mat'...
    '10psf_3aoa_CG3_1p75Hz_r1.mat'...
    '10psf_3aoa_CG3_1p75Hz_r2.mat'...
    '10psf_3aoa_CG3_2p0Hz_r1.mat'...
    '10psf_3aoa_CG3_2p0Hz_r2.mat'...
    '10psf_3aoa_CG3_2p5Hz_r1.mat'...
    '10psf_3aoa_CG3_3p0Hz_r1.mat'...
    '10psf_3aoa_CG3_3p0Hz_r2.mat'...
    '10psf_3aoa_DG3_1p0Hz_r1.mat'...
    '10psf_3aoa_DG3_1p0Hz_r2.mat'...
    '10psf_3aoa_DG3_1p5Hz_r1.mat'...
    '10psf_3aoa_DG3_1p5Hz_r2.mat'...
    '10psf_3aoa_DG3_2p0Hz_r1.mat'...
    '10psf_3aoa_DG3_2p0Hz_r2.mat'...
    '10psf_3aoa_DG3_2p5Hz_r3.mat'...
    '10psf_3aoa_DG3_2p5Hz_r1.mat'...
    '10psf_3aoa_DG3_2p5Hz_r2.mat'...
    '10psf_3aoa_DG3_3p0Hz_r1.mat'...
    '10psf_3aoa_DG3_3p0Hz_r2.mat'...
    '10psf_3aoa_DG4_2p5Hz_r1.mat'...
    '10psf_3aoa_DG5_2p5Hz_r1.mat'...
    '10psf_3aoa_DG6_2p5Hz_r1.mat'...
    '10psf_3aoa_DG7_2p5Hz_r1.mat'};


for s = 1:length(datafiles)
%% Model
dt = 1/100;
dt_data = 1/300;

load LARGE.mat

G = [0 0 0 0; 1 0 0 0; 0 0 0 0; 0 1 0 0; 0 0 0 0; 0 0 1 0; 0 0 0 0; 0 0 0 1];
H = zeros(3,4);

OBS = rank(obsv(A,C));                                                                      % Check observability
sysc_noise = ss(A,[B G],C,[D H]);                                                           % Continuous state-space model 
sys = c2d(sysc_noise,dt);                                                                   % Convert to discrete time ss with 1/100 time step 
c_scale = sqrtm(diag([10, 200, 1]));                                                        % Output scaling from "Ill-Conditioned Autocovariance Least Square Noise Covariance Identification for Colored Noise Gust Estimation," 

sysc_KF = ss(A,[B G],c_scale*C,[c_scale*D H]);                                              % Define scaled continuous state space plant model with process noise
sysKF = c2d(sysc_KF,dt);                                                                    % Convert to discrete time ss with 1/100 time step 
%% Wind Tunnel Data
load(['./Data/WindTunnel/' datafiles{s}]);

tt = data{4}.Values.B2_A.Time(900:3:1/dt_data*10)- data{4}.Values.B2_A.Time(900);
u = zeros(7,length(tt));                                                  

[b_aa,a_aa] = butter(25,50/(300/2));                                                               % Anti aliasing before down sample

y_aa = [filtfilt(b_aa,a_aa,data{4}.Values.WT_A1.Data)...
        filtfilt(b_aa,a_aa,data{4}.Values.WT_A2.Data)...
        filtfilt(b_aa,a_aa,data{4}.Values.WT_A3.Data)];
     
y_noise = c_scale*[y_aa(900:3:1/dt_data*10,1)-mean(y_aa(900:1/dt_data*10,1))...
                   y_aa(900:3:1/dt_data*10,2)-mean(y_aa(900:1/dt_data*10,2))...
                   y_aa(900:3:1/dt_data*10,3)-mean(y_aa(900:1/dt_data*10,3))]';

%% Estimator with Inital Q and R

Qest = 0.01*eye(4);                                                         % arbitrary guess of inital noise covariance
Rest = 0.1*eye(3);
[kest,L,P,M] = kalman(sys,Qest,Rest);                                       % define kalman filter


%% Play Data through Suboptimal Estimator

[y_est,t_est,x_est] = lsim(kest,[u;y_noise],tt);  

x_hat = [];
y_hat = [];
y_tilde = [];
inn = [];

x_hat =  x_est';


% measurements
zy = y_noise;
zy_hat = sys.C*x_hat+sys.D(:,1:7)*u;

% innovations
inn = zy-zy_hat;       


if s == 1
figure 
ax1(1) = subplot(4,1,1);
plot(tt,zy(1,:))
hold on; grid on;
plot(tt, zy_hat(1,:),'-.')
legend('Measurement','Suboptimal Measurment Estimate')
ylabel('WT1')
ylim([-1 1])

ax1(2) = subplot(4,1,2);
plot(tt,zy(2,:))
hold on; grid on;
plot(tt, zy_hat(2,:),'-.')
ylabel('WT2')
ylim([-1 1])

ax1(3) = subplot(4,1,3);
plot(tt,zy(3,:))
hold on; grid on;
plot(tt, zy_hat(3,:),'-.')
ylabel('WT3')
ylim([-1 1])

ax1(5) = subplot(4,1,4);
hold on; grid on;
plot(tt,inn(1,:))
ylabel('Innovations')
ylim([-1 1])
end

AA = sysKF.A;
BB = sysKF.B(:,1:7);
GG = sysKF.B(:,8:11);
CC = sysKF.C;
DD = sysKF.D(:,1:7);
HH = sysKF.D(:,8:11);

if ~isfolder('./Data/WindTunnel_processed');mkdir('./Data/WindTunnel_processed');end
save(['./Data/WindTunnel_processed/LW_sys_Noise_' num2str(s) ],'AA','BB','CC','GG','DD','HH','x_hat','zy','c_scale','Qest','Rest');
end

