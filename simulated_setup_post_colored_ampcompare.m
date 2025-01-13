close all; clear all

for amp = 1:3

%% Load linear model

load LARGE.mat

dt = 1/100;
G = [0 0 0 0; 1 0 0 0; 0 0 0 0; 0 1 0 0; 0 0 0 0; 0 0 1 0; 0 0 0 0; 0 0 0 1];
H = zeros(3,4);

OBS = rank(obsv(A,C));                                                                      % Check observability
sysc_noise = ss(A,[B G],C,[D H]);                                                           % Continuous state-space model 
sys = c2d(sysc_noise,dt);                                                                   % Convert to discrete time ss with 1/100 time step 
c_scale = sqrtm(diag([10, 200, 1]));                                                        % Output scaling from "Ill-Conditioned Autocovariance Least Square Noise Covariance Identification for Colored Noise Gust Estimation," 

%% Colored noise Kalman filter

Ag = sysc_noise.A;
Ag(:,9) = [sysc_noise.B(:,7)];                                                              % Appending gust input to A matrix
Ag(9,:) = [zeros(1,8) 1];                                                                   % Random bias state
Bg = [sysc_noise.B(:,[1:6 8:11]), zeros(8,1); [zeros(1,10),1]];
Cg = [sysc_noise.C, sysc_noise.D(:,7)];
Dg = [sysc_noise.D(:,[1:6 8:11]),zeros(3,1)];

syscg_KF = ss(Ag,Bg,c_scale*Cg,c_scale*Dg); 
sysgKF = c2d(syscg_KF,dt);         


%% Simulate Plant

Q = 0.01*diag([1,2,6,3]);                                                                                  % process and measurement noise to simulate
R = c_scale'*0.00001*diag([9,6,4])*c_scale;
                                                                            
tt = 0:dt:10;
u = zeros(7,length(tt));                                                                                   % construct input
gust_in = amp*1/2*sin(3*2*pi*tt); 

u(7,1:100) = zeros(1,100);
u(7,101:900) = gust_in(1:800);
u(7,901:1000) = zeros(1,100);

rng(1,'twister');
w = mvnrnd([0;0;0;0],Q,length(tt));                                                                        % process noise
v = mvnrnd(zeros(3,1),R,length(tt));                                                                       % measurement noise

[y,t,x] = lsim(sys,[u' w],tt);                                                                             % simulate the plant
y_noise = c_scale*(y') + v';                                                                               % add measurement noise

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
zy = y_noise;
zyg_hat = sysgKF.C*xg_hat+sysgKF.D(:,1:6)*ugKF;

% innovations
inng = zy-zyg_hat; 

%% Plots

for bb = 9
subplot(3,1,amp);
hold on; grid on;
plot(tt, xg_est(:,bb),'-.','Color',[0.4940 0.1840 0.5560])
plot(tt, u(7,:),'Color',[0 0.4470 0.7410])
ylabel('Gust');
ylim([-8 8])
if amp == 3
    xlabel('Time (sec)')
end
end
end
   

set(gcf,'Units','inches');
set(gcf,'Position',[0 1 12 5 ])
screenposition = get(gcf,'Position');
set(findall(gcf,'-property','FontSize'),'FontSize',15)
set(gcf,...
    'PaperPosition',[0 0 screenposition(3:4)],...
    'PaperSize',[screenposition(3:4)]);
print -dpng -painters -r400 Gust_sim_compare