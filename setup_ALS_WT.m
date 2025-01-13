function x = setup_ALS_WT(seed)


for N = [10 20 30 40 50 60 70 80 90 100 110 120 130 140 150]
    clearvars -global -except N seed datafiles
    clear ALS
    load(['./Data/WindTunnel_processed/LW_sys_Noise_' num2str(seed)]);

model.A = AA;
model.B = BB;
model.C = CC;
model.G = GG;
model.D = DD;
model.H = HH;
model.xhat0 = x_hat(:,1);
model.c_scale = c_scale;

data.yk = double(zy);
data.datapts = size(data.yk,2);
data.uk = zeros(size(BB,data.datapts));
data.start = 0;
data.xhatk = x_hat;
data.N = N;

estimator.Q = Qest;
estimator.R = Rest;

options.dt=  1/100;

ALS = als(data,model,estimator, options);
ALS.sdp_mrQ_diag

R = ALS.Rest_cell{1};
Q = ALS.Qest_cell{1};
A = ALS.A_scr;
b = ALS.Eyy;

if ~isfolder('./Results/WindTunnel_processed');mkdir('./Results/WindTunnel_processed');end
save(['./Results/WindTunnel_processed/M_lags' num2str(N) '_' num2str(seed)],'R','Q','A','model','b')
end

end

