# Estimation-of-Atmospheric-Disturbances
 ALS tunedcolored noise Kalman filter for estimation the atmospheric disturbances

 This repository contains the Matlab code to use dto generate the analysis in "Estimation of Atmospheric Disturbances
Encountered by a Large Commercial Aircraft" using the University of Washington's Kirsten Wind Tunnel Gust Load Alleviation Test bed.

The wind tunnel data and state space model used in this repository are from:   
`Hinson, K. A., & Morgansen, K. A. (2024). A Flexible Wing Model Uncertainty Evaluation Based on an Autocovariance Least Squares Tuned Optimal Estimate. 2024 American Control Conference (ACC), 2146–2151. https://doi.org/10.23919/ACC60939.2024.`

Wind Tunnel Test-bed state-space model: [LARGE.mat](https://github.com/kimberhinson/Estimation-of-Atmospheric-Disturbances/blob/main/LARGE.mat)

### Simulation

1. Generate simulated data sets, [simulated_setup.m](https://github.com/kimberhinson/Estimation-of-Atmospheric-Disturbances/blob/main/simulated_setup.m)
2. Run ALS for simulated data:
   - [run_als_sim.m](https://github.com/kimberhinson/Estimation-of-Atmospheric-Disturbances/blob/main/run_als_sim.m)
   - [setup_ALS_sim.m](https://github.com/kimberhinson/Estimation-of-Atmospheric-Disturbances/blob/main/setup_ALS_sim.m)
   - [als.m](https://github.com/kimberhinson/Estimation-of-Atmospheric-Disturbances/blob/main/als.m)
   - Plot identified nose covariances, [plot_lags_sim.m](https://github.com/kimberhinson/Estimation-of-Atmospheric-Disturbances/blob/main/plot_lags_sim.m)
     ![Q_sim_legend](https://github.com/user-attachments/assets/bdd5edfd-12a0-480a-9e45-280c37929150)
     ![R_sim_legend](https://github.com/user-attachments/assets/3e1972f4-604d-48f2-b2b7-2c6f9c426ae9)
   - Assess ALS conditioning, [svp_A.m](https://github.com/kimberhinson/Estimation-of-Atmospheric-Disturbances/blob/main/svp_A.m)
     ![Asv_sim](https://github.com/user-attachments/assets/6eb3d2ec-d4cb-4865-939a-48dbe7fccd18)
3. Evaluate performance of ALS tuned Kalman filters:
   - [simulated_setup_post.m](https://github.com/kimberhinson/Estimation-of-Atmospheric-Disturbances/blob/main/simulated_setup_post.m)
   - K-L Divergence metric, [KF_optimality_LTI.m](https://github.com/kimberhinson/Estimation-of-Atmospheric-Disturbances/blob/main/KF_optimality_LTI.m)
   - K-L Divergence plot, [E_opt.m](https://github.com/kimberhinson/Estimation-of-Atmospheric-Disturbances/blob/main/E_opt_plot.m)
     ![Eopt_sim_compare](https://github.com/user-attachments/assets/539ece68-90bb-4458-963e-ff4937a76498)
   - [simulated_setup_post_colored_ampcompare](https://github.com/kimberhinson/Estimation-of-Atmospheric-Disturbances/blob/main/simulated_setup_post_colored_ampcompare.m)
     ![Gust_sim_compare](https://github.com/user-attachments/assets/a89b7a2d-5d3b-4e87-98c9-1c302857c7f1)

### Wind Tunnel
1. Data is here:
2. Run ALS for wind tunnel data:
   - [run_als_WT.m](https://github.com/kimberhinson/Estimation-of-Atmospheric-Disturbances/blob/main/run_als_WT.m)
   - [setup_ALS_WT.m](https://github.com/kimberhinson/Estimation-of-Atmospheric-Disturbances/blob/main/setup_ALS_WT.m)
   - [als.m](https://github.com/kimberhinson/Estimation-of-Atmospheric-Disturbances/blob/main/als.m)
   - Plot identified nose covariances, [plot_lags_WT.m](https://github.com/kimberhinson/Estimation-of-Atmospheric-Disturbances/blob/main/plot_lags_WT.m)
   ![Q_WT_legend](https://github.com/user-attachments/assets/f3527a96-6fab-4f54-9f35-75f4d9db2fd6)
   ![R_WT_legend](https://github.com/user-attachments/assets/54a2aa01-1cb1-4441-8df8-8799adc42aa2)
3. Evaluate performance of ALS tuned Kalman filters
   - [WindTunnel_post_solored_gust.m](https://github.com/kimberhinson/Estimation-of-Atmospheric-Disturbances/blob/main/WindTunnel_post_colored_gust.m)
   ![Gust_est_WT_norm](https://github.com/user-attachments/assets/a2c22e00-8ced-417c-ab3c-3ec78067d882)



