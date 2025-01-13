function [average_autocovariance, norm_diff, theoretical_covariance] = KF_optimality_LTI(e,A,C,Q,R,P)
%*********************************************************************
%This function checks the optimality of an LTI KF. Two conditions will be
%checked: 

%1) the autocovariance of the innovations is zero
%2) the covariance of the innovations is C*P*C'+R

%The inputs to this function are:
%e -> the innovations (n-by-length of time series)
%A,C -> The system matrices
%Q,R -> The noise covariance matrices
%P -> The state-error covariance matrix (steady state)
%*********************************************************************

%grab the sizes of the matrices
n = size(A,1);%number of states
m = size(C,1);%number of outputs
Nd = size(e,2);%number of samples


%define variables to be used:
N = 150;%number of lags
epsilon = 1e-3;%values smaller than this will be considered zero

%preallocate
sum_temp = zeros(m,m);
size_of_b_hatN = zeros(1,N);

%STEP 1: calculate the estimate of the autocovariance of the innovation
for jj = 1:N
    for ii = 1:Nd-jj
        sum_temp = e(:,ii)*e(:,ii+jj)'+ sum_temp;
    end 
        b_hatN = (1/(Nd-jj))*sum_temp;%average autocovariance for the jjth lag
        size_of_b_hatN(jj) = norm(b_hatN,1);
        sum_temp = zeros(m,m);%reset to zero for the next lag        
end
average_autocovariance = mean(size_of_b_hatN);

%STEP 2: calculate the average covariance of the innovations
sum_temp = zeros(m,m);
for ii = 1:Nd
    sum_temp = sum_temp + e(:,ii)*e(:,ii)';
end
estimated_covariance = (1/Nd)*sum_temp;

%theoretical covariance
theoretical_covariance = C*P*C' + R;

%calculate the norm of the difference between estimated and theoretical
%norm_diff = norm(estimated_covariance - theoretical_covariance,1);
%norm_diff = norm((estimated_covariance - theoretical_covariance)./theoretical_covariance,1);

u0 = [mean(e(1,:)); mean(e(2,:)); mean(e(3,:))];
u1 = zeros(3,1);

norm_diff = 1/2*(trace(inv(theoretical_covariance)*estimated_covariance)+(u1-u0)'*inv(theoretical_covariance)*(u1-u0)-3+log(det(theoretical_covariance)/det(estimated_covariance)));






