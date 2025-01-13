classdef als < handle
%% Modified ALS-SDP function  
%
% Function Inputs: 
% data_yk:              measurements 
% data_uk:              inputs
% data_datapts:         number of data points considered
% data_start:           data to be ignored in the beginning until initial condition is negligible - optional, default is 100
% N:                    window size
% model.A:
% model.B:              
% model.C: 
% model.G: 
% model.D
% model.H            
% estimator.Q:          Qw initial guess
% estimator.R:          Rv initial guess
% options.dt:           Disctrete model time step
% 
% Function outputs : 
% Qest_cell:            cell containing estimated Qw 
% Rest_cell:            cell containing estimated Rv
% Phi:                  vector containing least-squares portion of objective for each penalty on tr(Q); scaled so that phi = 1 when rho = 0
    
    properties
        Ain                                         % Abar --> A-ALC
        data_yk                                     % measurement data
        data_uk                                     % input data
        data_xhat_                                  % calculated state estimate data
        data_datapts                                % number of points considered
        data_start                                  % data to be ignored in the beginning until initial condition is negligible
        datatrun                                    % data_datapts-data_start
        Eyy
        estimator_Q                                 % Initial process noise matrix 
        estimator_R                                 % Initial measurement noise matrix
        estimator_L                                 % Initial Kalman gain
        ga                                          % number of columns of G
        gamma
        inntrun
        M1
        M2
        model_Aa                                    % A matrix of model
        model_Ba                                    % B matrix of model
        model_Ca                                    % C matrix of model
        model_Ga                                    % G matrix of model
		model_Da 									% D matrix of model
		model_Ha									% H matrix of model
        model_xhat0                                 % initial state
        model_c_scale
        N                                           % Number of lags
        na                                          % number of rows of G
        norm_qest                                   % 2-norm of qest matrix
        norm_rest                                   % 2-norm or rest matrix
        OO                                          % Observability matrix where A is Ain
        options_dt                                  % Discrete time model time step
        pa                                          % number of columns of B
        Phi
        Qest_cell                                   % Process noise estimate
        Rest_cell                                   % Measurement noise estimate
        A_scr
        Wm
       
    end
    
    methods
        function obj = als(data,model,estimator, options)                                       % Constructor method
            obj.gamma = 0;
            obj.data_yk            = data.yk;
            obj.data_uk            = data.uk;
            obj.data_datapts       = data.datapts;
            obj.data_start         = data.start;
            obj.N                  = data.N;                
            obj.model_Aa           = model.A;
            obj.model_Ba           = model.B;
            obj.model_Ca           = model.C;
            obj.model_Ga           = model.G;
            obj.model_xhat0        = model.xhat0;
            obj.model_c_scale      = model.c_scale;
            
            obj.estimator_Q        = estimator.Q;
            obj.estimator_R        = estimator.R; 
            obj.options_dt         = options.dt;

                [~,~,~,obj.estimator_L] = kalman(ss(obj.model_Aa,[obj.model_Ba obj.model_Ga],obj.model_Ca, [obj.model_Da obj.model_Ha], options.dt),obj.estimator_Q,obj.estimator_R,0);
                % rajamani thesis eq 3.2 & 3.3: xhat(k+1) = A*xhat(k)+B*u(k)+A*L(y(k)-C(xhat(k)), matlab 'kalman funciton 2nd output is A*L and 4th output is L
            
            [obj.na,obj.ga]        = size(obj.model_Ga);
            obj.pa                 = size(obj.model_Ca,1);
            
            if isfield(data,'xhatk')
                obj.data_xhat_     = data.xhatk;
            else
                xhat               = zeros(obj.na,obj.data_datapts);
                xhat_              = zeros(obj.na,obj.data_datapts);
                xhat_(1:n,1)       = model.xhat0;  
                for i = 1:obj.data_datapts
                    xhat(:,i) = xhat_(:,i) + obj.estimator_L*(obj.data_yk(:,i)-obj.model_Ca*xhat_(:,i));
                    xhat_(:,i+1) = obj.model_Aa*xhat(:,i) + obj.model_Ba*obj.data_uk(:,i);
                end
                obj.data_xhat_ = xhat_(:,1:end-1);
            end
            
        end
         
         function obj = sdp_mrQ_diag(obj)
            obj.inntrun = obj.data_yk(:,obj.data_start+1:end)-obj.model_Ca*obj.data_xhat_(:,obj.data_start+1:end);      % Initial innovations          
            obj.oneColumn_autocorrelation
         	obj.Wm = eye(length(obj.Eyy));                                                                              % Use identity-based weighting
			obj.LSconstant_matrix
            obj.als_diag
         end
        
        function obj = oneColumn_autocorrelation(obj)                                                                   % Calculation of Autocorrelations for one column ALS
            
             obj.datatrun = obj.data_datapts - obj.data_start;
             obj.Eyy = [];
             for i = 0:obj.N-1
                 temp = obj.inntrun(:,i+1:end)*obj.inntrun(:,1:end-i)';                                                 % Rajamani thesis equation 3.14
                 temp = temp./(obj.datatrun-i);
                 obj.Eyy = [obj.Eyy; temp];
             end

             obj.Eyy = obj.Eyy(:);
        end
        
        function obj = LSconstant_matrix(obj)                                                                           % Building the constant matrix for the LS problem
           
        obj.Ain = obj.model_Aa-obj.model_Aa*obj.estimator_L*obj.model_Ca;                                               % Abar, Rajamani thesis equation 3.4
        
        obj.OO = [];
        temp = eye(obj.na);
        for i = 1:obj.N
            obj.OO = [obj.OO ; obj.model_Ca*temp];                                                                      % Script O, Rajamani thesis equation 3.12
            temp = temp*obj.Ain;
        end

        obj.M1 = zeros(obj.na^2,obj.ga^2);
        i=1;
        for j = 1:obj.ga
            for k = 1:obj.ga
                II = zeros(obj.ga);
                II(k,j) = 1;
                t1 = dlyap(obj.Ain,obj.model_Ga*II*obj.model_Ga');                                                      % Qw part of Lyapunov equation, Rajamani eq 3.8
                obj.M1(:,i) = t1(:);
                i = i+1;
            end
        end
        obj.M2 = zeros(obj.na^2,obj.pa^2);
        i=1;
        for j = 1:obj.pa
            for k = 1:obj.pa
                II = zeros(obj.pa);
                II(k,j) = 1;
                t2 = dlyap(obj.Ain,obj.model_Aa*obj.estimator_L*II*obj.estimator_L'*obj.model_Aa');                     % Rv part of Lyapunov equation, Rajamani eq 3.8
                obj.M2(:,i) = t2(:);
                i = i+1;
            end
        end
        
        end
     
        
        function obj = als_diag(obj)                                                                                    % Diagonal ALS
        
             PSI = eye(obj.pa);
             for  i= 1:obj.N-1
               PSI = [PSI; -obj.model_Ca*obj.Ain^(i-1)*obj.model_Aa*obj.estimator_L];
             end
 
             OOtemp = kron(obj.model_Ca,obj.OO);
             PSItemp = kron(eye(obj.pa),PSI);
 
             LH1 = OOtemp*obj.M1;   
             LH2 = OOtemp*obj.M2 + PSItemp;           
 
             As_diag = [LH1(:,1:obj.ga+1:obj.ga^2) LH2(:,1:obj.pa+1:obj.pa^2)];
 
             % Testing the uniqueness of covariance estimates  
             Arank = rank(As_diag, 1e-4);
 
             [nr, nc] = size(As_diag);
 
             if nc > Arank
               fprintf('Warning: Covariance estimates are not unique!\n')
               %pause
             end
 
            Xest_diag = quadprog(As_diag'*obj.Wm*As_diag, -As_diag'*obj.Wm*obj.Eyy, -eye(obj.pa+obj.ga), zeros(obj.ga+obj.pa,1));             
            Rss = obj.model_c_scale'*diag(Xest_diag(obj.ga+1:end))*obj.model_c_scale;
             
             if prod(Xest_diag)==0 
               fprintf('Warning: Covariance estimate(s) is (are) at constraints! You may have bad data! \n')
             end
             obj.Qest_cell{1} = diag(Xest_diag(1:obj.ga));
             obj.Rest_cell{1} = inv(obj.model_c_scale')*diag(Xest_diag(obj.ga+1:end))*inv(obj.model_c_scale);
 
             obj.A_scr = As_diag;       
        end  
    end
end

