% function PSO_phopq_0619_t
% global lb ub
% % % 1. Wide limits: k1k2 combo; kbtpn1 = 0.47*kbtpn2; k1b = f*k1, kb1 = kb2; kd1=f*kd2; k2b = g*k2
% % %1    %g     kmdeg ktlnA ktlnE kbtpn1 kbtpn2 f1   K1   ktlnB  k1  km1 k2  km2 k3 km3  k4 k5 km5 k6  K2  f2   kbtpn3 ktlnY
% % lb = [-5    -3    -2.1   -2.5  0      -5.9	 1	 -0.65 0      -4  -4  -3  -20 -3 -4   -2 -3 -4  -3  -1  1.0  -6.2   -1.3];
% % ub = [0     -2.39 -1.1   -0.5  0      -4.5   1.4 0     0       1  -1  1.1 -20  1 -1    1  1 -1   1   0  2    -6.2   -1.3];
% %             %kb kd   f   
% % lb(24:26) = [-2 -4   0];
% % ub(24:26) = [ 1 -1   0];

% % 2. same conditions as above, but reduced constraints
%      %kpdeg kmdeg ktlnA ktlnE kbtpn1 kbtpn2  f1   K1     ktlnB  k1   km1  k2   km2 k3   km3  k4   k5   km5  k6   K2   f2   kbtpn3 ktlnY
% lb = [-5     -3    -1.6  -2.5  -1     -5.6   1.2  -0.65  1      -2.5 -1.2 -0.5 -20 -0.5 -4.5 -0.5 -1.3 -4   -2.2 -1   1.6  -6.2   -1.3];
% ub = [0      -2.39 -1.0  -1.6  -1     -4.5   1.52  0     1      0    0.5  1.1  -20  0.7 -2    1   0.4  -2   0.2  -0.2 2    -6.2   -1.3];
%           %   kb   kd   f   
% lb(24:26) = [-1.6  -3.4 0]; 
% ub(24:26) = [ 0.5  -1   0 ];

% % % % 3. 02/28/2019 nonlinear mg
% %       %g     kmdeg ktlnA ktlnE n      kbtpn2 f1   K1     km     k1   km1  k2   km2 k3   km3  k4   k5   km5  k6   K2   f2   kbtpn3 ktlnY
% % lb = [-4     -3    -1.6  -2.5  0      -5.6	 1.2  -0.65  -3     -2.5 -1.2 -0.5 -20 -0.5 -4.5 -0.5 -1.3 -4   -2.2 -1   1.6  -6.2   -1.3];
% % ub = [0      -2.39 -1.0  -1.6  0.7    -4.5   1.52  0     2.39   0    0.5  1.1  -20  0.7 -2    1   0.4  -2   0.2	 -0.2 2    -6.2   -1.3];
% %           %   kb   kd   f   
% % lb(24:26) = [-1.6  -3.4 -2];
% % ub(24:26) = [ 0.5  -1   0 ];


% % % % 4. for k2k5
% %      %g     kmdeg ktlnA ktlnE kbtpn1 kbtpn2 f1   K1     ktlnB  k1   km1  k2   km2  k3  km3  k4   k5   km5  k6   K2   f2   kbtpn3 ktlnY
% % lb = [-4    -3    0     0     -1     -5.6	1.2  -0.65  1      -2.5 -1.2 -0.5 -20 -0.5 -3   -0.5 -1.3 -3   -2.2 -1   1.6  -6.2   -1.3];
% % ub = [0     -2.39 2     2     -1     -4.5   1.52  0     1      0    0.5  1.1  -20  1   -1.5  1    1   -1.5 0.2	-0.2 2    -6.2   -1.3];
% %           %   kb   kd   f    
% % lb(24:26) = [-1.6  -3   0 ];
% % ub(24:26) = [ 0.5  -1   2 ];

% % % 5. k1k2k3km3k4k5km5k6; Wide limits kbtpn1 = 0.47*kbtpn2; k1b = f*k1, kb1= kb2; kd1=f*kd2;
%1    k2b kmdeg ktlnA ktlnE <phoPQ> kbtpn2 f1   K1    <mgrB> k1  km1 k2  km2 k3 km3  k4 k5 km5 k6  K2  f2   [free] [free]
lb = [-6  -3    -2.1  -2.5  0       -5.9   1	-0.65 0      -4  -4  -3  -5  -3 -4   -2 -3 -4  -3  -1  1.0  -6.2   -1.3];
ub = [-3  -2.39 -1.1  -0.5  0       -4.5   1.4  0     0       1  -1  1.1 -2   1 -1    1  1 -1   1   0  2    -6.2   -1.3];
            %kb kd k1b k3b km3b k4b k5b km5b k6b    
lb(24:32) = [-2 -4 -6  -4  -5   -5  -4  -5   -5];
ub(24:32) = [ 1 -1 -3   1   1    1   1   1    1];


% load('yfp data sets.mat')
% loc = @(x) phopq_error_0619_t(x,yfp_datasets);
% load('0211_k1k2combo.mat'); x = Solution(1,:);
% in = [Solution(1:4,1:26) zeros(4,6)];
nrep=70;
Solution=zeros(nrep+2,length(lb)+1); Solution(:,length(lb)+1) = 1e3;
Solution(nrep+1,:)=cat(2,lb, 1e5);
Solution(nrep+2,:)=cat(2,ub, 1e5);
options = optimoptions('particleswarm','SwarmSize',50,'Display','iter'); %,'InitialSwarmMatrix',in);%); %'PlotFcn',@pswplotbestf,
options.MaxIterations = 100; options.FunctionTolerance = 5e-3;
for j=1:nrep
%         [X, fval] = particleswarm(@PQB_error,length(lb),lb,ub,options);
        [X, fval] = particleswarm(@phopq_error_0619_t,length(lb),lb,ub,options);
%         [X, fval] = particleswarm(loc,length(lb),lb,ub,options);
%         [X, fval] = particleswarm(@PhoQmutant_steadystates,length(lb),lb,ub,options);
        Solution(j,:)=cat(2,X, fval);
        save('12082019_all_0619t','Solution')
end
Solution=sortrows(Solution,length(lb)+1);
% % load('0211_k1k2combo.mat'); x = Solution(1,:); x(27:32) = 0;
% % lb = [-4     -3    -1.6  -2.5  -1     -5.6	 1.2  -0.65  1      -2.5 -1.2 -0.5 -20 -0.5 -4.5 -0.5 -1.3 -4   -2.2 -1   1.6  -6.2   -1.3];
% % ub = [0      -2.39 -1.0  -1.6  -1     -4.5   1.52  0     1      0    0.5  1.1  -20  0.7 -2    1   0.4  -2   0.2	 -0.2 2    -6.2   -1.3];
% %           %   kb   kd   f     km2b k3b km3b k4b
% % lb(24:26) = [-1.6  -3.4 -2]; % -1  -1   -1   -1];
% % ub(24:26) = [ 0.5  -1   0 ]; % 0   -1   -1   -1];
% % lb(27:32) = -5*ones(1,6); ub(27:32) = 2*ones(1,6);
% % load('yfp data sets.mat')
% % loc = @(x) phopq_error_0619_t(x,yfp_datasets);
% % options2 = optimoptions('fmincon','PlotFcn',@optimplotfval,'Display','iter','OptimalityTolerance',1e-3);
% % [s, fval] = fmincon(loc,x,[],[],[],[],lb,ub,[],options2)

% load('0331_k1k2_phoqt281r'); x = Solution(1,1:32);
% options2 = optimoptions('fmincon','Algorithm','sqp','PlotFcn',@optimplotfval,'Display','iter'); %,'OptimalityTolerance',1e-8);
% [s, fval] = fmincon(@PhoQmutant_steadystates,x,[],[],[],[],lb,ub,[],options2)