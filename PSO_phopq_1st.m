% function PSO_phopq_1st
% global lb ub
%% 1. k2 = f(mg). k2b = g*k2;
% % %     %g kmdeg ktlnA ktlnE kbtpn1 kbtpn2 f1  K1    ktlnB  F  Km  k20 km2 k3 km3  k4 k5 km5 k6  K2  f2   n    ktlnY kb kd
% % lb = [-5 -3    -2.1  -2.5  0      -5.9	 1	 -0.65 0      0  -4  -3  -20 -3 -4   -2 -3 -4  -3  -1  1.0  0    -1.3  -2 -4];
% % ub = [ 0 -2.39 -1.1  -0.5  0      -4.5   1.4  0    0      3  0   1.1 -20  1 -1    1  1 -1   1   0  2    0.5  -1.3   1 -1];
% % 
% % nrep=50;
% % Solution=zeros(nrep+2,length(lb)+1); Solution(:,length(lb)+1) = 1e3;
% % Solution(nrep+1,:)=cat(2,lb, 1e5);
% % Solution(nrep+2,:)=cat(2,ub, 1e5);
% % options = optimoptions('particleswarm','SwarmSize',50,'Display','iter'); %,'InitialSwarmMatrix',in);%); %'PlotFcn',@pswplotbestf,
% % options.MaxIterations = 80; options.FunctionTolerance = 5e-3;
% % for j=1:nrep
% %         [X, fval] = particleswarm(@phopq_1st_error,length(lb),lb,ub,options);
% %         Solution(j,:)=cat(2,X, fval);
% %         save('11252019_phopq_1st_2','Solution')
% % end
% % Solution=sortrows(Solution,length(lb)+1);
% % 
% % % % % fmincon
% % % % load('0211_k1k2combo.mat'); x = Solution(2,:);
% % % % x([12 10 11 22]) = log10([0.01 1e3 0.5 2]);
% % % % options2 = optimoptions('fmincon','PlotFcn',@optimplotfval,'Display','iter','OptimalityTolerance',1e-2);
% % % % [s, fval] = fmincon(@phopq_1st_error,x,[],[],[],[],lb,ub,[],options2)

%% 2 k2,k5 = f(mg); k2b = k2*g;
%     %g kmdeg ktlnA ktlnE const kbtpn2 f1  K1    const F  Km  k20 km2 k3 km3  k4 k50 km5 k6  K2  f2  n   Km2 kb kd n2
% lb = [-5 -3    -2.1  -2.5  0     -5.9	1	-0.65 0     0  -4  -3  -20 -3 -4   -2 -2  -4  -3  -1  1.0 0   -4  -2 -4 0];
% ub = [ 0 -2.39 -1.1  -0.5  0     -4.5   1.4 0     0     3  0   1.1 -20  1 -1    1 1.3 -1   1   0  2   0.5  0   1 -1 0.5];
% 
% nrep=50;
% Solution=zeros(nrep+2,length(lb)+1); Solution(:,length(lb)+1) = 1e3;
% Solution(nrep+1,:)=cat(2,lb, 1e5);
% Solution(nrep+2,:)=cat(2,ub, 1e5);
% options = optimoptions('particleswarm','SwarmSize',50,'Display','iter'); %,'InitialSwarmMatrix',in);%); %'PlotFcn',@pswplotbestf,
% options.MaxIterations = 80; options.FunctionTolerance = 5e-3;
% for j=1:nrep
%         [X, fval] = particleswarm(@phopq_1st_error,length(lb),lb,ub,options);
%         Solution(j,:)=cat(2,X, fval);
%         save('12022019_phopq_1st_k2k5inp_2','Solution')
% end
% Solution=sortrows(Solution,length(lb)+1);

%% 3. k2 = f(Mg); k2b = g*k2; k3 km3  k4 k5 km5 k6
%     %g kmdeg ktlnA ktlnE kbtpn1 kbtpn2 f1  K1    ktlnB  F  Km  k20 km2 k3 km3  k4 k5 km5 k6  K2  f2   n    ktlnY kb kd  k3b km3b  k4b k5b km5b k6b
lb = [-5 -3    -2.1  -2.5  0      -5.9	 1	 -0.65 0      0  -4  -3  -20 -3 -4   -2 -3 -4  -3  -1  1.0  0    -1.3  -2 -4  -6  -4    -6  -6  -4   -6];
ub = [ 0 -2.39 -1.1  -0.5  0      -4.5   1.4  0    0      3  0   1.1 -20  1 -1    1  1 -1   1   0  2    0.5  -1.3   1 -1  -2   1    -1  -2   1   -1];

nrep=80;
Solution=zeros(nrep+2,length(lb)+1); Solution(:,length(lb)+1) = 1e3;
Solution(nrep+1,:)=cat(2,lb, 1e5);
Solution(nrep+2,:)=cat(2,ub, 1e5);
options = optimoptions('particleswarm','SwarmSize',50,'Display','iter'); %,'InitialSwarmMatrix',in);%); %'PlotFcn',@pswplotbestf,
options.MaxIterations = 100; options.FunctionTolerance = 5e-3;
for j=1:nrep
        [X, fval] = particleswarm(@phopq_1st_error,length(lb),lb,ub,options);
        Solution(j,:)=cat(2,X, fval);
        save('01172020_phopq_1st_all','Solution')
end
Solution=sortrows(Solution,length(lb)+1);

% % % fmincon
% % load('12172019_phopq_1st_all.mat'); x = Solution(1,1:end-1);
% % options2 = optimoptions('fmincon','PlotFcn',@optimplotfval,'Display','iter','OptimalityTolerance',1e-2);
% % [s, fval] = fmincon(@phopq_1st_error,x,[],[],[],[],lb,ub,[],options2)
