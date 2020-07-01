function err = phopq_1st_error(x)%,yfp_datasets)
%initialize error
e_dy = zeros(1,12);
X0 = zeros(1,17); X0del = X0;
%% Steady state(?) dose response simulation
% cultures grown overnight at 1mM Mg, diluted, resuspended in desired
% Mg and measured after 3.5 hours
% % %
k=x(26:31); x(31:36) = k;
% % %
x([27 28 29]) = [0 0 1]; x(30) = 0; 
% define all mutant parameters:
y = x; y(27) = 1; z = x; z([27 7]) = [1 0]; z2 = x; z2(7) = 0;% w = x; w(28) = 1; z3 = x; z3(14:19) = x(14:19)+k; z4 = z3; z4(27) = 1;
% clearvars f g fvc k w
% x: WT; % y: -mgrB; % z: -mgrB,-autoreg; % z2: -autoreg; % w: constitutive; % z3:-phosphatase; %z4 = -phosphatase;-mgrB
tss = [0 15]*3600;
x(29) = 1;
[~, X1] =ode15s(@phopq_1st, [0 12*3600], X0,{},x); % steady state at 1mM
x(29) = 30;
[~, Y30] =ode15s(@phopq_1st, tss, X1(end,:),{},x);
x(29) = 10;
[~, Y10]=ode15s(@phopq_1st, tss, X1(end,:),{},x);
x(29) = 3;
[~, Y3]=ode15s(@phopq_1st, tss, X1(end,:),{},x);
x(29) = 1;
[~, Y1]=ode15s(@phopq_1st, tss, X1(end,:),{},x);
x(29) = 0.3;
[~, Y03]=ode15s(@phopq_1st, tss, X1(end,:),{},x);
x(29) = 0.1;
[~, Y01]=ode15s(@phopq_1st, tss, X1(end,:),{},x);
x(29) = 0.03;
[~, Y003]=ode15s(@phopq_1st, tss, Y30(end,:),{},x);
% very low Mg steady state:
x(29) = 10^-4.5;
[~, Ymax]=ode15s(@phopq_1st, tss, Y30(end,:),{},x);
yfp_ss = [Y003(end,15)/Y003(end,17) Y01(end,15)/Y01(end,17) Y03(end,15)/Y03(end,17) Y1(end,15)/Y1(end,17) Y3(end,15)/Y3(end,17) Y10(end,15)/Y10(end,17) Y30(end,15)/Y30(end,17)]/(Y30(end,15)/Y30(end,17));
phopq_dr(:,1) = [0.03 0.1 0.3 1 3 10 30]; % Mg (mM)
phopq_dr(:,2) = [0.58 0.58 0.56 0.5 0.39 0.19 0.1]; % 
e_dy(11) = sum((yfp_ss' - phopq_dr(:,2)/phopq_dr(end,2)).^2); %/(norm(1 - phopq_dr(:,2)/phopq_dr(end,2),2)^2))*25;%/norm(phopq_dr(:,2)/phopq_dr(end,2),2);
e_dy(11) = e_dy(11)+((Ymax(end,15)/Ymax(end,17))/(Y30(end,15)/Y30(end,17))<6*(phopq_dr(1,2)/phopq_dr(end,2)))*25;
% err = sum(e_dy);
%% Dynamic simulation : WT
load('yfp data sets.mat')
yfp_3s = yfp_datasets(:,:,1);yfp_4s = yfp_datasets(:,:,2);
yfp_5s = yfp_datasets(:,:,3);yfp_6s = yfp_datasets(:,:,4);
yfp_7s = yfp_datasets(:,:,5);yfp_8s = yfp_datasets(:,:,6);
yfp_9s = yfp_datasets(:,:,7);yfp_10s = yfp_datasets(:,:,8);
yfp_11s = yfp_datasets(:,:,1);yfp_12s = yfp_datasets(:,:,2);
yfp_1s = yfp_datasets(:,:,11);yfp_2s = yfp_datasets(:,:,12);
tspan = [0 3]*3600; 
% 50 --> 0.01
x(29) = 2;
[~, X1] =ode15s(@phopq_1st, [0 8*3600], X0,{},x);
x(29) = 50;
[~, Y50] =ode15s(@phopq_1st, tspan, X1(end,:),{},x);
x(29) = 0.01;
[~, Y001]=ode15s(@phopq_1st, yfp_1s(:,1)*60, Y50(end,:),{},x); % for PmgrB
[~, Y001_pq]=ode15s(@phopq_1st, yfp_9s(:,1)*60, Y50(end,:),{},x); % for Pphopq
% 2 --> 0.01
x(29) = 2;
[~, Y2] =ode15s(@phopq_1st, [0 11]*3600 , X0,{},x);
x(29) = 0.01;
[~, Y2_001]=ode15s(@phopq_1st, yfp_2s(:,1)*60, Y2(end,:),{},x);
% 50 --> 2
x(29) = 2;
[~, Y50_2]=ode15s(@phopq_1st, yfp_5s(:,1)*60, Y50(end,:),{},x);
% 50 --> 10
x(29) = 10;
[~, Y10]=ode15s(@phopq_1st, yfp_6s(:,1)*60, Y50(end,:),{},x);
% 4 errors for 4 pmgrb-yfp output
e_dy(1) = norm((Y001(:,15)./Y001(:,17))/(Y50(end,15)/Y50(end,17)) - yfp_1s(:,2)/yfp_1s(1,2)); %/norm(yfp_1s(:,2)/yfp_1s(1,2));
e_dy(2) = norm((Y2_001(:,15)./Y2_001(:,17))/(Y50(end,15)/Y50(end,17)) - yfp_2s(:,2)/yfp_1s(1,2));%/norm(yfp_2s(:,2)/yfp_1s(1,2));
e_dy(5) = norm((Y50_2(:,15)./Y50_2(:,17))/(Y50(end,15)/Y50(end,17)) - yfp_5s(:,2)/yfp_1s(1,2));%/norm(yfp_5s(:,2)/yfp_1s(1,2));
e_dy(6) = norm((Y10(:,15)./Y10(:,17))/(Y50(end,15)/Y50(end,17)) - yfp_6s(:,2)/yfp_1s(1,2));%/norm(yfp_6s(:,2)/yfp_1s(1,2));
e_dy([5 6])=e_dy([5 6])*0.25; % lower weightage for 50-->2, 50-->10
% one for pphopq-yfp output
e_dy(9) = norm((Y001_pq(:,16)./Y001_pq(:,17))/(Y50(end,15)/Y50(end,17)) - yfp_9s(:,2)/yfp_1s(1,2));%/norm(yfp_9s(:,2)/yfp_1s(1,2));
%% Dynamic simulation: delta mgrB
y(29) = 2;
[~, X1] =ode15s(@phopq_1st, [0 8*3600], X0del,{},y);
y(29) = 50;
[~, X_delmgr] =ode15s(@phopq_1st, tspan, X1(end,:),{},y);
y(29) = 0.01;
[~, Y_delmgr]=ode15s(@phopq_1st, yfp_3s(:,1)*60, X_delmgr(end,:),{},y);
[~, Y_delmgr_pq]=ode15s(@phopq_1st, yfp_10s(:,1)*60, X_delmgr(end,:),{},y);
e_dy(3) = norm((Y_delmgr(:,15)./Y_delmgr(:,17))/(Y50(end,15)/Y50(end,17)) - yfp_3s(:,2)/yfp_1s(1,2));%/norm(yfp_3s(:,2)/yfp_1s(1,2));
e_dy(10) = norm((Y_delmgr_pq(:,16)./Y_delmgr_pq(:,17))/(Y50(end,15)/Y50(end,17)) - yfp_10s(:,2)/yfp_1s(1,2));%/norm(yfp_10s(:,2)/yfp_1s(1,2)); % error for del mgrB Phopq-yfp
%% -autoreg
z2(29) = 2;
[~, X1] =ode15s(@phopq_1st, [0 8*3600], X0,{},z2);
z2(29) = 50;
[~, X_delautoreg] =ode15s(@phopq_1st, tspan, X1(end,:),{},z2);
z2(29) = 0.01;
[~, Y_delautoreg]=ode15s(@phopq_1st, yfp_7s(:,1)*60, X_delautoreg(end,:),{},z2);
e_dy(7) = norm((Y_delautoreg(:,15)./Y_delautoreg(:,17))/(Y50(end,15)/Y50(end,17)) - yfp_7s(:,2)/yfp_1s(1,2));%/norm(yfp_7s(:,2)/yfp_1s(1,2));
%% -autoreg;delta mgrB
% mgrB, MgrB
z(29) = 2;
[~, X1] =ode15s(@phopq_1st, [0 8*3600], X0del,{},z);
z(29) = 50;
[~, X_delmgrdelautoreg] =ode15s(@phopq_1st,tspan, X1(end,:),{},z);
z(29) = 0.01;
[~, Y_delmgrdelautoreg]=ode15s(@phopq_1st, yfp_8s(:,1)*60, X_delmgrdelautoreg(end,:),{},z);
e_dy(8) = norm((Y_delmgrdelautoreg(:,15)./Y_delmgrdelautoreg(:,17))/(Y50(end,15)/Y50(end,17)) - yfp_8s(:,2)/yfp_1s(1,2));%/norm(yfp_8s(:,2)/yfp_1s(1,2));
%% Total error
err = sum(([e_dy]).^2);
end