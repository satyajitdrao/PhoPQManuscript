function err = phopq_error_0619_t(x)%,yfp_datasets)
%initialize error
e_dy = zeros(1,12);
X0 = zeros(1,19); 
X0del = X0;
%% Steady state(?) dose response simulation
% cultures grown overnight at 1mM Mg, diluted, resuspended in desired
% Mg and measured after 3.5 hours
% % % % % IF Optimizing PhoQT281R parameters % % % % %
% k = x(27:32);
% % % % % ELSE comment out % % % % %
% x([27 28 29]) = [0 0 1]; f = 10^x(26); g = 10^x(1); fvc = ones(1,12); fvc([1 3 6]) = [g f f]; x(30:41) = fvc; x(42) = 0; % corrected detailed balance; k1b = f*k1; kd1 = kd2*f; k2b = g*k2;
% x([27 28 29]) = [0 0 1]; f = 10^x(26); g = 10^x(1); fvc = ones(1,12); fvc([1 10]) = [g f]; x(30:41) = fvc; x(42) = 0; % k2*g, k5*f f>1
% x([27 28 29]) = [0 0 1]; g = 10^x(1); fvc = ones(1,12); fvc(1) = g; x(30:41) = fvc; x(42) = 0; % corrected detailed balance; k2b = g*k2;
fvc = ones(1,12); fvc([1 3 7:12]) = 10.^[x(1) x(26) x(27:32)]./10.^x([12 10 14:19]); x(30:41) = fvc; x(42) = 0; % corrected detailed balance; k2b = g*k2;
x([27 28 29]) = [0 0 1];
% define all mutant parameters:
y = x; y(27) = 1; y(9) = -Inf; % -mgrB or mgrB constitutive
z = x; z([27 7]) = [1 0]; z(9) = -Inf; %-mgrB;-autoreg
z2 = x; z2(7) = 0; %-autoreg
% z3 = x; z3(14:19) = x(14:19)+k; z4 = z3; z4(27) = 1; % z3:-phosphatase; %z4 = -phosphatase;-mgrB
% z3(7) = 0; z4(7) = 0; %deltaautoregphoqt281r
% x: WT; % y: -mgrB; % z: -mgrB,-autoreg; % z2: -autoreg; % w: constitutive; % z3:-phosphatase; %z4 = -phosphatase;-mgrB

tss = [0 15]*3600;
x(29) = 1;
[~, X1] =ode15s(@phopq_0619_t, [0 12*3600], X0,{},x); % steady state at 1mM
x(29) = 30;
[~, Y30] =ode15s(@phopq_0619_t, tss, X1(end,:),{},x);
x(29) = 10;
[~, Y10]=ode15s(@phopq_0619_t, tss, X1(end,:),{},x);
x(29) = 3;
[~, Y3]=ode15s(@phopq_0619_t, tss, X1(end,:),{},x);
x(29) = 1;
[~, Y1]=ode15s(@phopq_0619_t, tss, X1(end,:),{},x);
x(29) = 0.3;
[~, Y03]=ode15s(@phopq_0619_t, tss, X1(end,:),{},x);
x(29) = 0.1;
[~, Y01]=ode15s(@phopq_0619_t, tss, X1(end,:),{},x);
x(29) = 0.03;
[~, Y003]=ode15s(@phopq_0619_t, tss, Y30(end,:),{},x);
% very low Mg steady state:
x(29) = 10^-4.5;
[~, Ymax]=ode15s(@phopq_0619_t, tss, Y30(end,:),{},x);
yfp_ss = [Y003(end,11)/Y003(end,19) Y01(end,11)/Y01(end,19) Y03(end,11)/Y03(end,19) Y1(end,11)/Y1(end,19) Y3(end,11)/Y3(end,19) Y10(end,11)/Y10(end,19) Y30(end,11)/Y30(end,19)]/(Y30(end,11)/Y30(end,19));
phopq_dr(:,1) = [0.03 0.1 0.3 1 3 10 30]; % Mg (mM)
phopq_dr(:,2) = [0.58 0.58 0.56 0.5 0.39 0.19 0.1]; % 
e_dy(11) = sum((yfp_ss' - phopq_dr(:,2)/phopq_dr(end,2)).^2); %/(norm(1 - phopq_dr(:,2)/phopq_dr(end,2),2)^2))*25;%/norm(phopq_dr(:,2)/phopq_dr(end,2),2);
e_dy(11) = e_dy(11)+((Ymax(end,11)/Ymax(end,19))/(Y30(end,11)/Y30(end,19))<6*(phopq_dr(1,2)/phopq_dr(end,2)))*25;
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
[~, X1] =ode15s(@phopq_0619_t, [0 8*3600], X0,{},x);
x(29) = 50;
[~, Y50] =ode15s(@phopq_0619_t, tspan, X1(end,:),{},x);
x(29) = 0.01;
[~, Y001]=ode15s(@phopq_0619_t, yfp_1s(:,1)*60, Y50(end,:),{},x); % for PmgrB
[~, Y001_pq]=ode15s(@phopq_0619_t, yfp_9s(:,1)*60, Y50(end,:),{},x); % for Pphopq
% 2 --> 0.01
x(29) = 2;
[~, Y2] =ode15s(@phopq_0619_t, [0 11]*3600 , X0,{},x);
x(29) = 0.01;
[~, Y2_001]=ode15s(@phopq_0619_t, yfp_2s(:,1)*60, Y2(end,:),{},x);
% 50 --> 2
x(29) = 2;
[~, Y50_2]=ode15s(@phopq_0619_t, yfp_5s(:,1)*60, Y50(end,:),{},x);
% 50 --> 10
x(29) = 10;
[~, Y10]=ode15s(@phopq_0619_t, yfp_6s(:,1)*60, Y50(end,:),{},x);
% 4 errors for 4 pmgrb-yfp output
e_dy(1) = norm((Y001(:,11)./Y001(:,19))/(Y50(end,11)/Y50(end,19)) - yfp_1s(:,2)/yfp_1s(1,2)); %/norm(yfp_1s(:,2)/yfp_1s(1,2));
e_dy(2) = norm((Y2_001(:,11)./Y2_001(:,19))/(Y50(end,11)/Y50(end,19)) - yfp_2s(:,2)/yfp_1s(1,2));%/norm(yfp_2s(:,2)/yfp_1s(1,2));
e_dy(5) = norm((Y50_2(:,11)./Y50_2(:,19))/(Y50(end,11)/Y50(end,19)) - yfp_5s(:,2)/yfp_1s(1,2));%/norm(yfp_5s(:,2)/yfp_1s(1,2));
e_dy(6) = norm((Y10(:,11)./Y10(:,19))/(Y50(end,11)/Y50(end,19)) - yfp_6s(:,2)/yfp_1s(1,2));%/norm(yfp_6s(:,2)/yfp_1s(1,2));
% one for pphopq-yfp output
e_dy(9) = norm((Y001_pq(:,12)./Y001_pq(:,19))/(Y50(end,11)/Y50(end,19)) - yfp_9s(:,2)/yfp_1s(1,2));%/norm(yfp_9s(:,2)/yfp_1s(1,2));
%% Dynamic simulation: delta mgrB
y(29) = 2;
[~, X1] =ode15s(@phopq_0619_t, [0 8*3600], X0del,{},y);
y(29) = 50;
[~, X_delmgr] =ode15s(@phopq_0619_t, tspan, X1(end,:),{},y);
y(29) = 0.01;
[~, Y_delmgr]=ode15s(@phopq_0619_t, yfp_3s(:,1)*60, X_delmgr(end,:),{},y);
[~, Y_delmgr_pq]=ode15s(@phopq_0619_t, yfp_10s(:,1)*60, X_delmgr(end,:),{},y);
e_dy(3) = norm((Y_delmgr(:,11)./Y_delmgr(:,19))/(Y50(end,11)/Y50(end,19)) - yfp_3s(:,2)/yfp_1s(1,2));%/norm(yfp_3s(:,2)/yfp_1s(1,2));
e_dy(10) = norm((Y_delmgr_pq(:,12)./Y_delmgr_pq(:,19))/(Y50(end,11)/Y50(end,19)) - yfp_10s(:,2)/yfp_1s(1,2));%/norm(yfp_10s(:,2)/yfp_1s(1,2)); % error for del mgrB Phopq-yfp

%% -autoreg
z2(29) = 2;
[~, X1] =ode15s(@phopq_0619_t, [0 8*3600], X0,{},z2);
z2(29) = 50;
[~, X_delautoreg] =ode15s(@phopq_0619_t, tspan, X1(end,:),{},z2);
z2(29) = 0.01;
[~, Y_delautoreg]=ode15s(@phopq_0619_t, yfp_7s(:,1)*60, X_delautoreg(end,:),{},z2);
e_dy(7) = norm((Y_delautoreg(:,11)./Y_delautoreg(:,19))/(Y50(end,11)/Y50(end,19)) - yfp_7s(:,2)/yfp_1s(1,2));%/norm(yfp_7s(:,2)/yfp_1s(1,2));
%% -autoreg;delta mgrB
% mgrB, MgrB
z(29) = 2;
[~, X1] =ode15s(@phopq_0619_t, [0 8*3600], X0del,{},z);
z(29) = 50;
[~, X_delmgrdelautoreg] =ode15s(@phopq_0619_t,tspan, X1(end,:),{},z);
z(29) = 0.01;
[~, Y_delmgrdelautoreg]=ode15s(@phopq_0619_t, yfp_8s(:,1)*60, X_delmgrdelautoreg(end,:),{},z);
e_dy(8) = norm((Y_delmgrdelautoreg(:,11)./Y_delmgrdelautoreg(:,19))/(Y50(end,11)/Y50(end,19)) - yfp_8s(:,2)/yfp_1s(1,2));%/norm(yfp_8s(:,2)/yfp_1s(1,2));
%% PhoQ(T281R) phosphatase mutant
% % z3(29) = 2;
% % [~, X1mono] =ode15s(@phopq_0619_t, [0 8*3600], X0,{},z3);
% % z3(29) = 50;
% % [~, X_mono] =ode15s(@phopq_0619_t, [0 3*3600], X1mono(end,:),{},z3);
% % z3(29) = 0.01;
% % [~, Y_mono]=ode15s(@phopq_0619_t, yfp_11s(:,1)*60, X_mono(end,:),{},z3);
% % % PhoQ(t281R) delta mgrB
% % z4(29) = 2;
% % [~, X1monodelmgr] =ode15s(@phopq_0619_t, [0 8*3600], X0del,{},z4);
% % z4(29) = 50;
% % [~, X_monodelmgr] =ode15s(@phopq_0619_t, [0 3*3600], X1monodelmgr(end,:),{},z4);
% % z4(29) = 0.01;
% % [~, Y_monodelmgr]=ode15s(@phopq_0619_t, yfp_12s(:,1)*60, X_monodelmgr(end,:),{},z4);
% % e_dy(4) = norm((Y_monodelmgr(:,11)./Y_monodelmgr(:,19))/(Y50(1,11)/Y50(1,19)) - yfp_12s(:,2)/yfp_1s(1,2));
% % e_dy(12) = norm((Y_mono(:,11)./Y_mono(:,19))/(Y50(1,11)/Y50(1,19)) - yfp_11s(:,2)/yfp_1s(1,2));
% % %% Qualitative error for inducible phoPQ(normal vs ph mutant)
% % % switch to inducible phopq; set PhoQ(T281R) phosphorylation/dephosphorylation parameters
% % x(42) = 0; z3(42) = 0; % phoPQ not inducible; to compute wild type levels to normalize to
% % x(29) = 1; z3(29) = 1;
% % % wild type at 1mM for normalization
% % [~, Xwt] =ode15s(@phopq_0619_t, [0 12*3600], X0,{},x); normlz = [Xwt(end,11)/Xwt(end,19) sum(Xwt(end,[1 2 8 9 16 18]),2)];
% % a = Xwt(end,2); f1 = 10^x(7); K1 = 10^x(8); wtind = (1+f1*(a/K1)^2)/(1+(a/K1)^2);
% % x(42) = 1; z3(42) = 1;% phopq induced
% % % simulate over range of induction rates
% % indrange = [0.01 0.1 0.25 0.75 1.5 2 2.5 3 3.5 4];
% % pp = zeros(length(indrange),4);
% % x(5) = log10(wtind*10^x(6)*0.47*indrange(1)); z3(5) = log10(wtind*10^x(6)*0.47*indrange(1));
% % [~, X1] =ode15s(@phopq_0619_t, [0 12*3600], X0,{},x);
% % [~, X1mono] =ode15s(@phopq_0619_t, [0 12*3600], X0,{},z3);
% % for i = 1:length(indrange)
% %     x(5) = log10(wtind*10^x(6)*0.47*indrange(i)); z3(5) = log10(wtind*10^x(6)*0.47*indrange(i));
% %     [~, Xdr] = ode15s(@phopq_0619_t, [0 3.5*3600], X1(end,:),{},x);
% %     [~, Xdrmono] = ode15s(@phopq_0619_t, [0 3.5*3600], X1mono(end,:),{},z3);
% %     pp(i,:) = [Xdr(end,11)/Xdr(end,19) sum(Xdr(end,[1 2 8 9 16 18]),2) Xdrmono(end,11)/Xdrmono(end,19) sum(Xdrmono(end,[1 2 8 9 16 18]),2)];
% % end
% % % define error for this i/o operation
% % idx = [4 5 9];
% % err_qual = zeros(1,3);
% % if pp(idx(1),1)/pp(idx(1),3) < 3
% %     err_qual(1) = 25;
% % end
% %     err_qual(2) = 25*(pp(idx(2),1)/pp(idx(2),3) - 1);
% % if pp(idx(3),1)/pp(idx(3),3)>= 0.3
% %     err_qual(3) = 25;
% % end
% % x(42) = 0; z3(42) = 0;
%% Total error
err = sum(([e_dy]).^2);
end