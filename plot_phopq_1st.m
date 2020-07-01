clear all
%% parameters
% 1. from full model
%load('0211_k1k2combo.mat'); x = Solution(2,:); x([27 28 29]) = [0 0 1]; x(30) = 0; % phoPQ not inducible
%x([12 10 11 22]) = log10([0.01 1e3 0.5 2]); % input k2(mg)
% 2. PSO 11252019: k2 = f(mg); k2b = g*k2;
load('11252019_phopq_1st'); x = Solution(1,:); x([27 28 29]) = [0 0 1]; x(30) = 0; 
% 3. PSO 12022019: k2,k5 = f(mg); k2b = g*k2; k2 INVERSE with stimulus.
% load('12022019_phopq_1st_k2k5inp'); x = Solution(3,:); x([27 28 29]) = [0 0 1]; x(30) = 0;
% 4. PSO 12172019: k2 = f(mg); k2,k3,km3,k4,k5,km5,k6 all affected by MgrB
% load('12172019_phopq_1st_all'); x = Solution(1,:); k=x(26:31); x(31:36) = k; x([27 28 29]) = [0 0 1]; x(30) = 0; 
% define all mutant parameters:
y = x; y(27) = 1; z = x; z([27 7]) = [1 0]; z2 = x; z2(7) = 0; w = x; w(28) = 1;
% z3(7) = 0; z4(7) = 0; %deltaautoregphoqt281r
% x: WT; % y: -mgrB; %z2: -autoreg; %z: -mgrB,-autoreg; % w: constitutive; % z3:-phosphatase; %z4 = -phosphatase;-mgrB
% initial vectors for +/- mgrB
X0 = zeros(1,17); X0del = X0; % initial condition 0 for mgrB deletion

%% dose response
mgrange = 10.^(-4.5:0.1:1.7); % dosage range of mg
pp = zeros(length(mgrange),3);
x(29) = 1; % mg = 1mM
[~, X1] =ode15s(@phopq_1st, [0 8*3600], X0,{},x);
[~, X2] =ode15s(@phopq_1st, [0 8*3600], X0del,{},y); 
[~, X3] =ode15s(@phopq_1st, [0 8*3600], X0,{},z2);
trange = [0 15*3600];
for i = 1:length(mgrange)
    x(29) = mgrange(i); y(29) = mgrange(i); z2(29) = mgrange(i);
    [~, Xdr] = ode15s(@phopq_1st,trange , X1(end,:),{},x);
    [~, Xdrdelmgr] =ode15s(@phopq_1st, trange, X2(end,:),{},y);
    [~, Xdrdelautoreg] =ode15s(@phopq_1st, trange, X3(end,:),{},z2);
    pp(i,:) = [Xdr(end,15)/Xdr(end,17) Xdrdelmgr(end,15)/Xdrdelmgr(end,17) Xdrdelautoreg(end,15)/Xdrdelautoreg(end,17)];
    pp2(i,:) = [Xdr(end,16)/Xdr(end,17) Xdrdelmgr(end,16)/Xdrdelmgr(end,17) Xdrdelautoreg(end,16)/Xdrdelautoreg(end,17)];
    rrp(i,:) = [Xdr(end,2) Xdrdelmgr(end,2) Xdrdelautoreg(end,2)];
end
idx = length(mgrange);%find(mgrange>=30,1,'first');
figure(1); subplot(2,3,6); semilogx(mgrange, pp(:,1)/pp(idx,1),'b','linewidth',1.5); hold on; semilogx(mgrange, pp(:,2)/pp(idx,1),'g');
semilogx(mgrange, pp(:,3)/pp(idx,1),'k');
phopq_dr(:,1) = [0.03 0.1 0.3 1 3 10 30]; % Mg (mM)
phopq_dr(:,2) = [0.58 0.58 0.56 0.5 0.39 0.19 0.1]; % 
semilogx(phopq_dr(:,1),phopq_dr(:,2)/phopq_dr(end,2),'s')
xlabel('Signal, k_{-1} (representing [Mg^{2+}])'); ylabel('Normalized YFP:CFP');
xlim([mgrange(1) mgrange(end)])
legend('WT-sim','WT-expt') %,'\DeltamgrB','\Deltaautoreg'
title('Dose-response');
set(gca,'fontsize',14,'linewidth',1)

%% 1. WT
% from Salazar et al: Cells were grown overnight in 2mM.
% diluted and shifted to 50mM for 2-3h
% washed and shifted to 2 or 0.01mM and imaged for 4h
% 50 mM --> 0.01 mM
x(29) = 2;
[~, X1] =ode15s(@phopq_1st, [0 8*3600], X0,{},x);
x(29) = 50;
[~, X] =ode15s(@phopq_1st, [0 3*3600], X1(end,:),{},x);
x(29) = 0.01;
[t, Y]=ode15s(@phopq_1st, [0 20]*3600, X(end,:),{},x);
% 2mM --> 0.01 mM
x(29) = 2;
[~, X_lh2] =ode15s(@phopq_1st, [0 3*3600], X1(end,:),{},x);
x(29) = 0.01;
[t_lh2, Y_lh2]=ode15s(@phopq_1st, [0 20]*3600, X_lh2(end,:),{},x);
% 50mM --> 10mM
x(29) = 10;
[t_lh3, Y_lh3]=ode15s(@phopq_1st, [0 20]*3600, X(end,:),{},x);
% 50mM --> 2mM
x(29) = 2;
[t_lh4, Y_lh4]=ode15s(@phopq_1st, [0 20]*3600, X(end,:),{},x);
%% 2. delta mgrB, w/ autoreg
y(29) = 2;
[~, X1delmgr] =ode15s(@phopq_1st, [0 8*3600], X0del,{},y);
y(29) = 50;
[~, X_delmgr] =ode15s(@phopq_1st, [0 3*3600], X1delmgr(end,:),{},y);
y(29) = 0.01;
[t_delmgr, Y_delmgr]=ode15s(@phopq_1st, [0 20]*3600, X_delmgr(end,:),{},y);
%% 3. delta mgrB, w/o autoreg
z(29) = 2;
[~, X1delmgrdelautoreg] =ode15s(@phopq_1st, [0 8*3600], X0del,{},z);
z(29) = 50;
[~, X_delmgrdelautoreg] =ode15s(@phopq_1st, [0 3*3600], X1delmgrdelautoreg(end,:),{},z);
z(29) = 0.01;
[t_delmgrdelautoreg, Y_delmgrdelautoreg]=ode15s(@phopq_1st, [0 20]*3600, X_delmgrdelautoreg(end,:),{},z);
%% 4. del autoreg only
z2(29) = 2;
[~, X1delautoreg] =ode15s(@phopq_1st, [0 8*3600], X0,{},z2);
z2(29) = 50;
[~, X_delautoreg] =ode15s(@phopq_1st, [0 3*3600], X1delautoreg(end,:),{},z2);
z2(29) = 0.01;
[t_delautoreg, Y_delautoreg]=ode15s(@phopq_1st, [0 20]*3600, X_delautoreg(end,:),{},z2);
%% 5. Constitutive expression of mgrB
w(9) = w(6)+0.85; %10xkbtpn2
w(29) = 2;
[~, X1const] =ode15s(@phopq_1st, [0 8*3600], X0,{},w);
w(29) = 50;
[~, X_const] =ode15s(@phopq_1st, [0 3*3600], X1const(end,:),{},w);
w(29) = 0.01;
[t_const, Y_const]=ode15s(@phopq_1st, [0 20]*3600, X_const(end,:),{},w);

%% 6. Induction response
zi = x; zi(30) = 1; % inducible phopq;
zi(28) = 1; zi(9) = x(6)+0.3; % constitutive mgrB expression rate
% zi(4) = -2; % translation rate of mgrB to simulate Figure s3
zi(29) = 1; % set signal
trx = linspace(-7.5,-3,40); % set range of phoPQ transcription rates
for j = 1:length(trx)
zi(5) = trx(j);
[~, Xind] =ode15s(@phopq_1st, [0 50*3600], X0,{},zi);
rind(j,:) = [Xind(end,15)/Xind(end,17) sum(Xind(end,[1 2 8 9 10 11]))]; %Xind(end,2) 
end

x(29) = 1; [~, X_norm] =ode15s(@phopq_1st, [0 40*3600], X0,{},x);
r_norm = [X_norm(end,15)/X_norm(end,17) sum(X_norm(end,[1 2 8 9 10 11]))];

figure(77); plot(rind(:,2)/r_norm(2),rind(:,1)/r_norm(1)); hold on
xlabel('{\it phoPphoQ} induction ([PhoP]/[PhoP]_{WT})')
ylabel('Normalized P_{mgrB} output')

%% plotting
load('yfp data sets.mat')
yfp_3s = yfp_datasets(:,:,1);yfp_4s = yfp_datasets(:,:,2);
yfp_5s = yfp_datasets(:,:,3);yfp_6s = yfp_datasets(:,:,4);
yfp_7s = yfp_datasets(:,:,5);yfp_8s = yfp_datasets(:,:,6);
yfp_9s = yfp_datasets(:,:,7);yfp_10s = yfp_datasets(:,:,8);
yfp_11s = yfp_datasets(:,:,9);yfp_12s = yfp_datasets(:,:,10);
yfp_1s = yfp_datasets(:,:,11);yfp_2s = yfp_datasets(:,:,12);

figure(1);% PmgrB-yfp outputs
   subplot(2,3,1) % 50 --> 0.01; WT, delmgr
plot(t/60,(Y(:,15)./Y(:,17))/(Y(1,15)/Y(1,17)),'b-'); hold on;
plot(t_delmgr/60, (Y_delmgr(:,15)./Y_delmgr(:,17))/(Y(1,15)/Y(1,17)),'m-'); %Delta mgrB
legend('WT','\DeltamgrB')
    % delmgrB
plot(yfp_1s(:,1),yfp_1s(:,2)/yfp_1s(1,2),'--'); hold on
plot(yfp_3s(:,1),yfp_3s(:,2)/yfp_1s(1,2),'m--'); % expt delta mgrB
ylabel('[YFP:CFP]/[YFP:CFP]_{WT,50mM}')
xlabel('time(mins)');
xlim([0 250]); ylim([0 45])
title('50 \rightarrow 0.01 mM : P_{mgrB}')
    subplot(2,3,2) 
     % 2 --> 0.01
plot(t_lh2/60, (Y_lh2(:,15)./Y_lh2(:,17))/(Y(1,15)/Y(1,17)),'r-'); hold on;
     % 50 --> 2
plot(t_lh4/60, (Y_lh4(:,15)./Y_lh4(:,17))/(Y(1,15)/Y(1,17)),'g-'); 
     % 50 --> 10
plot(t_lh3/60, (Y_lh3(:,15)./Y_lh3(:,17))/(Y(1,15)/Y(1,17)),'k-'); 
legend('2\rightarrow0.01','50\rightarrow2','50\rightarrow10')
plot(yfp_5s(:,1),yfp_5s(:,2)/yfp_1s(1,2),'g--'); hold on
plot(yfp_6s(:,1),yfp_6s(:,2)/yfp_1s(1,2),'k--'); hold on
plot(yfp_2s(:,1),yfp_2s(:,2)/yfp_1s(1,2),'r--'); hold on
xlabel('time(mins)');
axis([0 250 0 20])
title('WT: P_{mgrB}')
    subplot(2,3,3) % constitutive mgrB
plot(yfp_4s(:,1), yfp_4s(:,2)/yfp_1s(1,2),'--'); hold on
plot(t_const/60, (Y_const(:,15)./Y_const(:,17))/(Y(1,15)/Y(1,17)),'b-')
xlim([0 250]); ylim([0 20])
title('constitutive mgrB:P_{mgrB}')
xlabel('time(mins)');
    subplot(2,3,4) % phopq-yfp outputs
plot(t/60, (Y(:,16)./Y(:,17))/(Y(1,15)/Y(1,17)),'b-'); hold on
    % Pphopq-delmgr
plot(t_delmgr/60, (Y_delmgr(:,16)./Y_delmgr(:,17))/(Y(1,15)/Y(1,17)),'m-'); hold on
plot(yfp_9s(:,1), yfp_9s(:,2)/yfp_1s(1,2),'--')
plot(yfp_10s(:,1), yfp_10s(:,2)/yfp_1s(1,2),'m--')
ylabel('[YFP:CFP]/[YFP:CFP]_{WT,50mM}')
xlabel('time(mins)');
title('50-->0.01 mM : P_{phoPQ}')
xlim([0 250]); ylim([0 12])

    subplot(2,3,5)
plot(t_delmgrdelautoreg/60, (Y_delmgrdelautoreg(:,15)./Y_delmgrdelautoreg(:,17))/(Y(1,15)/Y(1,17)),'m-'); hold on % Delta mgrB 
                                                  % delta autoreg
plot(t_delautoreg/60, (Y_delautoreg(:,15)./Y_delautoreg(:,17))/(Y(1,15)/Y(1,17)),'k-'); %delta autoreg
plot(yfp_8s(:,1),yfp_8s(:,2)/yfp_1s(1,2),'m--'); % delmgrBdelautoreg
plot(yfp_7s(:,1),yfp_7s(:,2)/yfp_1s(1,2),'k--'); % delautoreg
%legend('\DeltamgrB\Deltaautoreg','\Deltaautoreg')%,'\DeltamgrB')
ylabel('YFP:CFP (fold)'); xlabel('time (mins)');
xlim([0 250]); %ylim([0 45]);
title('\Deltaautoreg: P_{mgrB}')