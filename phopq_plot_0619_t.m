clear all
%% parameters
% single parameter:
    % 1. k2
% load('09292019_k2_0619_t'); x = Solution(2,:); k = 0; x([27 28 29]) = [0 0 1]; f = 10^x(26); g = 10^x(1); fvc = ones(1,12); fvc(1) = g; x(30:41) = fvc;% corrected detailed balance(all equal assumption); k2b = g*k2;
% load('09292019_k2_0619_t2'); x = Solution(2,:); k = 0; x([27 28 29]) = [0 0 1]; f = 10^x(26); g = 10^x(1); fvc = ones(1,12); fvc(1) = g; x(30:41) = fvc;% corrected detailed balance(all equal assumption); k2b = g*k2;
% load('0214_k2.mat'); x = Solution(1,:); k=[-1.3 0 0 0 +1.3 -8]; x([27 28 29]) = [0 0 1]; f = 10^x(26); g = 10^x(1); fvc = ones(1,12); fvc(1) = g; x(30:41) = fvc;% corrected detailed balance(all equal assumption); k2b = g*k2;
    % 2. k3
% load('0214_KMT.mat'); x = Solution(1,:); k=[-1 0 0 -1 +2 -6]; x([27 28 29]) = [0 0 1]; f = 10^x(26); g = 10^x(1); fvc = ones(1,12); fvc(7) = f; x(30:41) = fvc; % corrected detailed balance(all equal assumption); k3b = f*k3;
    % 3. k5
% load('0214_KMP.mat'); x = Solution(1,:); k=0; x([27 28 29]) = [0 0 1];  f = 10^x(26); g = 10^x(1); fvc = ones(1,12); fvc(10) = 1/f; x(30:41) = fvc; % corrected detailed balance(all equal assumption); %k5b = k5/f;
    % 5. km5
% load('0214_km5.mat'); x = Solution(1,:); k=0; x([27 28 29]) = [0 0 1]; f = 10^x(26); g = 10^x(1); fvc = ones(1,12); fvc(11) = f; x(30:41) = fvc; % corrected detailed balance(all equal assumption); %km5b = km5*f;
% combos:
    % 1. k1 k2
% load('0214_k1k2combo.mat'); x = Solution(1,:); k=0; x([27 28 29]) = [0 0 1]; f = 10^x(26); g = 10^x(1); fvc = ones(1,12); fvc([1 3 6]) = [g f f]; x(30:41) = fvc; % corrected detailed balance; k1b = f*k1; kd1 = kd2*f; k2b = g*k2;
load('0211_k1k2combo.mat'); x = Solution(1,:); k=0; x([27 28 29]) = [0 0 1]; f = 10^x(26); g = 10^x(1); fvc = ones(1,12); fvc([1 3 6]) = [g f f]; x(30:41) = fvc; % corrected detailed balance; k1b = f*k1; kd1 = kd2*f; k2b = g*k2;
    % 2. k2 k3
% load('0214_k2KMTcombo.mat'); x = Solution(1,:); k=0; x([27 28 29]) = [0 0 1]; f = 10^x(26); g = 10^x(1); fvc = ones(1,12); fvc([1 7]) = [g f]; x(30:41) = fvc; % corrected detailed balance; k2b = g*k2; k3b = f*k3;
    % 3. k2 k5
% load('05242019_doseresp_k2k5.mat'); x = Solution(1,:); x([27 28 29]) = [0 0 1]; f = 10^x(26); g = 10^x(1); fvc = ones(1,12); fvc([1 10]) = [g f]; x(30:41) = fvc; k = 0; % k2*g, k5*f f>1
    % 4. k1k2k3km3k4k5km5k6
% load('12082019_all_0619t.mat'); x = Solution(4,:); fvc = ones(1,12); fvc([1 3 7:12]) = 10.^[x(1) x(26) x(27:32)]./10.^x([12 10 14:19]); 
% x(30:41) = fvc; x(42) = 0; % corrected detailed balance; k2b = g*k2;
% x([27 28 29]) = [0 0 1]; 

% % % % miscellaneous - test if needed, else archive.
% k1 k2 but nonlinear input function
% load('0228_k1k2_nonlinMg.mat'); x = Solution(1,:); k=0; x([27 28 29]) = [0 0 1]; f = 10^x(26); g = 10^x(1); fvc = ones(1,12); fvc([1 3 6]) = [g f f]; x(30:41) = fvc; % corrected detailed balance; k1b = f*k1; kb1 = kb2/f; k2b = g*k2;
% k1 k2, phoq mutant
% load('0321phoqt281r_initialmatrix.mat'); x = Solution(1,:); x(8) = log10(0.6); k = x(27:32);x([27 28 29]) = [0 0 1]; f = 10^x(26); g = 10^x(1); fvc = ones(1,12); fvc([1 3 6]) = [g f f]; x(30:41) = fvc; % corrected detailed balance; k1b = f*k1; kb1 = kb2/f; k2b = g*k2;
% load('0331_k1k2_steadystate'); x = Solution(1,:); x([7 8]) = [log10(25) log10(1)]; k = x(27:32);x([27 28 29]) = [0 0 1]; f = 10^x(26); g = 10^x(1); fvc = ones(1,12); fvc([1 3 6]) = [g f f]; x(30:41) = fvc; % corrected detailed balance; k1b = f*k1; kb1 = kb2/f; k2b = g*k2;

x(42) = 0; % phoPQ not inducible

% define all mutant parameters:
y = x; y(27) = 1; y(9) = -Inf; % -mgrB or mgrB constitutive
z = x; z([27 7]) = [1 0]; z(9) = -Inf; %-mgrB;-autoreg
z2 = x; z2(7) = 0; %-autoreg
z3 = x; z3(14:19) = x(14:19)+k; z4 = z3; z4(27) = 1; % z3:-phosphatase; %z4 = -phosphatase;-mgrB
% z3(7) = 0; z4(7) = 0; %deltaautoregphoqt281r

% initial vectors for +/- mgrB
X0 = zeros(1,19); X0del = X0;
%% dose response
mgrange = 10.^(-4.5:0.05:1.7); % dosage range of mg
pp = zeros(length(mgrange),3);
x(29) = 1; % mg = 1mM
[~, X1] =ode15s(@phopq_0619_t, [0 20*3600], X0,{},x);
[~, X2] =ode15s(@phopq_0619_t, [0 20*3600], X0del,{},y); 
[~, X3] =ode15s(@phopq_0619_t, [0 20*3600], X0,{},z2);
trange = [0 15*3600];
for i = 1:length(mgrange)
    x(29) = mgrange(i); y(29) = mgrange(i); z2(29) = mgrange(i);
    [~, Xdr] = ode15s(@phopq_0619_t,trange , X1(end,:),{},x);
    [~, Xdrdelmgr] =ode15s(@phopq_0619_t, trange, X2(end,:),{},y);
    [~, Xdrdelautoreg] =ode15s(@phopq_0619_t, trange, X3(end,:),{},z2);
    pp(i,:) = [Xdr(end,11)/Xdr(end,19) Xdrdelmgr(end,11)/Xdrdelmgr(end,19) Xdrdelautoreg(end,11)/Xdrdelautoreg(end,19)];
    pp2(i,:) = [Xdr(end,12)/Xdr(end,19) Xdrdelmgr(end,12)/Xdrdelmgr(end,19) Xdrdelautoreg(end,12)/Xdrdelautoreg(end,19)];
    rrp(i,:) = [Xdr(end,2) Xdrdelmgr(end,2) Xdrdelautoreg(end,2)];
end
idx = length(mgrange);%find(mgrange>=30,1,'first');
figure(1);subplot(2,3,6); semilogx(mgrange, pp(:,1)/pp(idx,1),'b'); hold on; semilogx(mgrange, pp(:,2)/pp(idx,1),'g');
semilogx(mgrange, pp(:,3)/pp(idx,1),'k');
phopq_dr(:,1) = [0.03 0.1 0.3 1 3 10 30]; % Mg (mM)
phopq_dr(:,2) = [0.58 0.58 0.56 0.5 0.39 0.19 0.1]; % 
semilogx(phopq_dr(:,1),phopq_dr(:,2)/phopq_dr(end,2),'s')
xlabel('Signal, k_{-1} (representing [Mg^{2+}])'); ylabel('Normalized YFP:CFP');
xlim([mgrange(1) mgrange(end)])
legend('WT','\DeltamgrB','\Deltaautoreg','WT-expt')
title('Dose-response');

%% sensitivity of rrp
% % % very crude way of calculating numerical derivative
% % ddmg = (rrp(2:end,:)-rrp(1:end-1,:))./(mgrange(2:end)-mgrange(1:end-1))';
% % loggain = ddmg.*mgrange(1:end-1)'./rrp(1:end-1,:);
% % figure(3); semilogx(mgrange(1:end-1),abs(loggain)); hold on
% % xlabel('Signal, k_{-1}/k_{-1}^0'); ylabel('|\partial log[PhoP~P]/\partial log(k_{-1})|');
% % legend('Wild type','No negative feedback, MgrB = 0','No autoregulation')

%% 1. WT
% from Salazar et al: Cells were grown overnight in 2mM.
% diluted and shifted to 50mM for 2-3h
% washed and shifted to 2 or 0.01mM and imaged for 4h
% 50 mM --> 0.01 mM
x(29) = 2;
[~, X1] =ode15s(@phopq_0619_t, [0 8*3600], X0,{},x);
x(29) = 50;
[~, X] =ode15s(@phopq_0619_t, [0 3*3600], X1(end,:),{},x);
x(29) = 0.01;
[t, Y]=ode15s(@phopq_0619_t, [0 20]*3600, X(end,:),{},x);
% 2mM --> 0.01 mM
x(29) = 2;
[~, X_lh2] =ode15s(@phopq_0619_t, [0 3*3600], X1(end,:),{},x);
x(29) = 0.01;
[t_lh2, Y_lh2]=ode15s(@phopq_0619_t, [0 20]*3600, X_lh2(end,:),{},x);
% 50mM --> 10mM
x(29) = 10;
[t_lh3, Y_lh3]=ode15s(@phopq_0619_t, [0 20]*3600, X(end,:),{},x);
% 50mM --> 2mM
x(29) = 2;
[t_lh4, Y_lh4]=ode15s(@phopq_0619_t, [0 20]*3600, X(end,:),{},x);
%% 2. delta mgrB, w/ autoreg
y(29) = 2;
[~, X1delmgr] =ode15s(@phopq_0619_t, [0 8*3600], X0del,{},y);
y(29) = 50;
[~, X_delmgr] =ode15s(@phopq_0619_t, [0 3*3600], X1delmgr(end,:),{},y);
y(29) = 0.01;
[t_delmgr, Y_delmgr]=ode15s(@phopq_0619_t, [0 20]*3600, X_delmgr(end,:),{},y);
%% 3. delta mgrB, w/o autoreg
z(29) = 2;
[~, X1delmgrdelautoreg] =ode15s(@phopq_0619_t, [0 8*3600], X0del,{},z);
z(29) = 50;
[~, X_delmgrdelautoreg] =ode15s(@phopq_0619_t, [0 3*3600], X1delmgrdelautoreg(end,:),{},z);
z(29) = 0.01;
[t_delmgrdelautoreg, Y_delmgrdelautoreg]=ode15s(@phopq_0619_t, [0 20]*3600, X_delmgrdelautoreg(end,:),{},z);
%% 4. del autoreg only
z2(29) = 2;
[~, X1delautoreg] =ode15s(@phopq_0619_t, [0 8*3600], X0,{},z2);
z2(29) = 50;
[~, X_delautoreg] =ode15s(@phopq_0619_t, [0 3*3600], X1delautoreg(end,:),{},z2);
z2(29) = 0.01;
[t_delautoreg, Y_delautoreg]=ode15s(@phopq_0619_t, [0 20]*3600, X_delautoreg(end,:),{},z2);
%% 5. Constitutive expression of mgrB
y(9) = x(6)+1; %10xkbtpn2
y(29) = 2;
[~, X1const] =ode15s(@phopq_0619_t, [0 8*3600], X0,{},y);
y(29) = 50;
[~, X_const] =ode15s(@phopq_0619_t, [0 3*3600], X1const(end,:),{},y);
y(29) = 0.01;
[t_const, Y_const]=ode15s(@phopq_0619_t, [0 20]*3600, X_const(end,:),{},y);
%% 6. Induction response at 1mM (WT)
xind = x;
xind(29) = 1; % signal level
xind(42) = 1;% phopq induced
xind(4) = -3; % varying basal MgrB, comment out if not needed
% simulate over range of induction rates
indrange = linspace(-7,-4,40);
pp = zeros(length(indrange),2);
for i = 1:length(indrange)
    xind(5) = indrange(i); 
    [~, Xdr] = ode15s(@phopq_0619_t, [0 20*3600], X0,{},xind);
    pp(i,:) = [Xdr(end,11)/Xdr(end,19) sum(Xdr(end,[1 2 8 9 16 18]),2)];
end
xind(42) = 0; 
% normalize to wild type phoPQ induction
x(29) = 1; [~, X_norm] =ode15s(@phopq_0619_t, [0 4*3600], X1(end,:),{},x);
ppnorm = [X_norm(end,11)/X_norm(end,19) sum(X_norm(end,[1 2 8 9 16 18]),2)];

figure(4); plot(pp(:,2)/ppnorm(2),pp(:,1)/ppnorm(1),'k'); hold on;
%% 6b. Induction response at 1mM (constitutive mgrB) For predictions
y(29) = 1; % signal level
y(9) = x(6) + 0.67; % induction level of mgrB = 10x wild type basal
y(42) = 1; % phopq induced
% simulate over range of induction rates
indrange = linspace(-7,-4,40);
pp_const = zeros(length(indrange),2);
for i = 1:length(indrange)
    y(5) = indrange(i); 
    [~, Xdr_const] = ode15s(@phopq_0619_t, [0 60*3600], X0,{},y);
    pp_const(i,:) = [Xdr_const(end,11)/Xdr_const(end,19) sum(Xdr_const(end,[1 2 8 9 16 18]),2)];
end
y(42) = 0;
figure(4); plot(pp_const(:,2)/ppnorm(2),pp_const(:,1)/ppnorm(1),'m'); hold on;
xlabel('{\it phoPphoQ} induction ([PhoP]/[PhoP]_{WT})')
ylabel('Normalized P_{mgrB} output')
legend('WT simulation','Expt','constitutive mgrB')
%% 7. PhoQ(T281R) - phosphatase lacking mutant
z3(17:19) = x(17:19)-5;
z4(17:19) = x(17:19)-5; z4(9) = -Inf;
z3(29) = 2;
[~, X1mono] =ode15s(@phopq_0619_t, [0 8*3600], X0,{},z3);
z3(29) = 50;
[~, X_mono] =ode15s(@phopq_0619_t, [0 3*3600], X1mono(end,:),{},z3);
z3(29) = 0.01;
[t_mono, Y_mono]=ode15s(@phopq_0619_t, [0 20]*3600, X_mono(end,:),{},z3);
% PhoQ(t281R) delta mgrB
z4(29) = 2;
[~, X1monodelmgr] =ode15s(@phopq_0619_t, [0 8*3600], X0del,{},z4);
z4(29) = 50;
[~, X_monodelmgr] =ode15s(@phopq_0619_t, [0 3*3600], X1monodelmgr(end,:),{},z4);
z4(29) = 0.01;
[t_monodelmgr, Y_monodelmgr]=ode15s(@phopq_0619_t, [0 20]*3600, X_monodelmgr(end,:),{},z4);

% % 1mM Mg +/-autoreg +/-PhoQ(phosphatase) barplot
% % z3(29) = 1; % phoqt281r +autoreg
% % [~, bar3] =ode15s(@phopq_0619_t, [0 12*3600], X0(end,:),{},z3);
% % z5(29) = 1; %phoqt281r -autoreg
% % [~, bar4] =ode15s(@phopq_0619_t, [0 12*3600], X0(end,:),{},z5);
% % x(29) = 1; %WT
% % [~, bar1] =ode15s(@phopq_0619_t, [0 12*3600], X0,{},x);
% % z2(29) = 1;
% % [~, bar2] =ode15s(@phopq_0619_t, [0 12*3600], X0,{},z2);
% % 
% % figure; bar([bar1(end,11)/bar1(end,19) bar2(end,11)/bar2(end,19) bar3(end,11)/bar3(end,19) bar4(end,11)/bar4(end,19)])
%% plotting all simulations against experimental data
load('yfp data sets.mat')
yfp_3s = yfp_datasets(:,:,1);yfp_4s = yfp_datasets(:,:,2);
yfp_5s = yfp_datasets(:,:,3);yfp_6s = yfp_datasets(:,:,4);
yfp_7s = yfp_datasets(:,:,5);yfp_8s = yfp_datasets(:,:,6);
yfp_9s = yfp_datasets(:,:,7);yfp_10s = yfp_datasets(:,:,8);
yfp_11s = yfp_datasets(:,:,9);yfp_12s = yfp_datasets(:,:,10);
yfp_1s = yfp_datasets(:,:,11);yfp_2s = yfp_datasets(:,:,12);

figure(1);% PmgrB-yfp outputs
   subplot(2,3,1) % 50 --> 0.01; WT, delmgr
plot(t/60,(Y(:,11)./Y(:,19))/(Y(1,11)/Y(1,19)),'b-'); hold on;
plot(t_delmgr/60, (Y_delmgr(:,11)./Y_delmgr(:,19))/(Y(1,11)/Y(1,19)),'m-'); %Delta mgrB
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
plot(t_lh2/60, (Y_lh2(:,11)./Y_lh2(:,19))/(Y(1,11)/Y(1,19)),'r-'); hold on;
     % 50 --> 2
plot(t_lh4/60, (Y_lh4(:,11)./Y_lh4(:,19))/(Y(1,11)/Y(1,19)),'g-'); 
     % 50 --> 10
plot(t_lh3/60, (Y_lh3(:,11)./Y_lh3(:,19))/(Y(1,11)/Y(1,19)),'k-'); 
legend('2\rightarrow0.01','50\rightarrow2','50\rightarrow10')
plot(yfp_5s(:,1),yfp_5s(:,2)/yfp_1s(1,2),'g--'); hold on
plot(yfp_6s(:,1),yfp_6s(:,2)/yfp_1s(1,2),'k--'); hold on
plot(yfp_2s(:,1),yfp_2s(:,2)/yfp_1s(1,2),'r--'); hold on
xlabel('time(mins)');
axis([0 250 0 20])
title('WT: P_{mgrB}')
    subplot(2,3,3) % constitutive mgrB
plot(yfp_4s(:,1), yfp_4s(:,2)/yfp_1s(1,2),'--'); hold on
plot(t_const/60, (Y_const(:,11)./Y_const(:,19))/(Y(1,11)/Y(1,19)),'b-')
xlim([0 250]); ylim([0 20])
title('constitutive mgrB:P_{mgrB}')
xlabel('time(mins)');
    subplot(2,3,4) % phopq-yfp outputs
plot(t/60, (Y(:,12)./Y(:,19))/(Y(1,11)/Y(1,19)),'b-'); hold on
    % Pphopq-delmgr
plot(t_delmgr/60, (Y_delmgr(:,12)./Y_delmgr(:,19))/(Y(1,11)/Y(1,19)),'m-'); hold on
plot(yfp_9s(:,1), yfp_9s(:,2)/yfp_1s(1,2),'--')
plot(yfp_10s(:,1), yfp_10s(:,2)/yfp_1s(1,2),'m--')
ylabel('[YFP:CFP]/[YFP:CFP]_{WT,50mM}')
xlabel('time(mins)');
title('50 \rightarrow 0.01 mM : P_{phoPQ}')
xlim([0 250]); ylim([0 12])

    subplot(2,3,5)
plot(t_delmgrdelautoreg/60, (Y_delmgrdelautoreg(:,11)./Y_delmgrdelautoreg(:,19))/(Y(1,11)/Y(1,19)),'m-'); hold on % Delta mgrB 
                                                  % delta autoreg
plot(t_delautoreg/60, (Y_delautoreg(:,11)./Y_delautoreg(:,19))/(Y(1,11)/Y(1,19)),'k-'); %delta autoreg
plot(yfp_8s(:,1),yfp_8s(:,2)/yfp_1s(1,2),'m--'); % delmgrBdelautoreg
plot(yfp_7s(:,1),yfp_7s(:,2)/yfp_1s(1,2),'k--'); % delautoreg
%legend('\DeltamgrB\Deltaautoreg','\Deltaautoreg')%,'\DeltamgrB')
ylabel('YFP:CFP (fold)'); xlabel('time (mins)');
xlim([0 250]); %ylim([0 45]);
title('\Deltaautoreg: P_{mgrB}')

% % figure(4); subplot(2,1,1) % PhoQ(T281R)
% % plot(yfp_11s(:,1), yfp_11s(:,2)/yfp_1s(1,2),'g--'); hold on
% % plot(yfp_12s(:,1), yfp_12s(:,2)/yfp_1s(1,2),'LineStyle','--','Color',[0 0.5 0]);
% % plot(t_mono/60, (Y_mono(:,11)./Y_mono(:,19))/(Y(1,11)/Y(1,19)),'g');
% % plot(t_monodelmgr/60, (Y_monodelmgr(:,11)./Y_monodelmgr(:,19))/(Y(1,11)/Y(1,19)),'Color',[0 0.5 0]);
% % xlim([0 400]);ylim([0 60]); ylabel('YFP (fold)'); xlabel('time (mins)');
% % title('PhoQ(T281R)+/-mgrB')

% % strg = 'kpdeg kmdeg ktlnA ktlnE kbtpn1 kbtpn2 f1  K1    ktlnB  k1  km1 k2  km2 k3 km3  k4 k5 km5 k6  K2     f2   kbtpn3 ktlnY kb kd k2b km2b k3b km3b k4b';
% % vars = strsplit(strg,' ');
% % figure
% % for i = 1:size(sol,2)-1
% %     if sol(end-1,i) ~= sol(end,i)
% %     subplot(5,6,i)
% %     hist(sol(1:end-2,i))
% %     set(gca,'XLim',sol([end-1 end],i))
% %     xlabel(vars{i})
% %     end
% % end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Pp vs PT and adaptation details
% % This section deals with saturation characteristics of PhoPQ-MgrB module
% % PhoP,PhoQ total remain constant, MgrB remains a function of PhoP~P
% % % clear all
% % % single parameter:
% % % load('0214_k2.mat'); x = Solution(1,:); k=[-1.3 0 0 0 +1.3 -8]; x([27 28 29]) = [0 0 1]; f = 10^x(26); g = 10^x(1); fvc = ones(1,12); fvc(1) = g; x(30:41) = fvc;% corrected detailed balance(all equal assumption); k2b = g*k2;
% % % load('0214_KMT.mat'); x = Solution(1,:); k=[-1 0 0 -1 +2 -6]; x([27 28 29]) = [0 0 1]; f = 10^x(26); g = 10^x(1); fvc = ones(1,12); fvc(7) = f; x(30:41) = fvc; % corrected detailed balance(all equal assumption); k3b = f*k3;
% % % load('0214_KMP.mat'); x = Solution(1,:); k=0; x([27 28 29]) = [0 0 1];  f = 10^x(26); g = 10^x(1); fvc = ones(1,12); fvc(10) = 1/f; x(30:41) = fvc; % corrected detailed balance(all equal assumption); %k5b = k5/f;
% % % load('0214_km5.mat'); x = Solution(1,:); k=0; x([27 28 29]) = [0 0 1]; f = 10^x(26); g = 10^x(1); fvc = ones(1,12); fvc(11) = f; x(30:41) = fvc; % corrected detailed balance(all equal assumption); %km5b = km5*f;
% % % combos:
% % load('0214_k1k2combo.mat'); x = Solution(1,:); k=0; x([27 28 29]) = [0 0 1]; f = 10^x(26); g = 10^x(1); fvc = ones(1,12); fvc([1 3 5]) = [g f 1/f]; x(30:41) = fvc; % corrected detailed balance; k1b = f*k1; kb1 = kb2/f; k2b = g*k2;
% % % load('0211_k1k2combo.mat'); x = Solution(1,:); k=0; x([27 28 29]) = [0 0 1]; f = 10^x(26); g = 10^x(1); fvc = ones(1,12); fvc([1 3 5]) = [g f 1/f]; x(30:41) = fvc; % corrected detailed balance; k1b = f*k1; kb1 = kb2/f; k2b = g*k2;
% % % load('0214_k2KMTcombo.mat'); x = Solution(1,:); k=0; x([27 28 29]) = [0 0 1]; f = 10^x(26); g = 10^x(1); fvc = ones(1,12); fvc([1 7]) = [g f]; x(30:41) = fvc; % corrected detailed balance; k2b = g*k2; k3b = f*k3;
% % % load('0228_k1k2_nonlinMg.mat'); x = Solution(1,:); k=0; x([27 28 29]) = [0 0 1]; f = 10^x(26); g = 10^x(1); fvc = ones(1,12); fvc([1 3 5]) = [g f 1/f]; x(30:41) = fvc; % corrected detailed balance; k1b = f*k1; kb1 = kb2/f; k2b = g*k2;
% % % load('0321phoqt281r_initialmatrix.mat'); x = Solution(1,:); k = x(27:32);x([27 28 29]) = [0 0 1]; f = 10^x(26); g = 10^x(1); fvc = ones(1,12); fvc([1 3 5]) = [g f 1/f]; x(30:41) = fvc; % corrected detailed balance; k1b = f*k1; kb1 = kb2/f; k2b = g*k2;
% % % load('0331_k1k2_steadystate'); x = Solution(2,:); k = x(27:32);x([27 28 29]) = [0 0 1]; f = 10^x(26); g = 10^x(1); fvc = ones(1,12); fvc([1 3 5]) = [g f 1/f]; x(30:41) = fvc; % corrected detailed balance; k1b = f*k1; kb1 = kb2/f; k2b = g*k2;
% % x(42) = 0; % phoPQ not inducible
% % % define all mutant parameters:
% % y = x; y(27) = 1; z = x; z([27 7]) = [1 0]; z2 = x; z2(7) = 0; w = x; w(28) = 1; z3 = x; z3(14:19) = x(14:19)+k; z4 = z3; z4(27) = 1; z5 = z3; z5(7) = 0;
% % % x: WT; % y: -mgrB; %z2: -autoreg; %z: -mgrB,-autoreg; % w: constitutive; % z3:-phosphatase; %z4 = -phosphatase;-mgrB
% % % initial vectors for +/- mgrB
% % X0 = zeros(1,19); X0([3 4]) = 0.4*10^x(6)/10^x(2); %initialize to ~ basal levels of mgrB, phoP/Q mRNA
% % X0del = X0; X0del(4) = 0; % initial condition 0 for mgrB deletion
% % y(21) = 0.8; y(11)=-0.5; % equivalent TCS by dose-response (0214_k1k2_combo#1)
% % x(29) = 1; x(42) = 1; z3(42) = 1; y(42) = 1; mg = 1;
% % mgrange = (-8:0.1:-2);
% % pp = zeros(length(mgrange),6);
% % x(5) = -8; z3(5) = -8; y(5) = -8;
% % [~, X1] =ode15s(@phopq_0619_t, [0 60*3600], X0,{},x);
% % [~, X1mono] =ode15s(@phopq_0619_t, [0 60*3600], X0,{},z3);
% % [~, X1delmgr] =ode15s(@phopq_0619_t, [0 60*3600], X0del,{},y);
% % x(29) = mg; y(29) = mg; z3(29) = mg;
% % kpdeg_eff = 3.1e-4;
% % for i = 1:length(mgrange)
% %     x(5) = mgrange(i); z3(5) = mgrange(i); y(5) = mgrange(i);
% %     [~, Xdr] = ode15s(@phopq_0619_t, [0 80*3600], X1(end,:),{},x);
% %     [~, Xdrmono] = ode15s(@phopq_0619_t, [0 80*3600], X1mono(end,:),{},z3);
% %     [~, Xdrdelmgr] = ode15s(@phopq_0619_t, [0 80*3600], X1delmgr(end,:),{},y);
% % % % %     pp(i,:) = [Xdr(end,11)/Xdr(end,19) sum(Xdr(end,[1 2 8 9 16 18]),2) Xdrmono(end,11)/Xdrmono(end,19) sum(Xdrmono(end,[1 2 8 9 16 18]),2) Xdrdelmgr(end,11)/Xdrdelmgr(end,19) sum(Xdrdelmgr(end,[1 2 8 9 16 18]),2)];
% %     pp(i,:) = [Xdr(end,2) sum(Xdr(end,[1 2 8 9 16 18]),2) Xdrmono(end,2) sum(Xdrmono(end,[1 2 8 9 16 18]),2) Xdrdelmgr(end,2) sum(Xdrdelmgr(end,[1 2 8 9 16 18]),2)];
% % % %     Btot(i,:) = [Xdr(end,4) Xdrmono(end,4) Xdrdelmgr(end,4)]*10^x(4)/kpdeg_eff;
% % end
% % figure; loglog(pp(:,2),pp(:,1),'k'); hold on; loglog(pp(:,6),pp(:,5),'c');
% % title({'input-output relation of PhoPQ phosphotransfer/phosphatase cycle';' under inducible \textit{phopq} promoter'},'Interpreter','latex')
% % legend('50mM, +mgrB','50mM, -mgrB','1mM, +mgrB','1mM, -mgrB');
% % xlabel('PhoP_{Total}'); ylabel('PhoP~P')
% % dydx = diff(pp(:,1))./diff(pp(:,2)); lg = pp(2:end,2).*(dydx./pp(2:end,1)); %lg(1) = [];
% % dydxdel = diff(pp(:,5))./diff(pp(:,6)); lgdel = pp(2:end,6).*dydxdel./pp(2:end,5); %lgdel(1) = [];
% % figure; subplot(2,1,2); semilogx(pp(2:end,2),lg,'k--',pp(2:end,4),lgdel,'c--'); xlabel('PhoP_T'); %ylabel('$\frac{\partial log([PhoP_P])}{\partial log([PhoP_T])}$','interpreter','latex','rotation',0,'fontsize',20)
% % legend('PhoP~P w/ NF','PhoP~P w/o NF','Log gain w/ NF','Log gain w/o NF')

% % x(5) = -4.7; [t1, int1] =ode15s(@phopq_0619_t, [0 60*3600], X0,{},x);
% % x(5) = -4.2; [t2, int2] =ode15s(@phopq_0619_t, [0 60*3600], int1(end,:),{},x);
% % loglog(sum(int1(end,[1 2 8 9 16 18])),int1(end,2),'ko','linewidth',1.8)
% % loglog(sum(int2(:,[1 2 8 9 16 18]),2),int2(:,2),'k','linewidth',1.8)
% % figure; subplot(3,1,1); plot(t2/3600, int2(:,5),'color','k'); hold on; ylabel('Kinase')
% % plot(t2/3600, g*int2(:,14),'color','k','linestyle','--'); legend('Q*','\lambdaQB*'); xlim([-0.1 2.5]);
% % subplot(3,1,2); plot(t2/3600, int2(:,6),'color','k'); hold on; ylabel('Phosphatase')
% % plot(t2/3600, int2(:,17),'color','k','linestyle','--'); legend('Q','QB'); xlim([-0.1 2.5])
% % subplot(3,1,3); plot(t2/3600, (int2(:,5)+g*int2(:,14))./(int2(:,6)+int2(:,17))); xlim([-0.1 2.5])
% % ylabel('(Q*+\lambdaQB*)/(Q+QB)'); xlabel('time (h)')
% % y(5) = -4.7; [t1, int1] =ode15s(@phopq_0619_t, [0 60*3600], X0del,{},y);
% % y(5) = -4.2; [t2, int2] =ode15s(@phopq_0619_t, [0 60*3600], int1(end,:),{},y);
% % loglog(sum(int1(end,[1 2 8 9 16 18])),int1(end,2),'ko','linewidth',1.8)
% % loglog(sum(int2(:,[1 2 8 9 16 18]),2),int2(:,2),'k','linewidth',1.8)
% % x(42) = 0; z3(42) = 0; y(42) = 0;
%% wide range dose response - on and off
% % clear all
% % % load('PSO0619_0729'); x = sol(1,:); % incorrect detailed balance
% % load('02112019_pso.mat'); x = Solution(1,:); % corrected detailed balance, f =/= 1
% % % load('0214_k2.mat'); x = Solution(1,:); % corrected detailed balance f = 1
% % % load('0214_KMT.mat'); x = Solution(1,:); % corrected detailed balance f = 1
% % % load('0214_KMP.mat'); x = Solution(1,:); % corrected detailed balance f = 1
% % x([27 28 29]) = [0 0 1];
% % y = x; y(27) = 1; z2 = x; z2(7) = 0;
% % z3 = x; z3(17) = -3; z(19) = -8;
% % 
% % mgrange = fliplr(logspace(-3,3,20));
% % pp = zeros(length(mgrange),4);
% % X0 = 0.01*ones(1,19);
% % x(29) = 1e3;
% % [~, X1] =ode15s(@phopq_0619_t, [0 80*3600], X0,{},x);
% % [~, X3] =ode15s(@phopq_0619_t, [0 80*3600], X0,{},z2);
% % [~, X4] =ode15s(@phopq_0619_t, [0 80*3600], X0,{},z3);
% % X0([4 10 14 15 16 17 18]) = 0;
% % [~, X2] =ode15s(@phopq_0619_t, [0 80*3600], X0,{},y); 
% % for i = 1:length(mgrange)
% %     x(29) = mgrange(i);y(29) = mgrange(i); z2(29) = mgrange(i); z3(29) = mgrange(i);
% %     [~, Xdr] = ode15s(@phopq_0619_t, [0 60*3600], X1(end,:),{},x);
% %     [~, Xdrdelmgr] =ode15s(@phopq_0619_t, [0 60*3600], X2(end,:),{},y);
% %     [~, Xdrdelautoreg] =ode15s(@phopq_0619_t, [0 60*3600], X3(end,:),{},z2);
% %     [~, Xph] = ode15s(@phopq_0619_t, [0 60*3600], X4(end,:),{},z3);
% %     pp(i,:) = [Xdr(end,11)/Xdr(end,19) Xdrdelmgr(end,11)/Xdrdelmgr(end,19) Xdrdelautoreg(end,11)/Xdrdelautoreg(end,19) Xph(end,11)/Xph(end,19)];
% % end
% % mgrange2 = fliplr(mgrange);
% % for i = 1:length(mgrange)
% %     x(29) = mgrange2(i);y(29) = mgrange2(i); z2(29) = mgrange2(i); z3(29) = mgrange2(i);
% %     [~, Xdr] = ode15s(@phopq_0619_t, [0 60*3600], Xdr(end,:),{},x);
% %     [~, Xdrdelmgr] =ode15s(@phopq_0619_t, [0 60*3600], Xdrdelmgr(end,:),{},y);
% %     [~, Xdrdelautoreg] =ode15s(@phopq_0619_t, [0 60*3600], Xdrdelautoreg(end,:),{},z2);
% %     [~, Xph] = ode15s(@phopq_0619_t, [0 60*3600], Xph(end,:),{},z3);
% % 
% %     ppoff(i,:) = [Xdr(end,11)/Xdr(end,19) Xdrdelmgr(end,11)/Xdrdelmgr(end,19) Xdrdelautoreg(end,11)/Xdrdelautoreg(end,19) Xph(end,11)/Xph(end,19)];
% % end
% % 
% % figure; loglog(mgrange, pp(:,1),mgrange2,ppoff(:,1)); hold on;
% % loglog(mgrange, pp(:,2),mgrange2,ppoff(:,2)); 
% % loglog(mgrange, pp(:,3),mgrange2,ppoff(:,3)); 
% % loglog(mgrange, pp(:,4),mgrange2,ppoff(:,4)); 
% % % dose response of YFP vs IPTG at any given Mg concentration
% % iptg = [(-7:0.05:-3) fliplr((-7:0.05:-3))];
% % yfp_ptac = zeros(size(iptg));
% % w(29) = 0.01;
% % for i = 1:length(iptg)
% % w(22) = iptg(i);
% % [~, X_const] =ode15s(@phopq_0619_t, [0 8*3600], X1const(end,:),{},w);
% % yfp_ptac(i) = X_const(end,11)/X_const(end,19);
% % end
% % idx1 = find(iptg == -3);
% % figure; plot(iptg(1:idx1(1)),yfp_ptac(1:idx1(1))); hold on;
% % plot(iptg(idx1(2):end),yfp_ptac(idx1(2):end),'g');
% % xlabel('Constitutive mgrB transcription rate \muMs^{-1}')
% % ylabel('Steady state YFP (PhoP~P activity)')
% % title('Dose response YFP vs IPTG at 0.01 mM Mg^{++}')
% % % end dose response 
%% Partial adaptation quantification
% % clear all
% % load('02112019_pso.mat'); x = Solution(1,:); % corrected detailed balance, f =/= 1
% % % load('0214_k2.mat'); x = Solution(1,:); % corrected detailed balance f = 1
% % % load('0214_KMT.mat'); x = Solution(1,:); % corrected detailed balance f = 1
% % % load('0214_KMP.mat'); x = Solution(1,:); % corrected detailed balance f = 1
% % % load('0214_k1k2combo.mat'); x = Solution(1,:); % corrected detailed balance f = 1
% % x([27 28 29]) = [0 0 1];
% % 
% % X0 = 0.01*ones(1,19);
% % % 50 mM --> 0.01 mM
% % mg = logspace(2,-2,20);
% % for i = 1:length(mg)
% %     for j = 1:length(mg)
% %         if i>j
% %     x(29) = mg(i);
% % [~, X1] =ode15s(@phopq_0619_t, [0 8*3600], X0,{},x);
% %     x(29) = mg(j);
% % [tx, X] =ode15s(@phopq_0619_t, [0 8*3600], X1(end,:),{},x);
% % % YFP:CFP partial adaptation
% % rat(i,j) = max(X(:,11)./X(:,19))/(X(end,11)./X(end,19));
% % % mRNA PmgrB partial adaptation
% % ratm(i,j) = max(X(:,4))/(X(end,4));
% % figure(7); subplot(2,1,1); plot(tx/3600, X(:,4)); hold on; title('PmgrB mRNA')
% % figure(7); subplot(2,1,2); plot(tx/3600, X(:,11)./X(:,19)); hold on; title('YFP/CFP')
% %         else
% %             rat(i,j) = 1;
% %             ratm(i,j) = 1;
% %         end
% % % pause
% %     end
% % end
% % 
% % % for i = 1:20; for j = 1:20; if(i>j); rat(i,j) = 1; end;end;end;
% % % for i = 1:20; for j = 1:20; if(i>j); ratm(i,j) = 1; end;end;end;
% % 
% % figure; 
% % subplot(1,2,1); surface(mg,mg,rat); set(gca,'XScale','log','YScale','log');
% % xlabel('[Mg^{2+}] - Post'); ylabel('[Mg^{2+}] - Pre'); title('YFP:CFP Partial adaptation')
% % 
% % subplot(1,2,2); surface(mg,mg,ratm); set(gca,'XScale','log','YScale','log');
% % xlabel('[Mg^{2+}] - Post'); ylabel('[Mg^{2+}] - Pre'); title('PmgrB mRNA Partial adaptation')

%% response to non step stimuli
% % % clear all
% % % %% parameters
% % % global fun
% % % load('0211_k1k2combo.mat'); x = Solution(2,:); k=0; x([27 28 29]) = [0 0 1]; f = 0.1*10^x(26); g = 10^x(1); fvc = ones(1,12); fvc([1 3 5]) = [g f 1/f]; x(30:41) = fvc; % corrected detailed balance; k1b = f*k1; kb1 = kb2/f; k2b = g*k2;
% % % x(42) = 0; % phoPQ not inducible
% % % y = x; y(27) = 1; z2 = x; z2(7) = 0;
% % % X0 = zeros(1,19); X0([3 4]) = 0.4*10^x(6)/10^x(2); %initialize to ~ basal levels of mgrB, phoP/Q mRNA
% % % X0del = X0; X0del(4) = 0; % initial condition 0 for mgrB deletion
% % % %% 1. WT
% % % % from Salazar et al: Cells were grown overnight in 2mM.
% % % % diluted and shifted to 50mM for 2-3h
% % % % washed and shifted to 2 or 0.01mM and imaged for 4h
% % % % 50 mM --> 0.01 mM
% % % fun = @(t) 50;
% % % [~, X] =ode15s(@phopq_0619_t, [0 20*3600], X0,{},x);
% % % [~, Xdelmgrb] =ode15s(@phopq_0619_t, [0 20*3600], X0del,{},y);
% % % [~, Xdelauto] =ode15s(@phopq_0619_t, [0 20*3600], X0del,{},z2);
% % % fun = @(t) 0.0001*(t<3600)+50*(t>3600);%0.01+0.0095*sin(t/3600);% 0.01;
% % % [t, Y]=ode15s(@phopq_0619_t, [0 20]*3600, X(end,:),{},x);
% % % [tdelmgrb, Ydelmgrb]=ode15s(@phopq_0619_t, [0 20]*3600, Xdelmgrb(end,:),{},y);
% % % [tdelauto, Ydelauto]=ode15s(@phopq_0619_t, [0 20]*3600, Xdelauto(end,:),{},z2);
% % % figure(1); hold on; plot(t/60, Y(:,2)); plot(tdelmgrb/60, Ydelmgrb(:,2)); plot(tdelauto/60, Ydelauto(:,2))

%% open loop gains
% % clear all
% % load('0211_k1k2combo.mat'); x = Solution(2,:); x([27 28 29]) = [0 0 1]; f = 0.1*10^x(26); g = 10^x(1); fvc = ones(1,12); fvc([1 3 6]) = [g f f]; x(30:41) = fvc; % corrected detailed balance; k1b = f*k1; kd1 = kd2*f; k2b = g*k2;
% % X0 = zeros(1,19); X0([3 4]) = 0.4*10^x(6)/10^x(2); %initialize to ~ basal levels of mgrB, phoP/Q mRNA
% % mgrange = 10.^linspace(2,-4,80);
% % for j = 1:length(mgrange)
% % 
% % x(28) = 0; x(42) = 0; x(29) = mgrange(j); % wild type at 1mM for normalization
% % [~, Xwt] =ode15s(@phopq_0619_t, [0 20*3600], X0,{},x);
% % f1 = 10^x(7); K1 = 10^x(8); K2 =10^x(20); f2 = 10^x(21);
% % % simulate over range of induction rates
% % h = 0.001;indrange = 0.98:h:1.02;
% % pp = zeros(length(indrange),2);
% % pp2 = zeros(length(indrange),2);
% % for i = 1:length(indrange)
% %     x(42) = 1; x(28) = 0;
% %     a = Xwt(end,2)*indrange(i); 
% %     wtind = (1+f1*(a/K1)^2)/(1+(a/K1)^2);
% %     x(5) = log10(wtind*10^x(6)*0.47);
% %     [~, Xdr] = ode15s(@phopq_0619_t, [0 20*3600], X0,{},x);
% %     pp(i,:) = [Xdr(end,2) sum(Xdr(end,[1 2 8 9 16 18]),2)];
% %     
% %     x(42) = 0; x(28) = 1; 
% %     wtind2 = (1+f2*(a/K2)^2)/(1+(a/K2)^2);
% %     x(9) = log10(wtind2*10^x(6));
% %     [~, Xdr2] = ode15s(@phopq_0619_t, [0 20*3600], X0,{},x);
% %     pp2(i,:) = [Xdr2(end,2) sum(Xdr2(end,[1 2 8 9 16 18]),2)];
% %     
% %     x(42) = 1; x(28) = 1;
% %     x(9) = log10(wtind2*10^x(6)); x(5) = log10(wtind*10^x(6)*0.47);
% %     [~, Xdr3] = ode15s(@phopq_0619_t, [0 20*3600], X0,{},x);
% %     pp3(i,:) = [Xdr3(end,2) sum(Xdr3(end,[1 2 8 9 16 18]),2)];
% % end
% % OLG(j) = mean(diff(pp(:,1))/(h*Xwt(end,2)));
% % OLG2(j) = mean(diff(pp2(:,1))/(h*Xwt(end,2)));
% % OLG3(j) = mean(diff(pp3(:,1))/(h*Xwt(end,2)));
% % end
% % 
% % figure; semilogx(mgrange, OLG, mgrange,OLG2,mgrange,OLG3)
% % title('partial and complete open loop gains'); ylabel('\partial[PhoP~P]/\partialR_0'); xlabel('Signal k_{-1}/k_{-1}^0'); legend('PF open','NF open','both open')

%% Impact of each free parameter between SK and SKX
% % load('12082019_all_0619t.mat'); 
% % J = [1 26 27:32; 12 10 14:19]';
% % for j = 1:5
% % for k = 1:8
% %     x = Solution(j,:); 
% %     x(J(k,1)) = x(J(k,2));
% %     err(j,k) = phopq_error_0619_t(x);
% % end
% % end
% % figure; bar((err(:,1:end-1)./err(:,end))');
% % ylim([0 200]); set(gca,'XTickLabel',{'k_1','k_2','k_3','k_{-3}','k_{4}','k_5','k_{-5}','k_6'})
% % ylabel('Error fold change')
