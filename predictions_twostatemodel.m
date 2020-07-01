clear all
%% parameters
% combos:
    % 1. k1 k2
% load('0214_k1k2combo.mat'); x = Solution(1,:); k=0; x([27 28 29]) = [0 0 1]; f = 10^x(26); g = 0*10^x(1); fvc = ones(1,12); fvc([1 3 6]) = [g f f]; x(30:41) = fvc; % corrected detailed balance; k1b = f*k1; kd1 = kd2*f; k2b = g*k2;
load('0211_k1k2combo.mat'); x = Solution(2,:); k=0; x([27 28 29]) = [0 0 1]; f = 0.1*10^x(26); g = 10^x(1); fvc = ones(1,12); fvc([1 3 6]) = [g f f]; x(30:41) = fvc; % corrected detailed balance; k1b = f*k1; kd1 = kd2*f; k2b = g*k2;

x(42) = 0; % phoPQ not inducible

% define all mutant parameters:
y = x; y(27) = 1; y(9) = -Inf; % -mgrB or mgrB constitutive
z = x; z([27 7]) = [1 0]; z(9) = -Inf; %-mgrB;-autoreg
z2 = x; z2(7) = 0; %-autoreg
z3 = x; z3(14:19) = x(14:19)+k; z4 = z3; z4(27) = 1; % z3:-phosphatase; %z4 = -phosphatase;-mgrB
% z3(7) = 0; z4(7) = 0; %deltaautoregphoqt281r

% initial vectors for +/- mgrB
X0 = zeros(1,19); X0del = X0;
%% dose -response predictions:
% % x(21)=log10(25); %f2
% % x(20) = -0.5; %K2
% % x(4)=-1.935; %ktlnE
y(9) = x(6);
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
figure;
semilogx(mgrange, pp(:,1)/pp(idx,1),'b');
hold on; semilogx(mgrange, pp(:,2)/pp(idx,1),'g');
semilogx(mgrange, pp(:,3)/pp(idx,1),'k');
phopq_dr(:,1) = [0.03 0.1 0.3 1 3 10 30]; % Mg (mM)
phopq_dr(:,2) = [0.58 0.58 0.56 0.5 0.39 0.19 0.1]; % 
semilogx(phopq_dr(:,1),phopq_dr(:,2)/phopq_dr(end,2),'s')
xlabel('Signal, k_{-1} (representing [Mg^{2+}])'); ylabel('Normalized YFP:CFP');
xlim([mgrange(1) mgrange(end)])
legend('WT','\DeltamgrB','\Deltaautoreg','WT-expt')
title('Dose-response');
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
% normalize to wild type phoPQ induction
x(29) = 1; [~, X_norm] =ode15s(@phopq_0619_t, [0 40*3600], X0,{},x);
ppnorm = [X_norm(end,11)/X_norm(end,19) sum(X_norm(end,[1 2 8 9 16 18]),2)];

figure(4); plot(pp_const(:,2)/ppnorm(2),pp_const(:,1)/ppnorm(1),'m'); hold on;
xlabel('{\it phoPphoQ} induction ([PhoP]/[PhoP]_{WT})')
ylabel('Normalized P_{mgrB} output')
legend('WT simulation','Expt','constitutive mgrB')

%% oscillation predictions with constitutive mgrB
y(9) = x(6);
% % y(7) = 0; % autoreg off
y(29) = 2;
[~, X1delmgr] =ode15s(@phopq_0619_t, [0 10*3600], X0del,{},y);
y(29) = 50;
[~, X_delmgr] =ode15s(@phopq_0619_t, [0 10*3600], X1delmgr(end,:),{},y);
y(29) = 1;
[t_delmgr, Y_delmgr]=ode15s(@phopq_0619_t, [0 20]*3600, X_delmgr(end,:),{},y);
figure(26);
plot(t_delmgr/60, (Y_delmgr(:,11)./Y_delmgr(:,19)),'m-'); hold on;%Delta mgrB
xlabel('time (min)'); ylabel('YFP:CFP'); title('50 \rightarrow 1 mM; constitutive mgrB')

figure; plot(t_delmgr/60,Y_delmgr(:,[3])) % mPQ
hold on;
plot(t_delmgr/60,Y_delmgr(:,[2])) %PhoP-P
plot(t_delmgr/60,sum(Y_delmgr(:,[5 7 8]),2)) % Q kin
plot(t_delmgr/60,sum(Y_delmgr(:,[6 9]),2)) % Q ph
plot(t_delmgr/60,sum(Y_delmgr(:,[14 15 16]),2)) % QB kin
plot(t_delmgr/60,sum(Y_delmgr(:,[17 18]),2)) % QB - ph

legend('m_{phoPQ}','PhoP-P','Q_{kin}','Q_{ph}','QB_{kin}','QB_{ph}')
set(gca,'YScale','log','linewidth',1)
xlabel('time (mins)'); xlim([-10 400])
%% open loop gain of oscillating system
y(29) = 10; y(9) = x(6);
[t0, Y0]=ode15s(@phopq_0619_t, [0 100]*3600, X0,{},y);
y(29) = 0.2; 
[t_delmgr_ss, Y_delmgr_ss]=ode15s(@phopq_0619_t, [0 100]*3600, Y0(end,:),{},y);

x(29) = 10;
[~, Y0wt]=ode15s(@phopq_0619_t, [0 100]*3600, X0,{},x);
x(29) = 0.2; 
[t_delmgr_wt, Y_delmgr_wt]=ode15s(@phopq_0619_t, [0 100]*3600, Y0wt(end,:),{},x);



figure; plot(t_delmgr_ss/3600, Y_delmgr_ss(:,2))
plot(t_delmgr_ss/3600, sum(Y_delmgr_ss(:,[10 14:18]),2)); hold on;
plot(t_delmgr_ss/3600, sum(Y_delmgr_ss(:,[5:9 14:18]),2))
plot(t_delmgr_ss/3600, sum(Y_delmgr_ss(:,[5:9]),2))

figure;
plot(t_delmgr_ss/3600,sum(Y_delmgr_ss(:,[5 7 8]),2)); hold on; % Q kin
plot(t_delmgr_ss/3600,sum(Y_delmgr_ss(:,[6 9]),2)) % Q ph
plot(t_delmgr_ss/3600,sum(Y_delmgr_ss(:,[14 15 16]),2)) % QB kin
plot(t_delmgr_ss/3600,sum(Y_delmgr_ss(:,[17 18]),2)) % QB - ph
legend('Q_{kin}','Q_{ph}','QB_{kin}','QB_{ph}')

pprange = linspace(min(Y_delmgr_ss(400:end,2)),max(Y_delmgr_ss(400:end,2)),100);

for i = 1:length(pprange)
delmgrpp= pprange(i); %Y_delmgr_ss(end,2);
kbtpn1 = 0.47*10^x(6);f1 = 10^x(7);K1 = 10^x(8);
y(42) = 1; y(5) = log10(kbtpn1*(1+f1*(delmgrpp/K1)^2)/(1+(delmgrpp/K1)^2));
[~, Y_olg]=ode15s(@phopq_0619_t, [0 100]*3600, X0,{},y);
ppout(i) = Y_olg(end,2);
end
diffpp = (ppout(2:end)-ppout(1:end-1))/(pprange(2)-pprange(1))

% delmgrpp= delmgrpp*1.01; y(5) = log10(kbtpn1*(1+f1*(delmgrpp/K1)^2)/(1+(delmgrpp/K1)^2));
% [~, Y_olg_plus]=ode15s(@phopq_0619_t, [0 100]*3600, X0,{},y);
% delmgrpp= delmgrpp*0.99; y(5) = log10(kbtpn1*(1+f1*(delmgrpp/K1)^2)/(1+(delmgrpp/K1)^2));
% [~, Y_olg_minus]=ode15s(@phopq_0619_t, [0 100]*3600, X0,{},y);
[Y_olg(end,2) Y_olg_plus(end,2) Y_olg_minus(end,2)]
((Y_olg_plus(end,2)-Y_olg_minus(end,2))/(0.2*Y_delmgr_ss(end,2)))*(Y_delmgr_ss(end,2)/Y_olg(end,2))


%%

y(9) = x(6);
% mgrange = 10.^(-4.5:0.05:1.7); % dosage range of mg
mgrange = logspace(0.31,-1,10);
pp = zeros(length(mgrange),3);
y(29) = 10; % mg = 1mM
[~, X2] =ode15s(@phopq_0619_t, [0 20*3600], X0del,{},y); 
trange = [0 50*3600];
for i = 1:length(mgrange)
    y(29) = mgrange(i); 
    [ton, Xdrdelmgr] =ode15s(@phopq_0619_t, trange, X2(end,:),{},y);
figure(22); plot(ton/3600,Xdrdelmgr(:,2)); hold on;
pause
    rrp(i) = Xdrdelmgr(end,2);
end