%% steady state solutions to PhoPQ-MgrB post-translational reactions coupled with transcriptional negative feedback of MgrB.
clear all
% general parameters
% k2 = 0.01; 
% km2 = 0.001 ;
% k3= 1;
% km3= 0.001;
% k4= 1;
% k5 = 1;
% km5 = 0.001;
% k6 = 1;
% kb2 = 1;
% kd2 = 0.001; KD = kd2/kb2;
% KMP = (km5+k6)/k5; KMT = (km3+k4)/k3; Ct = km2*KMT/k4; Cp = k2*KMP/k6;
% l = 0.1;
% K2 = 5*(Ct+Cp);
% BTmin = 5*(Ct+Cp); xt0 = BTmin; % set basal X-total
% f2 = 60;
% par = [k2 km2 k3 km3 k4 k5 km5 k6 kb2 kd2 l f2 K2];

% parameters from equivalent two state model
k2 = 2; 
km2 = 0.0001 ;
k3= 5;
km3= 0.01;
k4= 4;
k5 = 0.1;
km5 = 0.0007;
k6 = 1.5;
kb2 = .4;
kd2 = 0.0008; KD = kd2/kb2;
KMP = (km5+k6)/k5; KMT = (km3+k4)/k3; Ct = km2*KMT/k4; Cp = k2*KMP/k6;
l = 2.7e-4;
K2 = 0.2;
BTmin = .14;
xt0 = BTmin; % set basal X-total
f2 = 50;
par = [k2 km2 k3 km3 k4 k5 km5 k6 kb2 kd2 l f2 K2];

% solution process
s0 = zeros(1,11); % initial condition
rrt = logspace(-3,3,30); % range of RR-total values over which steady state evaluated

% %         % This indented section for computing RR~P as a function of
% %         % independent X-total. For graph of RR_P vs XT module
% %         xtr = linspace(1,f2,40)*BTmin;
% %         for m = 1:length(xtr)
% %             s0 = zeros(1,11); s0(1) = 0.1; s0(3) = s0(1)/50; s0(5) = xtr(m);
% %             [~, X] =ode15s(@phopq_1st_PPI, [0 100*3600], s0,{},par); 
% %             rrp_module(m) = X(end,2);
% %         end
% %         figure; 
% %         loglog(rrp_module/(Cp+Ct),xtr/(Cp+Ct),'r','linewidth',1.5); hold on;
% %         rrp_module2 = logspace(-4,2,40)*K2; loglog(rrp_module2/(Cp+Ct), (BTmin*(1+f2*(rrp_module2/K2).^2)./(1+(rrp_module2/K2).^2))/(Cp+Ct),'k','linewidth',1.5)
% %         % steady state analytical solution at a given RR-total
% %         fun = @(y0) indc(y0,0,BTmin,par); anxt = fmincon(fun,BTmin); anrrp = Cp*(1+l*anxt/KD)/(1+anxt/KD); loglog(anrrp*ones(1,20)/(Cp+Ct),linspace(1,f2,20)*BTmin/(Cp+Ct),'k--')

% steady state solution solving post-translational module with
% transcriptional feedback.
for k = 1:length(rrt)
    fun = @(y0) indc(y0,rrt(k),BTmin,par);
%     XT(k) = fmincon(fun,xt0,[],[],[],[],BTmin,f2*BTmin,[],options); 
%     XT(k) = fsolve(fun,xt0,options);
%     xt0 = XT(k);
    XT(k) = fminbnd(fun,BTmin,f2*BTmin);
    s0 = zeros(1,11); s0(1) = rrt(k); s0(3) = s0(1)/40; s0(5) = XT(k);
    [~, Xind] =ode15s(@phopq_1st_PPI, [0 100*3600], s0,{},par);    
    rrp(k,:) = Xind(end,:);
end
figure(9);
% yyaxis left
% subplot(1,2,1);
loglog(rrt,rrp(:,2),'linewidth',2,'linestyle','-'); hold on %rrp vs RRT 
ylabel('[PhoP~P]'); xlabel('[PhoP]_{Total}')
% yyaxis right
% subplot(1,2,2);
% loglog(rrt,rrp(:,5)./(sum(rrp(:,[5 6 7 10 11]),2)),'color',[1 0.5 0],'linewidth',2) % X/X-total vs RRT
% ylabel('[B]/[B]_T')
% loglog(rrt, XT./(rrt/50))

function ERR = indc(y0,RRT,BTmin, par)
    X0 = zeros(1,11);
    X0([1 3 5]) = [RRT RRT/40 y0];
    [~, Xind] =ode15s(@phopq_1st_PPI, [0 100*3600], X0,{},par);    
    Pp = Xind(end,2);
    
    km2 =par(2);k3= par(3);km3= par(4);k4= par(5); KMT = (km3+k4)/k3; Ct = km2*KMT/k4;
    k2 = par(1); k5 = par(6); km5 = par(7); k6 = par(8); KMP = (km5+k6)/k5; Cp = k2*KMP/k6;
    f2 = par(12); K2 = par(13);
    
%     KD = par(10)/par(9); l = par(11);
%     Pp = Cp*(1+l*y0/KD)/(1+y0/KD);
    ERR = (y0-BTmin*(1+f2*(Pp/K2).^2)./(1+(Pp/K2).^2)).^2;
end