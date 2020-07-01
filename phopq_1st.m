function dY = phopq_1st(t, Y, x)
% Generic version without explicit michaelis menten
% single PhoQ form that both phosphorylates and dephosphorylates
% dimensionless factors/inputs
fvc = ones(1,8); fvc(1) = 10^x(1); 
% % % comment out for anything but 1217_all parameters
% fvc(3:8) = 10.^x(31:36)./10.^x(14:19);
% % % comment out for anything but 1217_all parameters
mutE = x(27);
mgrBconst = x(28);
mg = x(29); 
phopqind = x(30);  % inducible phoPQ %
% % % % % % k2 increase with stimulus
k20 = 10^x(12); F = 10^x(10); Km = 10^x(11); n = 10^x(22);
k2 = k20*F/(1+(mg/Km)^n); 
    k5 = 10^x(17);
% % % % % % k2 decrease with stimulus, k5 decrease with stimulus
% k20 = 10^x(12); F = 10^x(10); Km = 10^x(11); n = 10^x(22);
% k2 = k20*(1+F*(mg/Km).^n)./(1+(mg/Km).^n);
% k50 = 10^x(17); Km2 = 10^x(23); n2= 10^x(26); k5 = k50*(mg/Km2).^n2./(1+(mg/Km2).^n2);
% % % % % % %

kpdeg_eff = 3.1e-4;
kmdeg = 10^x(2);
ktlnA = 10^x(3);
ktlnE = 10^x(4);
kbtpn1 = 0.47*10^x(6); 
kbtpn2 = 10^x(6);
f1 = 10^x(7);
K1 = 10^x(8);
ktlnB = ktlnA/40;
km2 = 10^x(13);
k3= 10^x(14);
km3= 10^x(15);
k4= 10^x(16);

km5 = 10^x(18);
k6 = 10^x(19);
K2 =10^x(20);
f2 = 10^x(21);
ktlnY = 10^x(23);
kb1 = 10^x(24);
kd1 = 10^x(25);
k2b = k2*fvc(1);
km2b = km2*fvc(2);
k3b = k3*fvc(3);
km3b = km3*fvc(4);
k4b = k4*fvc(5); 
k5b = k5*fvc(6);
km5b =km5*fvc(7);
k6b = k6*fvc(8);

P=Y(1);     %PhoP
Pp=Y(2);    %PhoP~P
Q=Y(3);     %PhoQ - ph 
Qp=Y(4);    %PhoQ~P -kin
MgrB = Y(5); %MgrB
QB = Y(6);   %PhoQ.MgrB -ph
QBp = Y(7); %PhoQ~P.MgrB -kin
QsP = Y(8); % PhoQ~P.PhoP
QP = Y(9);  % PhoQ.PhoP~P
QsBP = Y(10); %PhoQ~P.MgrB.PhoP
QBP = Y(11); %PhoQ.MgrB.PhoP~P
mPQ = Y(12);
mB = Y(13);
mYFP = Y(14);
yfp = Y(15);
yfp2 = Y(16);
cfp = Y(17);
%---------------
    dY(1) = (ktlnA*mPQ-k3*Qp*P + km3*QsP +k6*QP - k3b*QBp*P + km3b*QsBP +k6b*QBP- kpdeg_eff*P); % PhoP
    dY(2) = (k4*QsP - k5*Q*Pp + km5*QP- k5b*QB*Pp + km5b*QBP +k4b*QsBP -kpdeg_eff*Pp); % PhoP~P
    dY(3) = ktlnB*mPQ -kb1*Q*MgrB + kd1*QB - k5*Q*Pp + km5*QP + k6*QP - k2*Q + km2*Qp + k4*QsP - kpdeg_eff*Q; % PhoQ
    dY(4) = + k2*Q - km2*Qp + km3*QsP - k3*Qp*P-kpdeg_eff*Qp; % PhoQ~P
    dY(5) = ktlnE*mB -kb1*Q*MgrB + kd1*QB - kpdeg_eff*MgrB; %MgrB
    dY(6) = kb1*Q*MgrB - kd1*QB - k5b*QB*Pp + km5b*QBP + k6b*QBP - k2b*QB + km2b*QBp + k4b*QsBP - kpdeg_eff*QB; %QB [PhoQ.MgrB]
    dY(7) = + k2b*QB - km2b*QBp - k3b*QBp*P + km3b*QsBP - kpdeg_eff*QBp; %[PhoQ~P.MgrB]
    dY(8) = k3*Qp*P - km3*QsP - k4*QsP-kpdeg_eff*QsP; % [PhoQ~P.PhoP]
    dY(9) = k5*Q*Pp - km5*QP - k6*QP-kpdeg_eff*QP; % [PhoQ.PhoP~P]
    dY(10) = k3b*QBp*P - km3b*QsBP - k4b*QsBP-kpdeg_eff*QsBP; %QsBP [PhoQ~P.MgrB.PhoP]
    dY(11) = k5b*QB*Pp - km5b*QBP - k6b*QBP-kpdeg_eff*QBP; %QBP - [PhoQ.MgrB.PhoP~P]
dY(12)= phopqind*10^x(5)+(1-phopqind)*(kbtpn1*(1+f1*(Pp/K1)^2)/(1+(Pp/K1)^2))-kmdeg*mPQ;%mphoPphoQ % autoregulated
dY(13)=  mgrBconst*10^x(9) + (1-mgrBconst)*(1-mutE)*(kbtpn2*((1+f2*(Pp/K2)^2)/(1+(Pp/K2)^2))) - kmdeg*mB; %mgrB
dY(14) = kbtpn2*((1+f2*(Pp/K2)^2)/(1+(Pp/K2)^2)) - kmdeg*mYFP; %mRNA YFP from PmgrB
dY(15) = ktlnY*mYFP - kpdeg_eff*yfp; % YFP from P-mgrB
dY(16) = ktlnY*mPQ - kpdeg_eff*yfp2; % YFP from Pphopq
dY(17) = ktlnY*kbtpn1/kmdeg - kpdeg_eff*cfp;
dY=dY';
end