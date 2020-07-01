function dY = phopq_0619_t(t, Y, x)
% General comprehensive network - changes from 0605 - 
% PhoQ binds MgrB in any form
fvc = x(30:41);
mutE = x(27); % mgrB deletion marker
mg = x(29);
phopqind = x(42);  % inducible phoPQ %
    % nonlinear Mg %
    %     km = 10^x(9);
    %     n = 10^x(5);
    %     km10 = 0;
    %     km1max = 10^x(11);
    %     km1 = km10 + km1max*(mg^n/(mg^n + km^n)); %alternative, nonlinear mg dependence
    % %                 mg = logspace(-3,2,20); km1 = km10 + km1max*(mg.^n./(mg.^n + km^n)); figure; loglog(mg, km1)
    % nonlinear Mg %
kpdeg_eff = 3.1e-4;
kmdeg =10^x(2);
ktlnA =10^x(3); 
ktlnE = 10^x(4); 
kbtpn1 = 0.47*10^x(6); 
kbtpn2 = 10^x(6);
f1 = 10^x(7);
K1 = 10^x(8);
ktlnB = ktlnA/40;
k1 = 10^x(10);
km1 = 10^x(11)*mg;
k2 = 10^x(12);
km2 = 10^x(13);
k3= 10^x(14);
km3= 10^x(15);
k4= 10^x(16);
k5 = 10^x(17);
km5 = 10^x(18);
k6 = 10^x(19);
K2 =10^x(20);
f2 = 10^x(21);
ktlnY = 10^-1.3;
kb2 = 10^x(24); % binding SK*, X
kd2 = 10^x(25); % SKX* diss
kb1 = kb2; % binding SK, X
kd1 = kd2*fvc(3); % SKX diss 
% assume kd1/kd2 = k1b/k1, binding rates equal
k2b = k2*fvc(1);
km2b = km2*fvc(2);
k1b = k1*fvc(3);
km1b =km1;
k3b = k3*fvc(7);
km3b = km3*fvc(8);
k4b = k4*fvc(9); 
k5b = k5*fvc(10);
km5b =km5*fvc(11);
k6b = k6*fvc(12);

% % % Important note: Thermodynamic constraint % % %
% (k1b/km1b)/(k1/km1) = 1/((kb1/kd1)/(kb2/kd2))
% wild type markers

P=Y(1);     %PhoP
Pp=Y(2);    %PhoP~P
mPQ=Y(3);   %mRNA phoPQ
mB=Y(4);    %mRNA mgrB
Qs=Y(5);    %PhoQ* -kin
Q=Y(6);     %PhoQ - ph 
Qp=Y(7);    %PhoQ~P -kin
QsP = Y(8); %PhoQ~P.PhoP -kin
QP = Y(9);  %PhoQ.PhoP~P -ph
MgrB = Y(10); %MgrB
yfp = Y(11);  %YFP from PmgrB
yfp2 = Y(12); %YFP from Pphopq
mYFP = Y(13); %YFP mRNA from PmgrB
QsB = Y(14);  %PhoQ*.MgrB -kin
QBp = Y(15);  %PhoQ~P.MgrB -kin
QsBP = Y(16); %PhoQ~P.PhoP.MgrB -kin
QB = Y(17);   %PhoQ.MgrB -ph
QBP = Y(18);  %PhoQ.PhoP~P.MgrB -ph
cfp = Y(19);  %CFP

% kpdeg_eff = 1.5*(mg<=0.01)*((t<240*60)*(2.0576 - 5.673e-5*t)*1e-4 + (t>240*60)*1.24e-4) + 1.5*(mg>0.01)*2.0576e-4;
% kpdeg_eff = 3.1e-5 + (3.1e-4-3.1e-5)*mg./(0.08 + mg);
%---------------
dY(1) = (ktlnA*mPQ-k3*Qp*P + km3*QsP +k6*QP - k3b*QBp*P + km3b*QsBP +k6b*QBP- kpdeg_eff*P); % PhoP
dY(2) = (k4*QsP - k5*Q*Pp + km5*QP- k5b*QB*Pp + km5b*QBP +k4b*QsBP -kpdeg_eff*Pp); % PhoP~P
dY(6) = + km1*Qs - k1*Q - kb1*Q*MgrB+ kd1*QB- k5*Q*Pp + km5*QP + k6*QP - kpdeg_eff*Q; % PhoQ
dY(5) = ktlnB*mPQ + k1*Q - km1*Qs - k2*Qs + km2*Qp + k4*QsP -kb2*Qs*MgrB + kd2*QsB - kpdeg_eff*Qs; % PhoQ*
% New synthesized PhoQ - kinase
dY(7) = + k2*Qs - km2*Qp + km3*QsP - k3*Qp*P-kpdeg_eff*Qp; % PhoQ*~P
dY(10) = ktlnE*mB -kb2*Qs*MgrB+kd2*QsB -kb1*Q*MgrB + kd1*QB - kpdeg_eff*MgrB; %MgrB
dY(17) = kb1*Q*MgrB - kd1*QB - k1b*QB + km1b*QsB - k5b*QB*Pp + km5b*QBP + k6b*QBP - kpdeg_eff*QB; %QB [PhoQ.MgrB]
dY(14) = + kb2*Qs*MgrB - kd2*QsB +k1b*QB - km1b*QsB - k2b*QsB + km2b*QBp + k4b*QsBP - kpdeg_eff*QsB; %QB* - [PhoQ*.MgrB]
dY(15) = + k2b*QsB - km2b*QBp - k3b*QBp*P + km3b*QsBP - kpdeg_eff*QBp; %QB*~P - [PhoQ*~P.MgrB]
dY(8) = k3*Qp*P - km3*QsP - k4*QsP-kpdeg_eff*QsP; % [PhoQ*~P.PhoP]
dY(9) = k5*Q*Pp - km5*QP - k6*QP-kpdeg_eff*QP; % [PhoQ.PhoP~P]
dY(16) = k3b*QBp*P - km3b*QsBP - k4b*QsBP-kpdeg_eff*QsBP; %QsBP [PhoQ*~P.MgrB.PhoP]
dY(18) = k5b*QB*Pp - km5b*QBP - k6b*QBP-kpdeg_eff*QBP; %QBP - [PhoQ.MgrB.PhoP~P]
dY(3)= phopqind*10^x(5)+(1-phopqind)*(kbtpn1*(1+f1*(Pp/K1)^2)/(1+(Pp/K1)^2))-kmdeg*mPQ;%mphoPphoQ % autoregulated
dY(4)=  mutE*10^x(9) + (1-mutE)*(kbtpn2*((1+f2*(Pp/K2)^2)/(1+(Pp/K2)^2))) - kmdeg*mB; %mgrB
dY(13) = kbtpn2*((1+f2*(Pp/K2)^2)/(1+(Pp/K2)^2)) - kmdeg*mYFP; %mRNA YFP from PmgrB
dY(11) = ktlnY*mYFP - kpdeg_eff*yfp; % YFP from P-mgrB
dY(12) = ktlnY*mPQ - kpdeg_eff*yfp2; % YFP from Pphopq
dY(19) = ktlnY*kbtpn1/kmdeg - kpdeg_eff*cfp;
dY=dY';
end