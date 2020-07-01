y0 = zeros(1,9);
load('0211_k1k2combo.mat'); x = Solution(1,:); k=0; x([27 28 29]) = [0 0 1]; f = 10^x(26); g = 10^x(1); fvc = ones(1,12); fvc([1 3 6]) = [g f f]; x(30:41) = fvc; % corrected detailed balance; k1b = f*k1; kd1 = kd2*f; k2b = g*k2;
x(42) = 0; % phoPQ not inducible
y = x; y(27) = 1; y(9) = x(6);
y(29) = 50; [~,Y50] = ode15s(@PQBshort2,[0 20]*3600, y0,{},y);
y(29) = 1; [t,Y] = ode15s(@PQBshort2,[0 20]*3600, Y50(end,:),{},y);
figure(1); plot(t/60, Y(:,[2,3,4,7,9]));

% function dY = PQBshort2(~,Y,x)
% % Shorter cousin of @PQB assuming all 4 catalytic reaction quasi steady
% % state. Tracks total SK-kin, SK-ph, SKX-kin and SKX-ph. Substrate bound
% % states calculated within function. Other than these 4 variables - RR, RRP
% % and X tracked.
% % global variable prd can be set to 1 to remove synthesis or degradation.
% % global prd
% prd = 1;
% kpdeg_eff = 3.1e-4;
% 
% % x(13) = -Inf;
% fvc = x(30:41);
% mg = x(29);
% k1 = 10^x(10);
% km1 = 10^x(11)*mg;
% k2 = 10^x(12);
% km2 = 10^x(13);
% k3= 10^x(14);
% km3= 10^x(15);
% k4= 10^x(16);
% k5 = 10^x(17);
% km5 = 10^x(18);
% k6 = 10^x(19);
% kb2 = 10^x(24);
% kd2 = 10^x(25);
% k2b = k2*fvc(1);
% km2b = km2*fvc(2);
% k1b = k1*fvc(3);
% km1b =km1*fvc(4);
% kb1 = kb2*fvc(5);
% kd1 = kd2*fvc(6);
% k3b = k3*fvc(7);
% km3b = km3*fvc(8);
% k4b = k4*fvc(9); 
% k5b = k5*fvc(10);
% km5b =km5*fvc(11);
% k6b = k6*fvc(12);
% kmdeg =10^x(2);
% ktlnA =10^x(3); %*10^x(2);
% ktlnE = 10^x(4); %*10^x(2);
% kbtpn2 = 10^x(6);
% kbtpn1 = 0.47*10^x(6);
% f1 = 10^x(7);
% K1 = 10^x(8);
% K2 =10^x(20);
% f2 = 10^x(21);
% KMT = (km3+k4+kpdeg_eff)/k3;KMTb = (km3b+k4b+kpdeg_eff)/k3b;
% KMP = (km5+k6+kpdeg_eff)/k5;KMPb = (km5b+k6b+kpdeg_eff)/k5b;
% % % % Important note: Thermodynamic constraint % % %
% % (k1b/km1b)/(k1/km1) = 1/((kb1/kd1)/(kb2/kd2))
% % wild type markers
% 
% RR=Y(1);     %PhoP
% RRP=Y(2);    %PhoP~P
% SKkin=Y(3);    %PhoQ* -kin
% SKph=Y(4);     %PhoQ - ph 
% SKP = Y(5);    %PhoQ~P -kin
% % X = Y(6); %MgrB
% SKXkin = Y(7);  %PhoQ*.MgrB -kin
% SKXP = Y(8);  %PhoQ~P.MgrB -kin
% SKXph = Y(9);   %PhoQ.MgrB -ph
% 
% SKXRRP = SKXph*(RRP/KMPb)/(1+RRP/KMPb);
% SKXPRR = 0;%SKXkin*(RR/KMTb)/(1+RR/KMTb);
% SKRRP = SKph*(RRP/KMP)/(1+RRP/KMP);
% SKPRR = (SKkin)*(RR/KMT)/(1+(1+k4/k2)*RR/KMT);
% SKs = SKPRR*k4/k2;
% SK = SKph-SKRRP;
% SKXs = SKXkin;
% SKX = SKXph-SKXRRP;
% 
% 
% prP = prd*(ktlnA*(kbtpn1*((1+f1*(RRP/K1)^2)/(1+(RRP/K1)^2)))/kmdeg); prQ = prP/40;
% prB = prd*ktlnE*10^x(9)/kmdeg; %prd*(ktlnE*(kbtpn2*((1+f2*(RRP/K2)^2)/(1+(RRP/K2)^2)))/kmdeg);
% X = ktlnE*10^x(9)/(kmdeg*kpdeg_eff); 
% %---------------
% dY(1) = prP -k4*SKPRR + k6*SKRRP + k6b*SKXRRP - kpdeg_eff*RR; % PhoP
% dY(2) = (k4*SKPRR - k6*SKRRP - k6b*SKXRRP - kpdeg_eff*RRP); % PhoP~P
% dY(3) = 1*prQ + k1*SK - km1*SKs - kb2*SKs*X + kd2*SKXs - kpdeg_eff*SKkin; % PhoQ* chk
% dY(4) = + km1*SKs - k1*SK - kb1*SK*X+ kd1*SKX - kpdeg_eff*SKph; % PhoQ chk
% % dY(6) = prB -kb2*SKs*X+kd2*SKXs -kb1*SK*X + kd1*SKX - kpdeg_eff*X; %MgrB
% dY(7) = kb2*SKs*X - kd2*SKXs +k1b*SKX - km1b*SKXs - k2b*SKXs - kpdeg_eff*SKXkin; %QB* - [PhoQ*.MgrB] chk
% dY(9) = kb1*SK*X - kd1*SKX - k1b*SKX + km1b*SKXs - kpdeg_eff*SKXph; %QB [PhoQ.MgrB] chk
% dY=dY';
% end

function dY = PQBshort2(~,Y,x)
% Shorter cousin of @PQB assuming all 4 catalytic reaction quasi steady
% state. Tracks total SK-kin, SK-ph, SKX-kin and SKX-ph. Substrate bound
% states calculated within function. Other than these 4 variables - RR, RRP
% and X tracked.
% global variable prd can be set to 1 to remove synthesis or degradation.
% global prd
prd = 1;
kpdeg_eff = 3.1e-4;

% x(13) = -Inf;
fvc = x(30:41);
mg = x(29);
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
kb2 = 10^x(24);
kd2 = 10^x(25);
k2b = k2*fvc(1);
km2b = km2*fvc(2);
k1b = k1*fvc(3);
km1b =km1*fvc(4);
kb1 = kb2*fvc(5);
kd1 = kd2*fvc(6);
k3b = k3*fvc(7);
km3b = km3*fvc(8);
k4b = k4*fvc(9); 
k5b = k5*fvc(10);
km5b =km5*fvc(11);
k6b = k6*fvc(12);
kmdeg =10^x(2);
ktlnA =10^x(3); %*10^x(2);
ktlnE = 10^x(4); %*10^x(2);
kbtpn2 = 10^x(6);
kbtpn1 = 0.47*10^x(6);
f1 = 10^x(7);
K1 = 10^x(8);
K2 =10^x(20);
f2 = 10^x(21);
KMT = (km3+k4+kpdeg_eff)/k3;KMTb = (km3b+k4b+kpdeg_eff)/k3b;
KMP = (km5+k6+kpdeg_eff)/k5;KMPb = (km5b+k6b+kpdeg_eff)/k5b;



RRP=Y(2);    %PhoP~P
prP = prd*(ktlnA*(kbtpn1*((1+f1*(RRP/K1)^2)/(1+(RRP/K1)^2)))/kmdeg); prQ = prP/40;
prB = prd*ktlnE*10^x(9)/kmdeg; %prd*(ktlnE*(kbtpn2*((1+f2*(RRP/K2)^2)/(1+(RRP/K2)^2)))/kmdeg); % 
X = prB/kpdeg_eff; 
Pt = kbtpn1*ktlnA/(kmdeg*kpdeg_eff); %prP/kpdeg_eff; %
RR= Y(1); %Pt - RRP;%     %PhoP
SKkin=Y(3);    %PhoQ* -kin
SKph=Y(4);     %PhoQ - ph 
% X = Y(6); %MgrB
SKXkin = Y(7);  %PhoQ*.MgrB -kin
SKXph = Y(9);   %PhoQ.MgrB -ph

SKXRRP = SKXph*(RRP/KMPb)/(1+RRP/KMPb);
SKXPRR = 0;%SKXkin*(RR/KMTb)/(1+RR/KMTb);
SKRRP = SKph*(RRP/KMP)/(1+RRP/KMP);
SKPRR = (SKkin)*(RR/KMT)/(1+(1+k4/k2)*RR/KMT);
SKs = SKPRR*k4/k2;
SK = SKph-SKRRP;
SKXs = SKXkin;
SKX = SKXph-SKXRRP;


%---------------
dY(1) = prP -k4*SKPRR + k6*SKRRP + k6b*SKXRRP - kpdeg_eff*RR; % PhoP
dY(2) = (k4*SKPRR - k6*SKRRP - k6b*SKXRRP - kpdeg_eff*RRP); % PhoP~P
% dY(2) = (k4*SKPRR - k6*SKRRP - k6b*SKXRRP); % PhoP~P
dY(3) = 1*prQ + k1*SK - km1*SKs - kb2*SKs*X + kd2*SKXs - kpdeg_eff*SKkin; % PhoQ* chk
dY(4) = + km1*SKs - k1*SK - kb1*SK*X+ kd1*SKX - kpdeg_eff*SKph; % PhoQ chk
% dY(6) = prB -kb2*SKs*X+kd2*SKXs -kb1*SK*X + kd1*SKX - kpdeg_eff*X; %MgrB
dY(7) = kb2*SKs*X - kd2*SKXs +k1b*SKX - km1b*SKXs - k2b*SKXs - kpdeg_eff*SKXkin; %QB* - [PhoQ*.MgrB] chk
dY(9) = kb1*SK*X - kd1*SKX - k1b*SKX + km1b*SKXs - kpdeg_eff*SKXph; %QB [PhoQ.MgrB] chk
dY=dY';
end

% function dY = PQBshort2(~,Y,x)
% % Shorter cousin of @PQB assuming all 4 catalytic reaction quasi steady
% % state. Tracks total SK-kin, SK-ph, SKX-kin and SKX-ph. Substrate bound
% % states calculated within function. Other than these 4 variables - RR, RRP
% % and X tracked.
% % global variable prd can be set to 1 to remove synthesis or degradation.
% % global prd
% prd = 1;
% kpdeg_eff = 3.1e-4;
% 
% % x(13) = -Inf;
% fvc = x(30:41);
% mg = x(29);
% k1 = 10^x(10);
% km1 = 10^x(11)*mg;
% k2 = 10^x(12);
% km2 = 10^x(13);
% k3= 10^x(14);
% km3= 10^x(15);
% k4= 10^x(16);
% k5 = 10^x(17);
% km5 = 10^x(18);
% k6 = 10^x(19);
% kb2 = 10^x(24);
% kd2 = 10^x(25);
% k2b = k2*fvc(1);
% km2b = km2*fvc(2);
% k1b = k1*fvc(3);
% km1b =km1*fvc(4);
% kb1 = kb2*fvc(5);
% kd1 = kd2*fvc(6);
% k3b = k3*fvc(7);
% km3b = km3*fvc(8);
% k4b = k4*fvc(9); 
% k5b = k5*fvc(10);
% km5b =km5*fvc(11);
% k6b = k6*fvc(12);
% kmdeg =10^x(2);
% ktlnA =10^x(3); %*10^x(2);
% ktlnE = 10^x(4); %*10^x(2);
% kbtpn2 = 10^x(6);
% kbtpn1 = 0.47*10^x(6);
% f1 = 10^x(7);
% K1 = 10^x(8);
% K2 =10^x(20);
% f2 = 10^x(21);
% KMT = (km3+k4+kpdeg_eff)/k3;KMTb = (km3b+k4b+kpdeg_eff)/k3b;
% KMP = (km5+k6+kpdeg_eff)/k5;KMPb = (km5b+k6b+kpdeg_eff)/k5b;
% 
% 
% 
% RRP=Y(2);    %PhoP~P
% prP = prd*(ktlnA*(kbtpn1*((1+f1*(RRP/K1)^2)/(1+(RRP/K1)^2)))/kmdeg); prQ = prP/40;
% prB = prd*ktlnE*10^x(9)/kmdeg; %prd*(ktlnE*(kbtpn2*((1+f2*(RRP/K2)^2)/(1+(RRP/K2)^2)))/kmdeg);
% X = prB/kpdeg_eff; 
% Pt = prP/kpdeg_eff;
% RR= Pt - RRP;%Y(1);     %PhoP
% SKkin=Y(3);    %PhoQ* -kin
% % X = Y(6); %MgrB
% SKXph = Y(9);   %PhoQ.MgrB -ph
% 
% SKPRR = (SKkin)*(RR/KMT)/(1+(1+k4/k2)*RR/KMT);
% SKs = SKPRR*k4/k2;
% SKph= SKs*km1/(k1+kb2*X+kpdeg_eff); %Y(4);     %PhoQ - ph 
% SKXkin = SKs*kb2*X/(kd2+km1+kpdeg_eff);%Y(7);  %PhoQ*.MgrB -kin
% 
% SKXRRP = SKXph*(RRP/KMPb)/(1+RRP/KMPb);
% SKRRP = SKph*(RRP/KMP)/(1+RRP/KMP);
% SK = SKph-SKRRP;
% SKXs = SKXkin;
% SKX = SKXph-SKXRRP;
% 
% 
% %---------------
% dY(2) = (k4*SKPRR - k6*SKRRP - k6b*SKXRRP); % PhoP~P
% dY(3) = 1*prQ + k1*SK - km1*SKs - kb2*SKs*X + kd2*SKXs - kpdeg_eff*SKkin; % PhoQ* chk
% dY(4) = + km1*SKs - k1*SK - kb1*SK*X+ kd1*SKX - kpdeg_eff*SKph; % PhoQ chk
% % dY(6) = prB -kb2*SKs*X+kd2*SKXs -kb1*SK*X + kd1*SKX - kpdeg_eff*X; %MgrB
% dY(7) = kb2*SKs*X - kd2*SKXs +k1b*SKX - km1b*SKXs - k2b*SKXs - kpdeg_eff*SKXkin; %QB* - [PhoQ*.MgrB] chk
% dY(9) = kb1*SK*X - kd1*SKX - k1b*SKX + km1b*SKXs - kpdeg_eff*SKXph; %QB [PhoQ.MgrB] chk
% dY=dY';
% end