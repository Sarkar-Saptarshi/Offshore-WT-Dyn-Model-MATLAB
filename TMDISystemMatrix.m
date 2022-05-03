function [C_comb, f_comb] = TMDISystemMatrix(Comq, q_Nom, Controls, ElastoDyn, Twr, Bld, TMDI)
g = 9.80665;
% q_Sg   = Comq(1);
% q_Sw   = Comq(2);
% q_Hv   = Comq(3);
% q_R    = Comq(4);
% q_P    = Comq(5);
% q_Y    = Comq(6);
q_TFA1 = Comq(7);
q_TSS1 = Comq(8);
q_TFA2 = Comq(9);
q_TSS2 = Comq(10);
% q_yaw  = Comq(11);
% q_GeAz = Comq(12);
% q_DrTr = Comq(13);
% q_B1F1 = Comq(14);
% q_B1E1 = Comq(15);
% q_B1F2 = Comq(16);
% q_B2F1 = Comq(17);
% q_B2E1 = Comq(18);
% q_B2F2 = Comq(19);
% q_B3F1 = Comq(20);
% q_B3E1 = Comq(21);
% q_B3F2 = Comq(22);
q_TMD  = Comq(23);
 
% qd_Sg   = Comq(24);
% qd_Sw   = Comq(25);
% qd_Hv   = Comq(26);
qd_R    = Comq(27);
qd_P    = Comq(28);
qd_Y    = Comq(29);
qd_TFA1 = Comq(30);
qd_TSS1 = Comq(31);
qd_TFA2 = Comq(32);
qd_TSS2 = Comq(33);
% qd_yaw  = Comq(34);
% qd_GeAz = Comq(35);

% qd_DrTr = Comq(36);
% qd_B1F1 = Comq(37);
% qd_B1E1 = Comq(38);
% qd_B1F2 = Comq(39);
% qd_B2F1 = Comq(40);
% qd_B2E1 = Comq(41);
% qd_B2F2 = Comq(42);
% qd_B3F1 = Comq(43);
% qd_B3E1 = Comq(44);
% qd_B3F2 = Comq(45);
qd_TMD  = Comq(46);

% TMDI contribution to system matrices. To create these matrices
C_TMD       = zeros(23,23);
C_Inert     = zeros(23,23);
f_InertTwr  = zeros(23,1);
f_InertDamp = zeros(23,1);
f_ElasticTMD = zeros(23,1);
f_DampTMD    = zeros(23,1);
f_GravTMD    = zeros(23,1);
f_TMD        = zeros(23,1);

%% Read parameters
b        = TMDI.b;        
mTMD     = TMDI.mTMD;     
f1_TMD   = TMDI.f1_TMD;    
k_TMD    = TMDI.k_TMD;
zeta_TMD = TMDI.zeta_TMD;                 % same damping ratio with or without inerter
Idx_In   = TMDI.Idx;
H_I      = TMDI.H_I;
O1_TFA   = Twr.O1_TFA;
O2_TFA   = Twr.O2_TFA; 
O1_TSS   = Twr.O1_TSS;
O2_TSS   = Twr.O2_TSS;
s11_TFA = Twr.s11_TFA;
s12_TFA = Twr.s12_TFA;
s22_TFA = Twr.s22_TFA;
S11_TFA = Twr.S11_TFA;
S12_TFA = Twr.S12_TFA;
S22_TFA = Twr.S22_TFA;
s11_TSS = Twr.s11_TSS;
s12_TSS = Twr.s12_TSS;
s22_TSS = Twr.s22_TSS;
S11_TSS = Twr.S11_TSS;
S12_TSS = Twr.S12_TSS;
S22_TSS = Twr.S22_TSS;

TwrHt = ElastoDyn.TwrHt;
TowerBsHt = ElastoDyn.TowerBsHt;  

PtfmRefzt = ElastoDyn.PtfmRefzt;   

%% Coordinate systems
BlPitch = Controls(2:4);

[Z, A, D, ~, ~, ~, ~, ~, ~, ~, ~] = Coordinate_systems(q_Nom, BlPitch, ElastoDyn, Twr, Bld);

%% Position vectors
rZO_1 = q_TFA1 + q_TFA2;                                                     % S11_TFA, S22_TF etc are scalars
rZO_2 = PtfmRefzt + TwrHt - 0.5*(S11_TFA*q_TFA1^2 + S22_TFA ...
           *q_TFA2^2 +2*S12_TFA*q_TFA1*q_TFA2 + S11_TSS*q_TSS1^2 + S22_TSS*q_TSS2^2 ...
            + 2*S12_TSS*q_TSS1*q_TSS2);
rZO_3 = q_TSS1 + q_TSS2;
rZO   = rZO_1*A(1,:) + rZO_2*A(2,:) + rZO_3*A(3,:); 

rZHI_1 = O1_TFA(Idx_In)*q_TFA1 + O2_TFA(Idx_In)*q_TFA2;                                                     % S11_TFA, S22_TF etc are scalars
rZHI_2 = H_I + PtfmRefzt + TowerBsHt - 0.5*(s11_TFA(Idx_In)*q_TFA1^2 + s22_TFA(Idx_In) ...
           *q_TFA2^2 +2*s12_TFA(Idx_In)*q_TFA1*q_TFA2 + s11_TSS(Idx_In)*q_TSS1^2 + s22_TSS(Idx_In)*q_TSS2^2 ...
            + 2*s12_TSS(Idx_In)*q_TSS1*q_TSS2);
rZHI_3 = O1_TSS(Idx_In)*q_TSS1 + O1_TSS(Idx_In)*q_TSS2;
rZHI   = rZHI_1*A(1,:) + rZHI_2*A(2,:) + rZHI_3*A(3,:); 

rHIO   = rZO - rZHI;

% if strcmp(TMDI.Dir, 'SS')
%     rOD = q_TMD*A(3,:);
% elseif strcmp(TMD.Dir, 'FA')
%     rOD = q_TMDI*A(1,:);
% end

%% Velocity vectors
EwX = qd_R*Z(1,:) + qd_Y*Z(2,:) - qd_P*Z(3,:);

XvHI_1 = (1-O1_TFA(Idx_In))*qd_TFA1 + (1-O2_TFA(Idx_In))*qd_TFA2;
XvHI_2 = (S11_TFA - s11_TFA(Idx_In))*q_TFA1*qd_TFA1 + (S22_TFA - s22_TFA(Idx_In))*q_TFA2*qd_TFA2 ...
         + (S12_TFA - s12_TFA(Idx_In))*(qd_TFA1*q_TFA2 + q_TFA1*qd_TFA2) + (S11_TSS - s11_TSS(Idx_In))*q_TSS1*qd_TSS1 + (S22_TSS - s22_TSS(Idx_In))*q_TSS2*qd_TSS2 ...
         + (S12_TSS - s12_TSS(Idx_In))*(qd_TSS1*q_TSS2 + q_TSS1*qd_TSS2);
XvHI_3 = (1-O1_TSS(Idx_In))*qd_TSS1 + (1-O2_TSS(Idx_In))*qd_TSS2;

XvHI   = XvHI_1*A(1,:) - XvHI_2*A(2,:) + XvHI_3*A(3,:);

XvO_1 = qd_TFA1 + qd_TFA2;
XvO_2 = -(S11_TFA*q_TFA1*qd_TFA1 + S22_TFA*q_TFA2*qd_TFA2 + S12_TFA*(qd_TFA1*q_TFA2+q_TFA1*qd_TFA2) ...
        + S11_TSS*q_TSS1*qd_TSS1 + S22_TSS*q_TSS2*qd_TSS2 + S12_TSS*(qd_TSS1*q_TSS2+q_TSS1*qd_TSS2));
XvO_3 = qd_TSS1 + qd_TSS2;
XvO = XvO_1*A(1,:) + XvO_2*A(2,:) + XvO_3*A(3,:);

%% Partial velocities 
% Platform
EwX_R = Z(1,:);
EwX_P = -Z(3,:);
EwX_Y = Z(2,:);

% Inerter hook
EvHI_Sg = Z(1,:);
EvHI_Sw = -Z(3,:);
EvHI_Hv = Z(2,:);
EvHI_R  = cross(EwX_R,rZHI);
EvHI_P  = cross(EwX_P,rZHI);
EvHI_Y  = cross(EwX_Y,rZHI);
EvHI_TFA1 = O1_TFA(Idx_In)*A(1,:) - (s11_TFA(Idx_In)*q_TFA1 + s12_TFA(Idx_In)*q_TFA2)*A(2,:); 
EvHI_TSS1 = O1_TSS(Idx_In)*A(3,:) - (s11_TSS(Idx_In)*q_TSS1 + s12_TSS(Idx_In)*q_TSS2)*A(2,:);
EvHI_TFA2 = O2_TFA(Idx_In)*A(1,:) - (s22_TFA(Idx_In)*q_TFA2 + s12_TFA(Idx_In)*q_TFA1)*A(2,:);
EvHI_TSS2 = O2_TSS(Idx_In)*A(3,:) - (s22_TSS(Idx_In)*q_TSS2 + s12_TSS(Idx_In)*q_TSS1)*A(2,:);

% Tower top/Base plate
EvO_Sg = Z(1,:);
EvO_Sw = -Z(3,:);
EvO_Hv = Z(2,:);
EvO_R  = cross(EwX_R,rZO);
EvO_P  = cross(EwX_P,rZO);
EvO_Y  = cross(EwX_Y,rZO);
EvO_TFA1 = A(1,:) - (S11_TFA*q_TFA1 + S12_TFA*q_TFA2)*A(2,:);
EvO_TSS1 = A(3,:) - (S11_TSS*q_TSS1 + S12_TSS*q_TSS2)*A(2,:);
EvO_TFA2 = A(1,:) - (S22_TFA*q_TFA2 + S12_TFA*q_TFA1)*A(2,:);
EvO_TSS2 = A(3,:) - (S22_TSS*q_TSS2 + S12_TSS*q_TSS1)*A(2,:);

% TMD
EvD_Sg = EvO_Sg;
EvD_Sw = EvO_Sw;
EvD_Hv = EvO_Hv;
EvD_R  = EvO_R;
EvD_P  = EvO_P;
EvD_Y  = EvO_Y;
EvD_TFA1  = EvO_TFA1;
EvD_TSS1  = EvO_TSS1;
EvD_TFA2  = EvO_TFA2;
EvD_TSS2  = EvO_TSS2;

EvD_TMD  = zeros(1,3);
if strcmp(TMDI.Dir, 'SS')
    EvD_TMD = D(3,:);
elseif strcmp(TMDI.Dir, 'FA')
    EvD_TMD = D(1,:);
end

% Partial accelerations of TMD
dEvO_R = cross(EwX_R, XvO + cross(EwX,rZO));
dEvO_P = cross(EwX_P, XvO + cross(EwX,rZO));
dEvO_Y = cross(EwX_Y, XvO + cross(EwX,rZO));
dEvO_TFA1 = -(S11_TFA*qd_TFA1 + S12_TFA*qd_TFA2)*A(2,:) + cross(EwX, EvO_TFA1);
dEvO_TSS1 = -(S11_TSS*qd_TSS1 + S12_TSS*qd_TSS2)*A(2,:) + cross(EwX, EvO_TSS1);
dEvO_TFA2 = -(S22_TFA*qd_TFA2 + S12_TFA*qd_TFA1)*A(2,:) + cross(EwX, EvO_TFA2);
dEvO_TSS2 = -(S22_TSS*qd_TSS2 + S12_TSS*qd_TSS1)*A(2,:) + cross(EwX, EvO_TSS2);

dEvD_R = dEvO_R;
dEvD_P = dEvO_P;
dEvD_Y = dEvO_Y;
dEvD_TFA1 = dEvO_TFA1;
dEvD_TSS1 = dEvO_TSS1;
dEvD_TFA2 = dEvO_TFA2;
dEvD_TSS2 = dEvO_TSS2;
dEvD_TMD  = cross(EwX, EvD_TMD);

dEvD = dEvD_R*qd_R + dEvD_P*qd_P + dEvD_Y*qd_Y + dEvD_TFA1*qd_TFA1 + dEvD_TSS1*qd_TSS1 + dEvD_TFA2*qd_TFA2 + dEvD_TSS2*qd_TSS2 + dEvD_TMD*qd_TMD;
   
%% Partial TMDI forces
F_R    = b*cross(EwX_R, rHIO);
F_P    = b*cross(EwX_P, rHIO);
F_Y    = b*cross(EwX_Y, rHIO);
F_TFA1 = b*(  (1-O1_TFA(Idx_In))*A(1,:) - ( (S11_TFA-s11_TFA(Idx_In))*q_TFA1 + (S12_TFA - s12_TFA(Idx_In))*q_TFA2 )*A(2,:)  );
F_TSS1 = b*(  (1-O1_TSS(Idx_In))*A(3,:) - ( (S11_TSS-s11_TSS(Idx_In))*q_TSS1 + (S12_TSS - s12_TSS(Idx_In))*q_TSS2 )*A(2,:)  );
F_TFA2 = b*(  (1-O2_TFA(Idx_In))*A(1,:) - ( (S22_TFA-s22_TFA(Idx_In))*q_TFA2 + (S12_TFA - s12_TFA(Idx_In))*q_TFA1 )*A(2,:)  );
F_TSS2 = b*(  (1-O2_TSS(Idx_In))*A(1,:) - ( (S22_TSS-s22_TSS(Idx_In))*q_TSS2 + (S12_TSS - s12_TSS(Idx_In))*q_TSS1 )*A(2,:)  );
F_TMDI = b*(  EvD_TMD  );

%% Derivative of partial forces
dF_R    = b*(  cross(EwX_R, XvHI + cross(EwX, rHIO))  );
dF_P    = b*(  cross(EwX_P, XvHI + cross(EwX, rHIO))  );
dF_Y    = b*(  cross(EwX_Y, XvHI + cross(EwX, rHIO))  );
dF_TFA1 = b*( - ( (S11_TFA-s11_TFA(Idx_In))*qd_TFA1 + (S12_TFA - s12_TFA(Idx_In))*qd_TFA2 )*A(2,:)  ) + cross(EwX, F_TFA1);
dF_TSS1 = b*( - ( (S11_TSS-s11_TSS(Idx_In))*qd_TSS1 + (S12_TSS - s12_TSS(Idx_In))*qd_TSS2 )*A(2,:)  ) + cross(EwX, F_TSS1); 
dF_TFA2 = b*( - ( (S22_TFA-s22_TFA(Idx_In))*qd_TFA2 + (S12_TFA - s12_TFA(Idx_In))*qd_TFA1 )*A(2,:)  ) + cross(EwX, F_TFA2);
dF_TSS2 = b*( - ( (S22_TSS-s22_TSS(Idx_In))*qd_TSS2 + (S12_TSS - s12_TSS(Idx_In))*qd_TSS1 )*A(2,:)  ) + cross(EwX, F_TSS2);
dF_TMDI = cross(EwX, F_TMDI);

dFSum   = dF_R*qd_R + dF_P*qd_P + dF_Y*qd_Y + dF_TFA1*qd_TFA1 + dF_TSS1*qd_TSS1 + dF_TFA2*qd_TFA2 + dF_TSS2*qd_TSS2 + dF_TMDI*qd_TMD;

dFSum_D = dot(dFSum, EvD_TMD)*EvD_TMD;

%% Tower TMD

C_TMD(1,1)  = mTMD*dot(EvD_Sg,EvD_Sg);
C_TMD(1,2)  = mTMD*dot(EvD_Sg,EvD_Sw);
C_TMD(1,3)  = mTMD*dot(EvD_Sg,EvD_Hv);
C_TMD(1,4)  = mTMD*dot(EvD_Sg,EvD_R);
C_TMD(1,5)  = mTMD*dot(EvD_Sg,EvD_P);
C_TMD(1,6)  = mTMD*dot(EvD_Sg,EvD_Y);
C_TMD(1,7)  = mTMD*dot(EvD_Sg,EvD_TFA1);
C_TMD(1,8)  = mTMD*dot(EvD_Sg,EvD_TSS1);
C_TMD(1,9)  = mTMD*dot(EvD_Sg,EvD_TFA2);
C_TMD(1,10) = mTMD*dot(EvD_Sg,EvD_TSS2);
C_TMD(1,23) = mTMD*dot(EvD_Sg,EvD_TMD);

C_TMD(2,2)  = mTMD*dot(EvD_Sw,EvD_Sw);
C_TMD(2,3)  = mTMD*dot(EvD_Sw,EvD_Hv);
C_TMD(2,4)  = mTMD*dot(EvD_Sw,EvD_R);
C_TMD(2,5)  = mTMD*dot(EvD_Sw,EvD_P);
C_TMD(2,6)  = mTMD*dot(EvD_Sw,EvD_Y);
C_TMD(2,7)  = mTMD*dot(EvD_Sw,EvD_TFA1);
C_TMD(2,8)  = mTMD*dot(EvD_Sw,EvD_TSS1);
C_TMD(2,9)  = mTMD*dot(EvD_Sw,EvD_TFA2);
C_TMD(2,10) = mTMD*dot(EvD_Sw,EvD_TSS2);
C_TMD(2,23) = mTMD*dot(EvD_Sw,EvD_TMD);

C_TMD(3,3)  = mTMD*dot(EvD_Hv,EvD_Hv);
C_TMD(3,4)  = mTMD*dot(EvD_Hv,EvD_R);
C_TMD(3,5)  = mTMD*dot(EvD_Hv,EvD_P);
C_TMD(3,6)  = mTMD*dot(EvD_Hv,EvD_Y);
C_TMD(3,7)  = mTMD*dot(EvD_Hv,EvD_TFA1);
C_TMD(3,8)  = mTMD*dot(EvD_Hv,EvD_TSS1);
C_TMD(3,9)  = mTMD*dot(EvD_Hv,EvD_TFA2);
C_TMD(3,10) = mTMD*dot(EvD_Hv,EvD_TSS2);
C_TMD(3,23) = mTMD*dot(EvD_Hv,EvD_TMD);

C_TMD(4,4)  = mTMD*dot(EvD_R,EvD_R);
C_TMD(4,5)  = mTMD*dot(EvD_R,EvD_P);
C_TMD(4,6)  = mTMD*dot(EvD_R,EvD_Y);
C_TMD(4,7)  = mTMD*dot(EvD_R,EvD_TFA1);
C_TMD(4,8)  = mTMD*dot(EvD_R,EvD_TSS1);
C_TMD(4,9)  = mTMD*dot(EvD_R,EvD_TFA2);
C_TMD(4,10) = mTMD*dot(EvD_R,EvD_TSS2);
C_TMD(4,23) = mTMD*dot(EvD_R,EvD_TMD);

C_TMD(5,5)  = mTMD*dot(EvD_P,EvD_P);
C_TMD(5,6)  = mTMD*dot(EvD_P,EvD_Y);
C_TMD(5,7)  = mTMD*dot(EvD_P,EvD_TFA1);
C_TMD(5,8)  = mTMD*dot(EvD_P,EvD_TSS1);
C_TMD(5,9)  = mTMD*dot(EvD_P,EvD_TFA2);
C_TMD(5,10) = mTMD*dot(EvD_P,EvD_TSS2);
C_TMD(5,23) = mTMD*dot(EvD_P,EvD_TMD);

C_TMD(6,6)  = mTMD*dot(EvD_Y,EvD_Y);
C_TMD(6,7)  = mTMD*dot(EvD_Y,EvD_TFA1);
C_TMD(6,8)  = mTMD*dot(EvD_Y,EvD_TSS1);
C_TMD(6,9)  = mTMD*dot(EvD_Y,EvD_TFA2);
C_TMD(6,10) = mTMD*dot(EvD_Y,EvD_TSS2);
C_TMD(6,23) = mTMD*dot(EvD_Y,EvD_TMD);

C_TMD(7,7)  = mTMD*dot(EvD_TFA1,EvD_TFA1);
C_TMD(7,8)  = mTMD*dot(EvD_TFA1,EvD_TSS1);
C_TMD(7,9)  = mTMD*dot(EvD_TFA1,EvD_TFA2);
C_TMD(7,10) = mTMD*dot(EvD_TFA1,EvD_TSS2);
C_TMD(7,23) = mTMD*dot(EvD_TFA1,EvD_TMD); 

C_TMD(8,8)  = mTMD*dot(EvD_TSS1,EvD_TSS1);
C_TMD(8,9)  = mTMD*dot(EvD_TSS1,EvD_TFA2);
C_TMD(8,10) = mTMD*dot(EvD_TSS1,EvD_TSS2);
C_TMD(8,23) = mTMD*dot(EvD_TSS1,EvD_TMD);

C_TMD(9,9)  = mTMD*dot(EvD_TFA2,EvD_TFA2);
C_TMD(9,10) = mTMD*dot(EvD_TFA2,EvD_TSS2);
C_TMD(9,23) = mTMD*dot(EvD_TFA2,EvD_TMD);

C_TMD(10,10) = mTMD*dot(EvD_TSS2,EvD_TSS2);
C_TMD(10,23) = mTMD*dot(EvD_TSS2,EvD_TMD);

C_TMD(23,23) = mTMD*dot(EvD_TMD,EvD_TMD);

C_TMD  = triu(C_TMD)+triu(C_TMD,1)';

f_TMD(1) = -mTMD*dot(EvD_Sg,dEvD);
f_TMD(2) = -mTMD*dot(EvD_Sw,dEvD);
f_TMD(3) = -mTMD*dot(EvD_Hv,dEvD);
f_TMD(4) = -mTMD*dot(EvD_R,dEvD);
f_TMD(5) = -mTMD*dot(EvD_P,dEvD);
f_TMD(6) = -mTMD*dot(EvD_Y,dEvD);
f_TMD(7) = -mTMD*dot(EvD_TFA1,dEvD);
f_TMD(8) = -mTMD*dot(EvD_TSS1,dEvD);
f_TMD(9) = -mTMD*dot(EvD_TFA2,dEvD);
f_TMD(10) = -mTMD*dot(EvD_TSS2,dEvD);
f_TMD(23) = -mTMD*dot(EvD_TMD,dEvD);

f_GravTMD(1) = -mTMD*g*dot(EvD_Sg,Z(2,:)); 
f_GravTMD(2) = -mTMD*g*dot(EvD_Sw,Z(2,:));
f_GravTMD(3) = -mTMD*g*dot(EvD_Hv,Z(2,:));
f_GravTMD(4) = -mTMD*g*dot(EvD_R,Z(2,:));
f_GravTMD(5) = -mTMD*g*dot(EvD_P,Z(2,:));
f_GravTMD(6) = -mTMD*g*dot(EvD_Y,Z(2,:));
f_GravTMD(7) = -mTMD*g*dot(EvD_TFA1,Z(2,:));
f_GravTMD(8) = -mTMD*g*dot(EvD_TSS1,Z(2,:));
f_GravTMD(9) = -mTMD*g*dot(EvD_TFA2,Z(2,:));
f_GravTMD(10) = -mTMD*g*dot(EvD_TSS2,Z(2,:));
f_GravTMD(23) = -mTMD*g*dot(EvD_TMD,Z(2,:));

f_ElasticTMD(23) = -k_TMD*q_TMD;

f_DampTMD(23)    = -zeta_TMD*k_TMD/pi/f1_TMD*qd_TMD;

% end of EOM for Tower Side-to-Side TMD

%% Inerter C matrix
C_Inert(4, 4)  = dot ( EvD_R - EvHI_R, dot(F_R, EvD_TMD)   *EvD_TMD);
C_Inert(4, 5)  = dot ( EvD_R - EvHI_R, dot(F_P, EvD_TMD)   *EvD_TMD);
C_Inert(4, 6)  = dot ( EvD_R - EvHI_R, dot(F_Y, EvD_TMD)   *EvD_TMD);
C_Inert(4, 7)  = dot ( EvD_R - EvHI_R, dot(F_TFA1, EvD_TMD)*EvD_TMD);
C_Inert(4, 8)  = dot ( EvD_R - EvHI_R, dot(F_TSS1, EvD_TMD)*EvD_TMD);
C_Inert(4, 9)  = dot ( EvD_R - EvHI_R, dot(F_TFA2, EvD_TMD)*EvD_TMD);
C_Inert(4, 10) = dot ( EvD_R - EvHI_R, dot(F_TSS2, EvD_TMD)*EvD_TMD);
C_Inert(4, 23) = dot ( EvD_R - EvHI_R, dot(F_TMDI, EvD_TMD)*EvD_TMD);

C_Inert(5, 4)  = dot ( EvD_P - EvHI_P, dot(F_R, EvD_TMD)   *EvD_TMD);
C_Inert(5, 5)  = dot ( EvD_P - EvHI_P, dot(F_P, EvD_TMD)   *EvD_TMD);
C_Inert(5, 6)  = dot ( EvD_P - EvHI_P, dot(F_Y, EvD_TMD)   *EvD_TMD);
C_Inert(5, 7)  = dot ( EvD_P - EvHI_P, dot(F_TFA1, EvD_TMD)*EvD_TMD);
C_Inert(5, 8)  = dot ( EvD_P - EvHI_P, dot(F_TSS1, EvD_TMD)*EvD_TMD);
C_Inert(5, 9)  = dot ( EvD_P - EvHI_P, dot(F_TFA2, EvD_TMD)*EvD_TMD);
C_Inert(5, 10) = dot ( EvD_P - EvHI_P, dot(F_TSS2, EvD_TMD)*EvD_TMD);
C_Inert(5, 23) = dot ( EvD_P - EvHI_P, dot(F_TMDI, EvD_TMD)*EvD_TMD);

C_Inert(6, 4)  = dot ( EvD_Y - EvHI_Y, dot(F_R, EvD_TMD)   *EvD_TMD);
C_Inert(6, 5)  = dot ( EvD_Y - EvHI_Y, dot(F_P, EvD_TMD)   *EvD_TMD);
C_Inert(6, 6)  = dot ( EvD_Y - EvHI_Y, dot(F_Y, EvD_TMD)   *EvD_TMD);
C_Inert(6, 7)  = dot ( EvD_Y - EvHI_Y, dot(F_TFA1, EvD_TMD)*EvD_TMD);
C_Inert(6, 8)  = dot ( EvD_Y - EvHI_Y, dot(F_TSS1, EvD_TMD)*EvD_TMD);
C_Inert(6, 9)  = dot ( EvD_Y - EvHI_Y, dot(F_TFA2, EvD_TMD)*EvD_TMD);
C_Inert(6, 10) = dot ( EvD_Y - EvHI_Y, dot(F_TSS2, EvD_TMD)*EvD_TMD);
C_Inert(6, 23) = dot ( EvD_Y - EvHI_Y, dot(F_TMDI, EvD_TMD)*EvD_TMD);

C_Inert(7, 4)  = dot ( EvD_TFA1 - EvHI_TFA1, dot(F_R, EvD_TMD)   *EvD_TMD);
C_Inert(7, 5)  = dot ( EvD_TFA1 - EvHI_TFA1, dot(F_P, EvD_TMD)   *EvD_TMD);
C_Inert(7, 6)  = dot ( EvD_TFA1 - EvHI_TFA1, dot(F_Y, EvD_TMD)   *EvD_TMD);
C_Inert(7, 7)  = dot ( EvD_TFA1 - EvHI_TFA1, dot(F_TFA1, EvD_TMD)*EvD_TMD);
C_Inert(7, 8)  = dot ( EvD_TFA1 - EvHI_TFA1, dot(F_TSS1, EvD_TMD)*EvD_TMD);
C_Inert(7, 9)  = dot ( EvD_TFA1 - EvHI_TFA1, dot(F_TFA2, EvD_TMD)*EvD_TMD);
C_Inert(7, 10) = dot ( EvD_TFA1 - EvHI_TFA1, dot(F_TSS2, EvD_TMD)*EvD_TMD);
C_Inert(7, 23) = dot ( EvD_TFA1 - EvHI_TFA1, dot(F_TMDI, EvD_TMD)*EvD_TMD);

C_Inert(8, 4)  = dot ( EvD_TSS1 - EvHI_TSS1, dot(F_R, EvD_TMD)   *EvD_TMD);
C_Inert(8, 5)  = dot ( EvD_TSS1 - EvHI_TSS1, dot(F_P, EvD_TMD)   *EvD_TMD);
C_Inert(8, 6)  = dot ( EvD_TSS1 - EvHI_TSS1, dot(F_Y, EvD_TMD)   *EvD_TMD);
C_Inert(8, 7)  = dot ( EvD_TSS1 - EvHI_TSS1, dot(F_TFA1, EvD_TMD)*EvD_TMD);
C_Inert(8, 8)  = dot ( EvD_TSS1 - EvHI_TSS1, dot(F_TSS1, EvD_TMD)*EvD_TMD);
C_Inert(8, 9)  = dot ( EvD_TSS1 - EvHI_TSS1, dot(F_TFA2, EvD_TMD)*EvD_TMD);
C_Inert(8, 10) = dot ( EvD_TSS1 - EvHI_TSS1, dot(F_TSS2, EvD_TMD)*EvD_TMD);
C_Inert(8, 23) = dot ( EvD_TSS1 - EvHI_TSS1, dot(F_TMDI, EvD_TMD)*EvD_TMD);

C_Inert(9, 4)  = dot ( EvD_TFA2 - EvHI_TFA2, dot(F_R, EvD_TMD)   *EvD_TMD);
C_Inert(9, 5)  = dot ( EvD_TFA2 - EvHI_TFA2, dot(F_P, EvD_TMD)   *EvD_TMD);
C_Inert(9, 6)  = dot ( EvD_TFA2 - EvHI_TFA2, dot(F_Y, EvD_TMD)   *EvD_TMD);
C_Inert(9, 7)  = dot ( EvD_TFA2 - EvHI_TFA2, dot(F_TFA1, EvD_TMD)*EvD_TMD);
C_Inert(9, 8)  = dot ( EvD_TFA2 - EvHI_TFA2, dot(F_TSS1, EvD_TMD)*EvD_TMD);
C_Inert(9, 9)  = dot ( EvD_TFA2 - EvHI_TFA2, dot(F_TFA2, EvD_TMD)*EvD_TMD);
C_Inert(9, 10) = dot ( EvD_TFA2 - EvHI_TFA2, dot(F_TSS2, EvD_TMD)*EvD_TMD);
C_Inert(9, 23) = dot ( EvD_TFA2 - EvHI_TFA2, dot(F_TMDI, EvD_TMD)*EvD_TMD);

C_Inert(10, 4)  = dot ( EvD_TSS2 - EvHI_TSS2, dot(F_R, EvD_TMD)   *EvD_TMD);
C_Inert(10, 5)  = dot ( EvD_TSS2 - EvHI_TSS2, dot(F_P, EvD_TMD)   *EvD_TMD);
C_Inert(10, 6)  = dot ( EvD_TSS2 - EvHI_TSS2, dot(F_Y, EvD_TMD)   *EvD_TMD);
C_Inert(10, 7)  = dot ( EvD_TSS2 - EvHI_TSS2, dot(F_TFA1, EvD_TMD)*EvD_TMD);
C_Inert(10, 8)  = dot ( EvD_TSS2 - EvHI_TSS2, dot(F_TSS1, EvD_TMD)*EvD_TMD);
C_Inert(10, 9)  = dot ( EvD_TSS2 - EvHI_TSS2, dot(F_TFA2, EvD_TMD)*EvD_TMD);
C_Inert(10, 10) = dot ( EvD_TSS2 - EvHI_TSS2, dot(F_TSS2, EvD_TMD)*EvD_TMD);
C_Inert(10, 23) = dot ( EvD_TSS2 - EvHI_TSS2, dot(F_TMDI, EvD_TMD)*EvD_TMD);

C_Inert(23, 4)  = dot ( EvD_TMD, dot(F_R, EvD_TMD)   *EvD_TMD);
C_Inert(23, 5)  = dot ( EvD_TMD, dot(F_P, EvD_TMD)   *EvD_TMD);
C_Inert(23, 6)  = dot ( EvD_TMD, dot(F_Y, EvD_TMD)   *EvD_TMD);
C_Inert(23, 7)  = dot ( EvD_TMD, dot(F_TFA1, EvD_TMD)*EvD_TMD);
C_Inert(23, 8)  = dot ( EvD_TMD, dot(F_TSS1, EvD_TMD)*EvD_TMD);
C_Inert(23, 9)  = dot ( EvD_TMD, dot(F_TFA2, EvD_TMD)*EvD_TMD);
C_Inert(23, 10) = dot ( EvD_TMD, dot(F_TSS2, EvD_TMD)*EvD_TMD);
C_Inert(23, 23) = dot ( EvD_TMD, dot(F_TMDI, EvD_TMD)*EvD_TMD);

% disp(C_Inert([4:10,23],[4:10,23]))
% pause
%% Inerter force vectors
f_InertDamp(1)  = -dot(EvD_Sg, dFSum_D);
f_InertDamp(2)  = -dot(EvD_Sw, dFSum_D);
f_InertDamp(3)  = -dot(EvD_Hv, dFSum_D);
f_InertDamp(4)  = -dot(EvD_R, dFSum_D);
f_InertDamp(5)  = -dot(EvD_P, dFSum_D);
f_InertDamp(6)  = -dot(EvD_Y, dFSum_D);
f_InertDamp(7)  = -dot(EvD_TFA1, dFSum_D);
f_InertDamp(8)  = -dot(EvD_TSS1, dFSum_D);
f_InertDamp(9)  = -dot(EvD_TFA2, dFSum_D);
f_InertDamp(10) = -dot(EvD_TSS2, dFSum_D);
f_InertDamp(23) = -dot(EvD_TMD, dFSum_D);

f_InertTwr(1)  = dot(EvHI_Sg, dFSum_D);
f_InertTwr(2)  = dot(EvHI_Sw, dFSum_D);
f_InertTwr(3)  = dot(EvHI_Hv, dFSum_D);
f_InertTwr(4)  = dot(EvHI_R, dFSum_D);
f_InertTwr(5)  = dot(EvHI_P, dFSum_D);
f_InertTwr(6)  = dot(EvHI_Y, dFSum_D);
f_InertTwr(7)  = dot(EvHI_TFA1, dFSum_D);
f_InertTwr(8)  = dot(EvHI_TSS1, dFSum_D);
f_InertTwr(9)  = dot(EvHI_TFA2, dFSum_D);
f_InertTwr(10) = dot(EvHI_TSS2, dFSum_D);

%% Assembly
C_comb = C_TMD + C_Inert;
f_comb = f_TMD + f_GravTMD + f_ElasticTMD + f_DampTMD + f_InertDamp + f_InertTwr;


end

