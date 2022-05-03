% TCD Onshore Wind Turbine Code

clc; clearvars -except RMSStore; clear -globals;
clear SystemMatrices BaselineControllers;
pause(.5);
global GRAVACC FLUIDDENSITY;      
GRAVACC = 9.80655; FLUIDDENSITY = 1025;

if libisloaded('MoorDyn') || libisloaded('MoorApiwin64')
    if libisloaded('MoorApiwin64')
    calllib('MoorApiwin64','finish'); unloadlibrary('MoorApiwin64');
    else
    calllib('MoorDyn','LinesClose');      unloadlibrary MoorDyn;  
    end                                    % unload library (never forget to do this!)end
end

% Read FAST input data files
[Airfoils, Geometry] = ReadWindTurbineAeroDataInterp('rad');
[Blade, Tower]       = ReadWindTurbineStructuralData();
[ElastoDyn]          = ReadElastoDyn();
[Servo]              = ServoDyn();

% Create tower and blade structures
Twr = CreateTwr(Tower,ElastoDyn); Bld = CreateBld(ElastoDyn,Geometry,Blade); Platform = CreatePlatform();
Platform.Mooring     = 1;    % 1 for Moordyn, 2 for OpenMoor
Platform.WaveLoads   = 1;    % 1 to calculate wave load on spar using Morisson's equation, 0 to not.

% DOFs available
DOFsStr = {'Sg','Sw','Hv','R','P','Y','TFA1','TSS1','TFA2','TSS2','NacYaw','GeAz','DrTr','B1F1','B1E1','B1F2','B2F1','B2E1','B2F2','B3F1','B3E1','B3F2','TMD'};

if Platform.WaveLoads == 1
    waveopt.wave_profile = 'load';
    waveopt.wave_file = 'WaveHs0_75_Tp6_Dir0_NoCur0';
    wave = Wave(waveopt);   % Define Parameters for TMDI
else
    wave.wave_file = [];
end
m_dummy_TSS1 = 3.8507e+05;

TMDI.Switch    = true(1);   % false -> when TMD is off, true -> when TMD is on
if TMDI.Switch
    TMDI.MassRatio = 1/100;
    TMDI.beta      = .4;
else
    TMDI.MassRatio = 0/100;
    TMDI.beta      = 0/100;
end
TMDI.Idx       = 8;
TMDI.Dir       = 'FA';
if strcmp(TMDI.Dir, 'SS')
    TMDI.phi  = Twr.O1_TSS(TMDI.Idx);
elseif strcmp(TMDI.Dir, 'FA')
    TMDI.phi  = Twr.O1_TFA(TMDI.Idx);
end
[TMDI.TuningRat, TMDI.zeta_TMD] = OptTMDIParams(TMDI.MassRatio, TMDI.beta, TMDI.phi);
TMDI.b         = TMDI.beta*m_dummy_TSS1;
TMDI.mTMD      = TMDI.MassRatio*m_dummy_TSS1;
TMDI.f1_TMD    = 0.4885*TMDI.TuningRat;
TMDI.k_TMD     = (TMDI.f1_TMD*2*pi)^2*(TMDI.mTMD);
%TMDI.TMDzn     = (TMDI.mTMD/7850*3/pi/4)^(1/3);
TMDI.H_I       = Twr.TwrSec(TMDI.Idx);

% Defined the DOFs that are turned off
TurnedOffDOFsStr = {'TMD'};
if TMDI.Switch
    TurnedOffDOFsStr = {};
end
TurnedOffDOFs   = find(ismember(DOFsStr,TurnedOffDOFsStr));
DOFs.Avail      = length(DOFsStr);
DOFs.Active     = setdiff(1:length(DOFsStr),TurnedOffDOFs);
DOFs.ActNominal = setdiff(1:22,TurnedOffDOFs);
DOFs.nDOFs      = length(DOFs.Active);

nDOFs = length(DOFs.Active);
if find(ismember(TurnedOffDOFsStr,'GeAz'))
    q0 = zeros(1,2*nDOFs);
else
    q0 = zeros(1,2*nDOFs);
    q0((12 - sum(TurnedOffDOFs < 12)) + nDOFs) = ElastoDyn.RotSpeed;
end

% Get the TurbSim generated wind field and grid
[velocity, Wind.y, Wind.z, Wind.nz, Wind.ny, Wind.dz,...
          Wind.dy, Wind.dt, Wind.zHub, Wind.z1, Wind.SummVars] = readBLgrid('NTMClassC8_TI17_4');
Wind.t_TurbSim = (0:1:size(velocity,1)-1)*Wind.dt;

% if Steady
% [Wind.t_TurbSim, velocity] = SteadyWindNew(Wind.y, Wind.z, 11.4, 0.2);

gv = {Wind.t_TurbSim,[1 2 3], Wind.y, Wind.z};
Wind.Velocity = griddedInterpolant(gv,velocity,'linear');

WindNom.PittandPeters = true(0);
WindNom.PitchControl  = true(1);

WindNom.y = Wind.y;
WindNom.z = Wind.z;

% Initialize inflow angles
phi = InitializeInflowAngle(ElastoDyn,Bld,Wind);

% Mooring Lines Initialization
mooring_load_ptr = zeros(1,6);           % going to make a pointer so LinesCalc can modify FLines
mooring_load_ptr = libpointer('doublePtr',mooring_load_ptr);  % access returned value with FLines_p.value

if Platform.Mooring == 1
    [notfound,warnings] = loadlibrary('MoorDyn','MoorDyn');     % load MoorDyn DLL
    mooring_status = calllib('MoorDyn','LinesInit',zeros(6,1),zeros(6,1));   % initialize MoorDyn 
    if mooring_status ~= 0
        error('MoorDyn unable to load');
    else
        disp('MoorDyn loaded successfully');
    end
elseif Platform.Mooring == 2
    loadlibrary('MoorApiwin64','moorapi.h'); input_file = 'CaseOC3.xml'; calllib('MoorApiwin64','initialize',input_file);
else
    error('Choose a proper mooring dynamics model');
end
%               === X === X ===
Controls0 = [Servo.VS_RtTq, ElastoDyn.BlPitch, reshape(phi, [1 3*Bld.nb])];

Lin = 0;
t0  = 0;
tf  = 20;
deltat = 0.0125;
t = t0:deltat:tf;
n = (tf-t0)/deltat;
[q, Controls] = ode4(@(t,q,Controls) RHS(t, q, Controls, DOFs, ElastoDyn, Airfoils, Twr, Bld,...
                Platform, Wind, WindNom, wave, mooring_load_ptr, Servo, TMDI),t0,deltat,tf,q0,Controls0);

if Platform.Mooring == 1
    calllib('MoorDyn','LinesClose');       % close MoorDyn
    unloadlibrary MoorDyn;                 % unload library (never forget to do this!)
elseif Platform.Mooring == 2
    pause(1)
    calllib('MoorApiwin64','finish');
    unloadlibrary('MoorApiwin64');
end

xmin = 100/deltat +1;

[OopDefl1, IPDefl1, TipDzb1, TTDspFA, TTDspSS, LSSTipV, NacYaw, TwrClrnB1, TwrClrnB2,...
                    TwrClrnB3, BldTwrHrClr, Surge, Sway, Heave, Roll, Pitch, Yaw] = Results( t, q, Controls, DOFs, ElastoDyn, Twr, Bld );

figure(1)
hold on
plot(t,OopDefl1)
xlabel('Time')
ylabel('OopDefl1')

% figure(2)
% hold on
% plot(t,IPDefl1)
% xlabel('Time')
% ylabel('IPDefl1')
% 
% 
% figure(3)
% hold on
% plot(t,LSSTipV)
% xlabel('Time')
% ylabel('LSSTipV')
% 
% figure(4)
% hold on
% plot(t,NacYaw)
% xlabel('Time')
% ylabel('NacYaw')
 
figure(5)
hold on
plot(t,TTDspFA)
xlabel('Time')
ylabel('TTDspFA')
 
% figure(6)
% hold on
% plot(t,TTDspSS)
% xlabel('Time (sec)','FontSize',12)
% ylabel('Tower side-to-side displacement (m)','FontSize',12)
 
% figure(7)
% hold on
% plot(t, Surge)
% xlabel('Time')
% ylabel('Platform Surge')
% 
% figure(8)
% hold on
% plot(t, Sway)
% xlabel('Time')
% ylabel('Platform Sway')
% 
% figure(9)
% hold on
% plot(t, Heave)
% xlabel('Time')
% ylabel('Platform Heave')
% 
% figure(10)
% hold on
% plot(t, Roll)
% xlabel('Time')
% ylabel('Platform Roll')
% 
% figure(11)
% hold on
% plot(t, Pitch)
% xlabel('Time')
% ylabel('Platform Pitch')
% 
% figure(12)
% hold on
% plot(t, Yaw)
% xlabel('Time')
% ylabel('Platform Yaw')

% figure(13)
% plot(t, TwrClrnB1, t, TwrClrnB2, t, TwrClrnB3)
% xlabel('Time')
% ylabel('Tower Blade Clearence')
% TwrRad = interp1([0 87.6],[6 3.87],87.6-(63-(90-87.6)))/2;
% hline(TwrRad,'--r','Tower Radius')
% 
% figure(14)
% hold on
% plot(t, Controls(:,2), t, Controls(:,3), t, Controls(:,4))
% xlabel('Time')
% ylabel('Blade pitch angles')

% figure(15)
% hold on
% [f,Amp] = FS(t(xmin:end),TTDspFA(xmin:end)); plot(f,Amp); xlim([0 1.2]); % ylim([0 .6]);
% xlabel('Frequency (Hz)')
% ylabel('TTDspFA (f)')

% figure(16)
% hold on
% plot(t,q(:,23));
% xlabel('Time')
% ylabel('TMD stroke')
