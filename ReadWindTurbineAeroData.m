
function [Airfoils, Geometry] = ReadWindTurbineAeroData(angle)
% function required to be run before the BEM Simulink block can be used.
% This script reads the FAST input files obtain the blade geometrical 
% and aerodynamic properties.

% Enter the path of the folders that contains the files.
% Read airfoil data
if strcmp(angle,'deg')
    fac = 1;
elseif strcmp(angle, 'rad')
    fac = pi/180;
end

fid = fopen('./5MW_Baseline/AeroData/Cylinder1.dat','r');               % Airfoil 1
formatSpec = repmat('%f', 1, 4);
D = textscan(fid,formatSpec,'HeaderLines',14);
fclose(fid);
Cylinder1 = cell2mat(D);
Cylinder1(:,1) = Cylinder1(:,1)*fac;
Airfoils.Cylinder1Lift = griddedInterpolant(Cylinder1(:,1),Cylinder1(:,2),'linear');
Airfoils.Cylinder1Drag = griddedInterpolant(Cylinder1(:,1),Cylinder1(:,3),'linear');

fid = fopen('./5MW_Baseline/AeroData/Cylinder2.dat','r');               % Airfoil 2
formatSpec = repmat('%f', 1, 4);
D = textscan(fid,formatSpec,'HeaderLines',14);
fclose(fid);
Cylinder2 = cell2mat(D);
Cylinder2(:,1) = Cylinder2(:,1)*fac;
Airfoils.Cylinder2Lift = griddedInterpolant(Cylinder2(:,1),Cylinder2(:,2),'linear');
Airfoils.Cylinder2Drag = griddedInterpolant(Cylinder2(:,1),Cylinder2(:,3),'linear');

fid = fopen('./5MW_Baseline/AeroData/DU21_A17.dat','r');                % Airfoil 3
formatSpec = repmat('%f', 1, 4);
D = textscan(fid,formatSpec,'HeaderLines',14);
fclose(fid);
DU21_A17 = cell2mat(D);
DU21_A17(:,1) = DU21_A17(:,1)*fac;
Airfoils.DU21_A17Lift = griddedInterpolant(DU21_A17(:,1),DU21_A17(:,2),'linear');
Airfoils.DU21_A17Drag = griddedInterpolant(DU21_A17(:,1),DU21_A17(:,3),'linear');

fid = fopen('./5MW_Baseline/AeroData/DU25_A17.dat','r');                % Airfoil 4
formatSpec = repmat('%f', 1, 4);
D = textscan(fid,formatSpec,'HeaderLines',14);
fclose(fid);
DU25_A17 = cell2mat(D);
DU25_A17(:,1) = DU25_A17(:,1)*fac;
Airfoils.DU25_A17Lift = griddedInterpolant(DU25_A17(:,1),DU25_A17(:,2),'linear');
Airfoils.DU25_A17Drag = griddedInterpolant(DU25_A17(:,1),DU25_A17(:,3),'linear');

fid = fopen('./5MW_Baseline/AeroData/DU30_A17.dat','r');                % Airfoil 5
formatSpec = repmat('%f', 1, 4);
D = textscan(fid,formatSpec,'HeaderLines',14);
fclose(fid);
DU30_A17 = cell2mat(D);
DU30_A17(:,1) = DU30_A17(:,1)*fac;
Airfoils.DU30_A17Lift = griddedInterpolant(DU30_A17(:,1),DU30_A17(:,2),'linear');
Airfoils.DU30_A17Drag = griddedInterpolant(DU30_A17(:,1),DU30_A17(:,3),'linear');

fid = fopen('./5MW_Baseline/AeroData/DU35_A17.dat','r');                % Airfoil 6
formatSpec = repmat('%f', 1, 4);
D = textscan(fid,formatSpec,'HeaderLines',14);
fclose(fid);
DU35_A17 = cell2mat(D);
DU35_A17(:,1) = DU35_A17(:,1)*fac;
Airfoils.DU35_A17Lift = griddedInterpolant(DU35_A17(:,1),DU35_A17(:,2),'linear');
Airfoils.DU35_A17Drag = griddedInterpolant(DU35_A17(:,1),DU35_A17(:,3),'linear');

fid = fopen('./5MW_Baseline/AeroData/DU40_A17.dat','r');                % Airfoil 7
formatSpec = repmat('%f', 1, 4);
D = textscan(fid,formatSpec,'HeaderLines',14);
fclose(fid);
DU40_A17 = cell2mat(D);
DU40_A17(:,1) = DU40_A17(:,1)*fac;
Airfoils.DU40_A17Lift = griddedInterpolant(DU40_A17(:,1),DU40_A17(:,2),'linear');
Airfoils.DU40_A17Drag = griddedInterpolant(DU40_A17(:,1),DU40_A17(:,3),'linear');

fid = fopen('./5MW_Baseline/AeroData/NACA64_A17.dat','r');              % Airfoil 8
formatSpec = repmat('%f', 1, 4);
D = textscan(fid,formatSpec,'HeaderLines',14);
fclose(fid);
NACA64_A17 = cell2mat(D);
NACA64_A17(:,1) = NACA64_A17(:,1)*fac;
Airfoils.NACA64_A17Lift = griddedInterpolant(NACA64_A17(:,1),NACA64_A17(:,2),'linear');
Airfoils.NACA64_A17Drag = griddedInterpolant(NACA64_A17(:,1),NACA64_A17(:,3),'linear');

% Read Blade geometry
fid = fopen('./5MW_Baseline/AeroData/NRELOffshrBsline5MW_AeroDyn_blade.dat');
Info = textscan(fid,'%f %s','delimiter','\n','HeaderLines',3);
formatSpec = repmat('%f', 1, 7);
D = textscan(fid,formatSpec,Info{1},'HeaderLines',2);
fclose(fid);
Geometry.Blade = cell2mat(D);
Geometry.Blade(:,5) = Geometry.Blade(:,5)*fac;
% Geometry.Blade(1,:) = [];
% Geometry.Blade(end,:) = [];
Geometry.Blade(:,1) = Geometry.Blade(:,1);


% Read Tower geometry
fid = fopen('./5MW_Baseline/AeroData/NRELOffshrBsline5MW_Onshore_AeroDyn15.dat');
Info = textscan(fid,'%f %s','delimiter','\n','HeaderLines',48);
formatSpec = repmat('%f', 1, 3);
D = textscan(fid,formatSpec,Info{1},'HeaderLines',2);
Geometry.Tower= cell2mat(D);
fclose(fid);











