clear
clc

% mass = attograms = 10^-18 gram
% distance = nanometers
% time = nanoseconds
% energy = attogram-nanometer^2/nanosecond^2
% velocity = nanometers/nanosecond
% force = attogram-nanometer/nanosecond^2
% torque = attogram-nanometer^2/nanosecond^2
% temperature = Kelvin
% pressure = attogram/(nanometer-nanosecond^2)
% dynamic viscosity = attogram/(nanometer-nanosecond)
% charge = multiple of electron charge (1.0 is a proton)
% dipole = charge-nanometer
% electric field = volt/nanometer
% density = attograms/nanometer^dim

BeadSize=3;
WaterEta=1; %ag/nm/ns %0.001 kg/m/s
BeadCsi=6*pi*WaterEta*BeadSize/2; %ag/ns

Temp=300; %K
kB=1.38*10^-2; %ag*nm^2/ns^2/K 1.38*10^-23 Kg*m^2/s^2/K
kBT=kB*Temp;

Damp=20;
D=kBT/BeadCsi;
BeadMass=Damp*BeadCsi;

dt=Damp/100;
T=dt:0.1:10^4;
loglog(T,3*kBT/BeadMass*T.^2,'-'); hold on
loglog(T,6*D*T,'-');

SaveFolder='Parameter/';
mkdir(SaveFolder);
save([SaveFolder 'Parameter.mat'],'BeadSize','BeadCsi','Damp','Temp','kBT');
