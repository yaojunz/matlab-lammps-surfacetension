clear
clc

global Monomer Bond Atype BoxSize BeadSize BeadCsi BeadMass Temp kBT Damp

RandStream.setGlobalStream(RandStream('mt19937ar','seed',sum(100*clock)));

LoadFolder='Parameter/Parameter.mat';
load(LoadFolder);

BeadMass=Damp*BeadCsi;

Replicates=10;
A=6;

InFolder=['MediumSystem_Stoichiometry/In_Linker2/'];
OutFolder=['Out_SurfaceTension/'];
mkdir([InFolder OutFolder]);

index=0;
L1=8; 
L2=8;
Ratio=1:0.04:1.32;
NR=length(Ratio);
Linker=2;
for nr=1:NR
   	ratio=Ratio(nr);
  	NP=2500;
   	np1=round(NP*ratio/(ratio+1)/L1);
   	np2=round(NP/(ratio+1)/L2);
  	for rep=(1:Replicates)
    	index=index+1;
     	Filename=['L1_' num2str(L1) '_L2_' num2str(L2) '_N1_' num2str(np1) '_N2_' num2str(np2) '_A' num2str(A) '_Rep' num2str(rep)];
      	InitCondFilename=['Out_Record4/L1_' num2str(L1) '_L2_' num2str(L2) '_N1_' num2str(np1) '_N2_' num2str(np2) '_A' num2str(A) '_Rep' num2str(rep) '.restart'];
        InFilename=['SurfaceTension_Index_' num2str(index) '.in'];
      	OutFilename=[OutFolder Filename];
        InFileGenerate(InFolder,InFilename,InitCondFilename,OutFilename,A,rep,Linker)
    end
end


function []=InFileGenerate(InFolder,InFilename,InitCondFilename,OutFilename,A,rep,Linker)

global Monomer Bond Atype BoxSize BeadSize BeadCsi BeadMass Temp kBT Damp
Ai=-0*kBT;
Af=-A*kBT;
eps=kBT;
sigma=BeadSize;
Rc=BeadSize/2;
if Linker==1
    K=0.3*kBT;  % Lp=0.8;
    R0=7;   % Lc=0.38*20;
    nbin=7;
elseif Linker==2
    K=0.15*kBT; % Lp=0.8;
    R0=14;  % Lc=0.38*40;
    nbin=13;
end
Delta=sigma;
Rb=60;

TimeStep=Damp/100;
RunSteps=4*10^8;
Thermo=RunSteps;

fid=fopen([InFolder InFilename],'w');
fprintf(fid, ['processors 2 * *\n\n']);

fprintf(fid, ['read_restart ' InitCondFilename '\n\n']);

fprintf(fid, ['pair_style hybrid lj/cut ' num2str(sigma*2^(1/6)) ' soft ' num2str(Rc) '\n']);
fprintf(fid, ['pair_coeff * * lj/cut ' num2str(eps) ' ' num2str(sigma) ' ' num2str(sigma*2^(1/6)) '\n']);
fprintf(fid, ['pair_coeff 1 2 soft ' num2str(Af) ' ' num2str(Rc) '\n']);
fprintf(fid, ['pair_modify shift yes\n']);
fprintf(fid, ['special_bonds lj/coul 1.0 1.0 1.0\n\n']);

fprintf(fid, ['bond_style fene/expand\n']);
fprintf(fid, ['bond_coeff * ' num2str(K) ' ' num2str(R0) ' ' num2str(0) ' ' num2str(0) ' ' num2str(Delta) '\n\n']);

fprintf(fid, ['neighbor ' num2str(nbin) ' bin \n']);
fprintf(fid, ['neigh_modify every 1 delay 0\n\n']);

fprintf(fid, ['fix 1 all nve\n']);
fprintf(fid, ['fix 2 all langevin ' num2str(Temp) ' ' num2str(Temp) ' ' num2str(Damp) ' ' num2str(randi(10^7)) ' zero no\n']);
fprintf(fid, ['fix 3 all balance ' num2str(RunSteps/400) ' 1.05 shift x 10 1.05\n\n']);

fprintf(fid, ['thermo ' num2str(Thermo) '\n']);
fprintf(fid, ['timestep ' num2str(TimeStep) '\n\n']);

fprintf(fid, ['reset_timestep 0\n']);
fprintf(fid, ['dump 2 all xyz ' num2str(RunSteps/400) ' ' OutFilename '.xyz\n']);
fprintf(fid, ['fix 4 all ave/time 1 ' num2str(RunSteps/10) ' ' num2str(RunSteps/10) ' c_thermo_press[1] c_thermo_press[2] c_thermo_press[3] file ' OutFilename '.dump\n']);
fprintf(fid, ['run ' num2str(RunSteps) '\n\n']);

fprintf(fid, ['write_restart ' OutFilename '.restart\n\n']);

fclose(fid); 

end
