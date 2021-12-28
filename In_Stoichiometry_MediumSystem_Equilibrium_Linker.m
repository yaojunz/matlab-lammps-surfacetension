clear
clc

global Monomer Bond Atype BoxSize BeadSize BeadCsi BeadMass Temp kBT Damp

RandStream.setGlobalStream(RandStream('mt19937ar','seed',sum(100*clock)));

LoadFolder='Parameter/Parameter.mat';
load(LoadFolder);

BeadMass=Damp*BeadCsi;

Replicates=10;

rcut=1.5;

A=6;
% DeltaF=A2DeltaF(A,rcut);
% Tau=12*exp(DeltaF);

for Linker=2
    InFolder=['MediumSystem_Stoichiometry/In_Linker' num2str(Linker) '/'];
    mkdir(InFolder);
    OutFolder='Out_Record1/';
    mkdir([InFolder OutFolder]);
    index=0;

    L1=8; 
    L2=8;
    Ratio=1:0.04:1.32;
    NR=length(Ratio);
    
    for nr=1:NR
        ratio=Ratio(nr);
        NP=2500;
        np1=round(NP*ratio/(ratio+1)/L1);
        np2=round(NP/(ratio+1)/L2);
    
        LoadFolder='InitialState_MediumSystem_Stoichiometry/';
        load([LoadFolder 'L1_' num2str(L1) '_L2_' num2str(L2) ...
                        '_N1_' num2str(np1) '_N2_' num2str(np2) '.mat']);

        for rep=(1:Replicates)
            index=index+1;
            Filename=['L1_' num2str(L1) '_L2_' num2str(L2) '_N1_' num2str(np1) '_N2_' num2str(np2) '_A' num2str(A) '_Rep' num2str(rep)];
            InitCondFilename=['L1_' num2str(L1) '_L2_' num2str(L2) '_N1_' num2str(np1) '_N2_' num2str(np2) '.initial'];
            if rep==1
                InitCondGenerate(InFolder,InitCondFilename)
            end
            InFilename=['Record1_Index_' num2str(index) '.in'];
            OutFilename=[OutFolder Filename];
            InFileGenerate(InFolder,InFilename,InitCondFilename,OutFilename,A,rep,Linker)
        end
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
%     K=0.2*kBT;  % Lp=0.8;
%     R0=10.5;% Lc=0.38*30; 
%     nbin=10;
% elseif Linker==3
    K=0.15*kBT; % Lp=0.8;
    R0=14;  % Lc=0.38*40;
    nbin=13;
end
Delta=sigma;
Rb=60;

TimeStep=Damp/100;
RunSteps=10^8;%50*round(Tau)/TimeStep;
Thermo=RunSteps;

fid=fopen([InFolder InFilename],'w');
fprintf(fid, ['units nano \n']);
fprintf(fid, ['boundary p p p \n']);
fprintf(fid, ['atom_style bond\n\n']);

fprintf(fid, ['processors 2 * *\n\n']);

fprintf(fid, ['read_data ' InitCondFilename '\n\n']);

fprintf(fid, ['pair_style hybrid lj/cut ' num2str(sigma*2^(1/6)) ' soft ' num2str(Rc) '\n']);
fprintf(fid, ['pair_coeff * * lj/cut ' num2str(eps) ' ' num2str(sigma) ' ' num2str(sigma*2^(1/6)) '\n']);
fprintf(fid, ['pair_coeff 1 2 soft ' num2str(Ai) ' ' num2str(Rc) '\n']);
fprintf(fid, ['pair_modify shift yes\n']);
fprintf(fid, ['special_bonds lj/coul 1.0 1.0 1.0\n\n']);

fprintf(fid, ['bond_style fene/expand\n']);
fprintf(fid, ['bond_coeff * ' num2str(K) ' ' num2str(R0) ' ' num2str(0) ' ' num2str(0) ' ' num2str(Delta) '\n\n']);

fprintf(fid, ['variable A equal "ramp(' num2str(Ai) ',' num2str(Af) ')"\n\n']);

fprintf(fid, ['region wallx block -' num2str(Rb) ' ' num2str(Rb) ...
                                ' -' num2str(100) ' ' num2str(100) ...
                                ' -' num2str(100) ' ' num2str(100) ' open 3 open 4 open 5 open 6\n\n']);

fprintf(fid, ['neighbor ' num2str(nbin) ' bin \n']);
fprintf(fid, ['neigh_modify every 1 delay 0\n\n']);

fprintf(fid, ['fix 1 all nve\n']);
fprintf(fid, ['fix 2 all langevin ' num2str(Temp) ' ' num2str(Temp) ' ' num2str(Damp) ' ' num2str(randi(10^7)) ' zero no\n']);
fprintf(fid, ['fix 3 all adapt 1 pair soft a 1 2 v_A\n']);
fprintf(fid, ['fix 4 all wall/region wallx lj126 ' num2str(kBT) ' ' num2str(BeadSize) ' ' num2str(BeadSize*2^(1/6)) '\n']);
fprintf(fid, ['fix 5 all balance ' num2str(RunSteps/400) ' 1.05 shift x 10 1.05\n\n']);

fprintf(fid, ['thermo ' num2str(Thermo) '\n']);
fprintf(fid, ['timestep ' num2str(TimeStep) '\n\n']);

if rep==1
fprintf(fid, ['dump 1 all movie ' num2str(RunSteps/100) ' ' OutFilename '.mpeg type type zoom 4.5 box yes 0.01 view 85 85 size 1000 400 shiny 0.5\n']);
fprintf(fid, ['dump_modify 1 acolor 1 blue\n']); 
fprintf(fid, ['dump_modify 1 acolor 2 orange\n']); 
fprintf(fid, ['dump_modify 1 adiam 1 3\n']); 
fprintf(fid, ['dump_modify 1 adiam 2 3\n\n']);
end

fprintf(fid, ['run ' num2str(RunSteps/10) '\n']);
fprintf(fid, ['unfix 3\n']);
fprintf(fid, ['run ' num2str(RunSteps/10) '\n']);
fprintf(fid, ['unfix 4\n']);
fprintf(fid, ['run ' num2str(RunSteps) '\n\n']);

fprintf(fid, ['dump 2 all xyz ' num2str(RunSteps/100) ' ' OutFilename '.xyz\n']); 
fprintf(fid, ['run ' num2str(RunSteps) '\n\n']);
fprintf(fid, ['write_restart ' OutFilename '.restart\n']);

fclose(fid); 

end

function []=InitCondGenerate(InFolder,InitCondFilename)

global Monomer Bond Atype BoxSize BeadMass

Natom_type=2; %number of atom types;
Natom=length(Monomer);
Nbond_type=1; %number of bond types;
Nbond=length(Bond); %number of DNA bonds;
Nangl_type=0;  %number of angle types;
Nangl=0; %number of angles;
Ndihe_type=0; %number of dihedral types;
Ndihe=0; %number of dihedrals;
Nimpr_type=0; %number of improper types;
Nimpr=0; %number of impropers;

xlo=-BoxSize(1)/2; xhi=BoxSize(1)/2; %x boundary
ylo=-BoxSize(2)/2; yhi=BoxSize(2)/2; %y boundary
zlo=-BoxSize(3)/2; zhi=BoxSize(3)/2; %z boundary

V=zeros(3,Natom);

fid=fopen([InFolder InitCondFilename],'w');
fprintf(fid,'LAMMPS chain data file\n\n');
fprintf(fid,'%d atoms\n', Natom);
fprintf(fid,'%d bonds\n', Nbond);
fprintf(fid,'%d angles\n', Nangl);
fprintf(fid,'%d dihedrals\n', Ndihe);
fprintf(fid,'%d impropers\n\n', Nimpr);
fprintf(fid,'%d atom types\n', Natom_type);
fprintf(fid,'%d bond types\n', Nbond_type);
fprintf(fid,'%d angle types\n', Nangl_type);
fprintf(fid,'%d dihedral types\n', Ndihe_type);
fprintf(fid,'%d improper types\n\n', Nimpr_type);
fprintf(fid,'%8.5f %8.5f xlo xhi\n', xlo, xhi);
fprintf(fid,'%8.5f %8.5f ylo yhi\n', ylo, yhi);
fprintf(fid,'%8.5f %8.5f zlo zhi\n\n', zlo, zhi);

fprintf(fid,'Masses\n\n');
fprintf(fid,'%d %8.5f\n',1,BeadMass);
fprintf(fid,'%d %8.5f\n\n',2,BeadMass);

fprintf(fid,'Atoms\n\n');
for i=1:Natom
    Mole_type=Atype(i);
    Atom_type=Atype(i);
    fprintf(fid,[num2str(i) ' ' num2str(Mole_type) ' ' num2str(Atom_type) ' ' ...
                 num2str(Monomer(1,i)) ' ' num2str(Monomer(2,i)) ' ' num2str(Monomer(3,i)) '\n']);
end
fprintf(fid,'\n');

fprintf(fid,'Velocities\n\n');
for i=1:Natom
    fprintf(fid,[num2str(i) ' ' num2str(V(1,i)) ' ' num2str(V(2,i)) ' ' num2str(V(3,i)) '\n']);
end
fprintf(fid,'\n');

fprintf(fid,'Bonds\n\n');
for i=1:Nbond
    Bond_type=1;
    fprintf(fid,'%d %d %d %d\n',...
            i,Bond_type,Bond(1,i),Bond(2,i));
end
fprintf(fid,'\n');
fclose(fid); 
end

function DeltaF=A2DeltaF(A,rcut)
dr=10^-4;
r=dr/2:dr:(rcut-dr/2);
DeltaF=zeros(length(A),1);
for na=1:length(A)
    DeltaF(na,1)=log(sum(4*pi*r.^2.*exp(U_Soft(r,A(na),rcut)))*dr/sum(4*pi*r.^2*dr));
end
end

function y=U_Soft(r,A,rcut)
y=A*(1+cos(r*pi/rcut)).*(r<rcut);
end