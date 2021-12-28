clear
clc

Folder='InitialState_MediumSystem_Stoichiometry/';
mkdir(Folder)

fig=0;
Rescale=3;
BoxSize(1)=300;
BoxSize(2)=50;
BoxSize(3)=50;

L1=8; 
L2=8;
Ratio=1:0.04:1.32;%[1,8/7,8/6];

NR=length(Ratio);
    
for nr=1:NR
        
    ratio=Ratio(nr);
    
    NP=2500;
    NP1=round(NP*ratio/(ratio+1)/L1);
    NP2=round(NP/(ratio+1)/L2);

    NB=14;
    x=[-L1/2 L1/2];
    y=(1:NB)-mean(1:NB);
    z=(1:NB)-mean(1:NB);

    [X,Y,Z]=meshgrid(x,y,z);
    X=reshape(X,[],1);
    Y=reshape(Y,[],1);
    Z=reshape(Z,[],1);

    nx=length(x);
    ny=length(y);
    nz=length(z);

    Polymer1=zeros(3,L1,nx*ny*nz);

    for np=1:nx*ny*nz
        Polymer1(1,:,np)=(0:L1-1)+X(np,1);
        Polymer1(2,:,np)=0+Y(np,1);
        Polymer1(3,:,np)=0+Z(np,1);
    end

    Polymer2=Polymer1(:,1:L2,:);
    Polymer1=Polymer1(:,:,1:NP1);
    Polymer2=Polymer2(:,:,1:NP2);
    Polymer1=Polymer1*Rescale;
    Polymer2=Polymer2*Rescale;

    fig=fig+1;
    figure(fig)
    for np=1:NP1
        plot3(Polymer1(1,:,np),Polymer1(2,:,np),Polymer1(3,:,np),'b.-'); hold on
    end
    for np=1:NP2
        plot3(Polymer2(1,:,np),Polymer2(2,:,np),Polymer2(3,:,np),'ro-'); hold on
    end
    axis equal
    axis([-BoxSize(1)/2 BoxSize(1)/2 -BoxSize(2)/2 BoxSize(2)/2 -BoxSize(3)/2 BoxSize(3)/2])
    
    NM=L1*NP1+L2*NP2;
    Monomer=zeros(3,NM);
    NB=(L1-1)*NP1+(L2-1)*NP2;
    Bond=zeros(2,NB);
    Atype=zeros(1,NM);

    nm=0;
    nb=0;
    for np=1:NP1
        for l=1:L1
            nm=nm+1;
            Monomer(1,nm)=Polymer1(1,l,np);
            Monomer(2,nm)=Polymer1(2,l,np);
            Monomer(3,nm)=Polymer1(3,l,np);
            Atype(1,nm)=1;
            if l~=L1
                nb=nb+1;
                Bond(1,nb)=nm;
                Bond(2,nb)=nm+1;
            end
        end
    end
    for np=1:NP2
        for l=1:L2
            nm=nm+1;
            Monomer(1,nm)=Polymer2(1,l,np);
            Monomer(2,nm)=Polymer2(2,l,np);
            Monomer(3,nm)=Polymer2(3,l,np);
            Atype(1,nm)=2;
            if l~=L2
                nb=nb+1;
                Bond(1,nb)=nm;
                Bond(2,nb)=nm+1;
            end
        end
    end

    save([Folder 'L1_' num2str(L1) '_L2_' num2str(L2) ...
                '_N1_' num2str(NP1) '_N2_' num2str(NP2) '.mat'],'Monomer','Bond','Atype','BoxSize');
end

