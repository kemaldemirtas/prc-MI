% clear all; clc;

%% PDB reading:



fname1='6wx4.pdb';  %pdb adını yaz.


% Normal PDB için: 
[d1, d2, d3, d4, d5 ,d6, x_1, y_1, z_1, d7, d8, d9]=textread(fname1,'%s%f%s%s%s%f%f%f%f%f%f%s','headerlines',1);
 
A=strcmp(d3,'CA'); %To grep CA atoms only.

x_1=x_1(A);
y_1=y_1(A);
z_1=z_1(A);

x_1 = x_1(1:2:640);
y_1 = y_1(1:2:640);
z_1 = z_1(1:2:640);

x=[x_1']; 
y=[y_1'];
z=[z_1'];   


%  general
resnum=size(x,2); 
frame=size(x,1);

%% GNM Parameters
rcut_gnm=7.5; %A
ga=1;
t0=6;

%% preallocate arrays
A=zeros(resnum,resnum);
U=zeros(resnum,resnum);
S=zeros(resnum,resnum);
V=zeros(resnum,resnum);
w=zeros(resnum,frame);

%% GNM Calculation
for j=1:resnum
    for k=1:resnum
         distx = x(1,j)-x(1,k);   
         disty = y(1,j)-y(1,k);   
         distz = z(1,j)-z(1,k);
         
         r=sqrt(distx^2+disty^2+distz^2);
        % Kirchhoff
        if (r <= rcut_gnm && j~=k && r > 0.0001)
            A(j,k)=-1;
        else
            A(j,k)=0;                                
        end
    end
end

% detailed balance for connectivity matrix
diagonal=sum(A(:,:));
for j=1:resnum
    for i=1:resnum
        if i == j
            A(i,j)=-1*diagonal(i); % Kirchhoff Connectivity
        end
    end
end
    


%% --------------------------------------------------------------------

[U(:,:),S(:,:),V(:,:)]=svd(A(:,:)); %singular value decomposition
winit=diag(S(:,:));  %eigenvalues        
S1=pinv(S);
w=diag(S1(:,:));  %1/eigenvalue


      
m1=1;  %first mode
m2=319; %second mode

for j=1:resnum
        for i=1:resnum
            invcont_av10(i,j)=0;
            for     k=resnum-m2:resnum-m1
            invcont_av10(i,j)=invcont_av10(i,j)+U(i,k)*V(j,k)*w(k);      
            end
        end
end



%% GNM cross correlation normalization (C_ij)
for j=1:resnum
    for i=1:resnum
        crosscorr_av10(i,j)=invcont_av10(i,j)/(sqrt(invcont_av10(i,i)*invcont_av10(j,j)));
    end
end



%% mutual information calculation
for j=1:resnum
    for i=1:resnum
        kesir(i,j)=(invcont_av10(i,j)*invcont_av10(i,j))/(invcont_av10(i,i)*invcont_av10(j,j));
        ic(i,j)=1-kesir(i,j);
        mutual_information(i,j)=-0.5*log(ic(i,j));
    end
end

mutinf=real(mutual_information); 


%% MI MAP:
f1=figure
f1=imagesc(1:resnum,1:resnum,mutinf); 
set(gca,'YDir','normal');
set(gca,'FontSize',14, 'FontWeight', 'normal');
set(gca,'FontSize',24)
axis square;

C=colorbar;
caxis([0 0.01]); %can be changed wrt the molecule!!!
colormap(jet)
