% ========================================================================
% Introduction
% ========================================================================
% This code provides a simple demonstration of steady optical beam (STOB).
%
% Author: Yiqian Yang (yang-yq22@mails.tsinghua.edu.cn)
% =========================================================================

%% 
clc;
clear all;
%%%%

L         = 0.08;                     
N         = 1024;                       
x         = linspace(-L/2,L/2,N)';
y         = x;
dx        = x(2)-x(1);                   
[x,y]     = meshgrid(x,y);
wx        = [0:N/2-1  -N/2:-1]'*2*pi/L;
wy        = [0:N/2-1  -N/2:-1]'*2*pi/L;
[wx,wy]   = meshgrid(wx,wy);

lz        = 2;                
Nz        = 500;                         
z         = linspace(0,lz,Nz);
dz        = z(2)-z(1);                     
         
%%%%
lamda      = 532*10^(-9);               
k0         = 2*pi/lamda;                

aa  =1;                              
r=sqrt(x.*x+y.*y);
circ = @(x,y,rx,ry,fat)sqrt((x-rx).^2+(y-ry).^2) <= fat;

fat0=0.006;   
fat1=0.004; 
Gau_ini_1=aa*circ(x,y,0,0,fat0).*exp(-1i*k0*(fat1-r)*0.003); 

fat2=0.016;
Gau_ini_2=aa*exp(1i*k0*(fat1-r)*0.003).*(circ(x,y,0,0,fat2)-circ(x,y,0,0,0.008)); 
k=k0*(fat1-r)*0.003;
 
%%%%
A1=Gau_ini_2+Gau_ini_1;
figure(1);
imagesc(abs(A1).^2);
axis square;axis off;colormap(hot);

%%%%
LO         = exp(-1i*(wx.^2+wy.^2)*dz /(2*k0));  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% propagation
u1        = A1;
I1        =abs(u1).^2;
% figure(2); imagesc(I1);
% axis equal
pic=zeros(1024,1024,512);
pic(:,:,1)=I1;
EvI(:,1)  =I1(N/2,:);

for j = 2:Nz
	um        = ifft2(LO.*fft2(u1));
	NO        = exp(1i.* dz.*k0);
	u1        = um.* NO;
    I=abs(u1).^2;
    pic(:,:,j)=I;
    EvI(:,j)=I(N/2,:);

    pause(0.01)
    figure(5);imagesc(abs(u1).^2);colormap(hot);axis tight;axis off;axis square
end

figure(12);
subplot(2,4,1)
imagesc(pic(:,:,150));axis off;colormap(hot);
subplot(2,4,2)
imagesc(pic(:,:,250));axis off;colormap(hot);
subplot(2,4,3)
imagesc(pic(:,:,350));axis off;colormap(hot);
subplot(2,4,4)
imagesc(pic(:,:,450));axis off;colormap(hot);
subplot(2,4,[5,8])
imagesc(EvI);
axis off;colormap(hot);

