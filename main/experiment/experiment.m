% ========================================================================
% Introduction
% ========================================================================
% This code provides a simple demonstration of holograms generation of 
% steady optical beam (STOB) for spatial light modulator (SLM).
%
% Author: Yiqian Yang (yang-yq22@mails.tsinghua.edu.cn)
% =========================================================================

%% 
clear all;
clc;

%%   
Nx        = 3840;                                                          
Lx        = 3.74*1e-6*Nx;                                                  
Ny        = 2160;
Ly        = 3.74*1e-6*Ny;
x         = linspace(-Lx/2,Lx/2,Nx)';
y         = linspace(-Ly/2,Ly/2,Ny)';                 
[x,y]     = meshgrid(x,y);
r         = sqrt(x.^2+y.^2);
lambda    = 532*1e-9;                                                      
k0        = 2*pi/lambda;                                                    
 
%%
circ      = @(r,R)r<R;                                                     
cirNums   = 9;                                                             
R(1)      = Ly/sqrt(cirNums)/2;                                           
f(1)      = 0.3;                                                           
phase     = k0/(2*f(1))*(r.^2).*circ(r,R(1));                              
w(1)      = R(1)*f(1)/(sqrt((f(1))^2+(pi*R(1)^2/lambda)^2));               
RL(1)     = pi*w(1).^2/lambda;                                              

fprintf('RL')
disp(RL(1))

for i = 2:cirNums        
    
    f(i)  = f(i-1)+RL(i-1);
    R(i)  = sqrt(i)*R(1);   
    phase = phase + k0/(2*f(1))*(r.^2).*(circ(r,R(i))-circ(r,R(i-1)));   
    w(i)  = R(i)*f(i)/(sqrt((f(i))^2+(pi*R(i)^2/lambda)^2));                
    RL(i) = pi*w(i).^2/lambda; 
    fprintf('RL')
    disp(RL(i))
end

phase     = mod(phase,2*pi);
imagesc(phase);colormap(gray);axis off;
      
