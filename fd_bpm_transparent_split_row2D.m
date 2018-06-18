function fd_bpm_transparent_split_row2D
% Propagation of Gaussian pulse in a free space by Crank-Nicholson method
% Transparent boundary conditions are introduced
clear all
L_x=10.0;
L_y=10.0;
% transversal dimension (along x-axis and y-axis)
w_0=1.0;
% width of input Gaussian pulse
lambda = 0.6;
% wavelength
n=1.0;
% refractive index of the medium
k_0=2*pi/lambda;
% wavenumber

N_x=25;
N_y=25;
% points on x and y axis
Delta_x=L_x/(N_x-1);
Delta_y=L_y/(N_y-1);
% x and y axis spacing

h=5*(Delta_x+Delta_y);
% propagation step along z-axis
N_z=80;
% number of propagation steps
plotting=zeros(N_x*N_y,N_z);
nplot=zeros(N_x,N_y,N_z);
% storage for plotting
x=linspace(-0.5*L_x,0.5*L_x,N_x);
y=linspace(-0.5*L_y,0.5*L_y,N_y);
% [X,Y]=meshgrid(x,y);
% E_mesh=exp(-(X.^2+Y.^2)/(w0^2));
% coordinates along x-axis and y-axis
x = x';
E=exp(-(x/w_0).^2)*exp(-(y/w_0).^2);
% initial Gaussian field

z = 0; z_plot = zeros(N_z); 
% E=transform(E);
for r=1:N_z
z = z + h;
z_plot(r) = z + h;
vE=transform(E);
plotting(:,r)=abs(vE).^2;
for j=1:N_y
    v= step(Delta_x,k_0,h,n,E(j,:)');
    E(j,:)=v';
end
for j=1:N_x
    E(:,j)= step(Delta_x,k_0,h,n,E(:,j));
end
end


%
% C=(1:N_x)'*(1:N_y);

for k = 1:N_z
% choosing 2D plots every step
nplot(:,:,k)=v2m(plotting(:,k),N_x,N_y);
if k==1
%         surf(1:N_x,1:N_y,nplot(:,:,k),C);
        surf(1:N_x,1:N_y,nplot(:,:,k));
        caxis manual
        caxis([0,1])
        colorbar
        xlim([0 25])
        ylim([0 25])
        zlim([0 1])
        set(gca,'nextplot','replacechildren','visible','off')
        f = getframe(1);
        [im,map] = rgb2ind(f.cdata,256,'nodither');
        im(1,1,1,40) = 0;
    else
        surf(1:N_x,1:N_y,nplot(:,:,k));
        caxis manual
        caxis([0,1])
        colorbar
        xlim([0 25])
        ylim([0 25])
        zlim([0 1])
        f = getframe(1);
        im(:,:,1,k) = rgb2ind(f.cdata,map,'nodither');
    end
% size of tick marks on both axes
% pause(0.5)
end
imwrite(im,map,'fd_bpm_transparent_split_row2D.gif','DelayTime',0.1,'LoopCount',inf)
close all


function [vE]=transform(mE);
% transform matrix E to vector E
[m,n]=size(mE);
vE=zeros(n*m,1);
for i=1:n
    vE(m*(i-1)+1:i*m)=mE(:,i);
end


function [mE]=v2m(vE,m,n);
% transform matrix E to vector E
mE=zeros(m,n);
for i=1:n
    mE(:,i)=vE(m*(i-1)+1:i*m);
end
