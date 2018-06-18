function fd_bpm_free2D
% Propagation of Gaussian pulse in a free space by Crank-Nicholson method
% No boundary conditions are introduced
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
% points on x axis
Delta_x=L_x/(N_x-1);
Delta_y=L_y/(N_y-1);
% x axis spacing
h=5*(Delta_x+Delta_y);
% propagation step along z-axis
N_z=60;
% number of propagation steps
plotting=zeros(N_x*N_y,N_z);
nplot=zeros(N_x,N_y,N_z);
% storage for plotting
x=linspace(-0.5*L_x,0.5*L_x,N_x);
y=linspace(-0.5*L_y,0.5*L_y,N_y);
% coordinates along x-axis
x = x';
E=exp(-(x/w_0).^2)*exp(-(y/w_0).^2);
% initial Gaussian field
%
% beta = n*k_0. With this choice, last term in propagator vanishes
prefactor = 1/(2*n*k_0*Delta_x^2); main = ones(N_x*N_y,1); 
above =ones(N_x*N_y-1,1); 
P = prefactor*(diag(above,-1)-2*diag(main,0)+diag(above,1)); % matrix P
below=ones(N_x*N_y-N_x,1);
Q = prefactor*(diag(below,-N_x)-2*diag(main,0)+diag(below,N_x));
step_plus = eye(N_x*N_y) + 0.25i*h*P;
% step forward
step_minus =eye(N_x*N_y)-0.25i*h*Q;
% step backward

%
z = 0; z_plot = zeros(N_z); 
E=transform(E);
for r=1:N_z
z = z + h;
z_plot(r) = z + h;
plotting(:,r)=abs(E).^2;
E=step_plus\step_minus*E;
E=v2v(E,N_x,N_y);
E=step_plus\step_minus*E;
E=iv2v(E,N_x,N_y);
end;
%

for k = 1:N_z
% choosing 2D plots every step
nplot(:,:,k)=v2m(plotting(:,k),N_x,N_y);
if k==1
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
imwrite(im,map,'DancingPeaks.gif','DelayTime',0.1,'LoopCount',inf)
close all
%
% for k = 1:N_z/2:N_z
% % choosing 3D plots every 10-th step
% y = z_plot(k)*ones(size(x));
% % spread out along y-axis
% plot3(x,y,plotting(:,k),'LineWidth',1.5)
% hold on
% end
% grid on 
% xlabel('x (mm)','FontSize',14)
% ylabel('z (mm)','FontSize',14)
% % along propagation direction
% set(gca,'FontSize',14);
% % size of tick marks on both axes
% pause(2)
% close all

function [vE]=transform(mE);
% transform matrix E to vector E
[m,n]=size(mE);
vE=zeros(n*m,1);
for i=1:n
    vE(m*(i-1)+1:i*m)=mE(:,i);
end

function [new_E]=v2v(E,m,n);
% transform E order
new_E=E;
for i=1:m
    for j=1:n
        new_E((i-1)*n+j)=E((j-1)*m+i);
    end
end

function [new_E]=iv2v(E,m,n);
% inverse of transform E order
new_E=E;
for i=1:m
    for j=1:n
        new_E((j-1)*m+i)=E((i-1)*n+j);
    end
end


function [mE]=v2m(vE,m,n);
% transform matrix E to vector E
mE=zeros(m,n);
for i=1:n
    mE(:,i)=vE(m*(i-1)+1:i*m);
end
