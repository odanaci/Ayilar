% function fd_bpm_transparent2D
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

N_x=64;
N_y=N_x;
% points on x and y axis
Delta_x=L_x/(N_x-1);
Delta_y=L_y/(N_y-1);
% x and y axis spacing

h=5*(Delta_x+Delta_y);
% propagation step along z-axis
N_z=80;
z_range1=10;
z_range2=1.5;
N_z=round(z_range1/h); % number of propagation steps
% N_z=2*round(z_range2/Delta_x); % number of propagation steps
% number of propagation steps
plotting=zeros(N_x*N_y,N_z);
nplot=zeros(N_x,N_y,N_z);
% storage for plotting
x=linspace(-0.5*L_x,0.5*L_x,N_x);
y=linspace(-0.5*L_y,0.5*L_y,N_y);
% coordinates along x-axis and y-axis
x = x';

[X,Y]=meshgrid(x,y);
r=sqrt(X.^2 + Y.^2);

%E=exp(-(r/w_0).^2);
 E=exp(-(r/w_0).^2).*besselj(0,12*r);

% E=exp(-(x/w_0).^2)*exp(-(y/w_0).^2);

amax=max(abs(E(:)).^2);
% initial Gaussian field
%


z = 0; z_plot = zeros(N_z); 
E=transform(E);
for r=1:N_z
z = z + h;
z_plot(r) = z + h;
plotting(:,r)=abs(E).^2;
E = nstep(Delta_x,k_0,h,n,N_x,N_y,E);
E=v2v(E,N_x,N_y);
E = nstep(Delta_y,k_0,h,n,N_y,N_x,E);
E=iv2v(E,N_x,N_y);
end;
%

for k = 1:N_z
% choosing 2D plots every step
nplot(:,:,k)=v2m(plotting(:,k),N_x,N_y);
if k==1
        surf(1:N_x,1:N_y,nplot(:,:,k));
        xlim([0 25])
        ylim([0 25])
        zlim([0 1])
        shading interp;
        colormap jet;
        caxis manual
        caxis([0,amax])
        set(gca,'nextplot','replacechildren','visible','off')
        f = getframe(1);
        [im,map] = rgb2ind(f.cdata,256,'nodither');
        im(1,1,1,40) = 0;
    else
        surf(1:N_x,1:N_y,nplot(:,:,k));
        xlim([0 25])
        ylim([0 25])
        zlim([0 1])
        shading interp;
        colormap jet;
        caxis manual
        caxis([0,amax])
        f = getframe(1);
        im(:,:,1,k) = rgb2ind(f.cdata,map,'nodither');
    end
% size of tick marks on both axes
% pause(0.5)
end
imwrite(im,map,'gauss_trans_25.gif','DelayTime',0.1,'LoopCount',inf)
close all

for k = 1:ceil((z/h)/25):(z/h) % choosing 3D plots every 10-th step
 y = z_plot(k)*ones(size(x)); % spread out along y-axis
 plot3(x,y,squeeze(nplot(round(N_x/2),:,k)),'LineWidth',1.5)
 hold on
 end
 grid on
 xlabel('x','FontSize',14)
 ylabel('z','FontSize',14) % along propagation direction
 set(gca,'FontSize',14); % size of tick marks on both axes
 print('gauss_split_25.png','-dpng')
 % pause
 close all
 


% 
% function [vE]=transform(mE);
% % transform matrix E to vector E
% [m,n]=size(mE);
% vE=zeros(n*m,1);
% for i=1:n
%     vE(m*(i-1)+1:i*m)=mE(:,i);
% end
% 
% function [new_E]=v2v(E,m,n);
% % transform E order
% new_E=E;
% for i=1:m
%     for j=1:n
%         new_E((i-1)*n+j)=E((j-1)*m+i);
%     end
% end
% 
% function [new_E]=iv2v(E,m,n);
% % inverse of transform E order
% new_E=E;
% for i=1:m
%     for j=1:n
%         new_E((j-1)*m+i)=E((i-1)*n+j);
%     end
% end
% 
% 
% function [mE]=v2m(vE,m,n);
% % transform matrix E to vector E
% mE=zeros(m,n);
% for i=1:n
%     mE(:,i)=vE(m*(i-1)+1:i*m);
% end
% 
% function E_new = nstep(Delta_x,k_0,h,n,N_x,N_y,E_old)
% % Function propagates BPM solution along one step
% % determine size of the system
% %--- Defines operator P outside of a boundary
% prefactor = 1/(2*n*k_0*Delta_x^2); main = ones(N_x*N_y,1); 
% above =ones(N_x*N_y-1,1); 
% P = prefactor*(diag(above,-1)-2*diag(main,0)+diag(above,1)); % matrix P
% below=ones(N_x*N_y-N_x,1);
% Q = prefactor*(diag(below,-N_x)-2*diag(main,0)+diag(below,N_x));
% L_plus = eye(N_x*N_y) + 0.25i*h*P;
% % step forward
% L_minus =eye(N_x*N_y)-0.25i*h*Q;
% % step backward
% %
% %---- Implementation of boundary conditions
% %
% pref = 0.25i*h/(2*k_0*Delta_x^2);
% for j=1:N_y
% k1=1i/Delta_x*log(E_old((j-1)*N_x+2)/E_old((j-1)*N_x+1));
% k2=1i/Delta_x*log(E_old(N_x+j)/E_old(j));
% if real(k1)<0
% k1=1i*imag(k1);
% end;
% if real(k2)<0
%     k2=1i*imag(k2);
% end
% left = pref*exp(1i*k1*Delta_x);
% % left correction for next step
% bottom=pref*exp(1i*k2*Delta_x);
% % bottom correction for next step
% L_plus((j-1)*N_x+1,(j-1)*N_x+1) = L_plus((j-1)*N_x+1,(j-1)*N_x+1)+left;
% L_minus((j-1)*N_x+1,(j-1)*N_x+1) = L_minus((j-1)*N_x+1,(j-1)*N_x+1)-left;
% L_plus(j,j)=L_plus(j,j)+bottom;
% L_minus(j,j)=L_minus(j,j)-bottom;
% %
% k1=-1i/Delta_x*log(E_old(j*N_x)/E_old(j*N_x-1));
% k2=-1i/Delta_x*log(E_old((N_y-1)*N_x+j)/E_old((N_y-2)*N_x+j));
% if real(k1)<0
% k1=1i*imag(k1);
% end;
% if real(k2)<0
%     k2=1i*imag(k2);
% end
% right = pref*exp(1i*k1*Delta_x);
% % right correction for next step
% top=pref*exp(1i*k2*Delta_x);
% % top correction for next step
% L_plus(j*N_x,j*N_x) = L_plus(j*N_x,j*N_x) + right;
% L_minus(j*N_x,j*N_x) = L_minus(j*N_x,j*N_x) - right;
% L_plus((N_y-1)*N_x+j,(N_y-1)*N_x+j)=L_plus((N_y-1)*N_x+j,(N_y-1)*N_x+j)+top;
% L_minus((N_y-1)*N_x+j,(N_y-1)*N_x+j)=L_minus((N_y-1)*N_x+j,(N_y-1)*N_x+j)-top;
% % L_minus(N_x*(N_y-1)+j) = L_minus(N_x*(N_y-1)+j) - right;
% end
% %
% E_new = L_minus\L_plus*E_old;
% % determine new solution
