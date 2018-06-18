% File name: bpm_tbc.m
% Illustrates propagation of Gaussian pulse in a free space
% using BPM with transparent boundary conditions
% Operator P is determined in a separate function
close all;
clear all;

L_x=5.0; % transversal dimension (along x-axis)
w_0=1.0; % width of input Gaussian pulse
lambda = 0.6; % wavelength
beta0=12; %bg number

n=1.0; % refractive index of the medium
k_0=2*pi/lambda; % wavenumber
N_x=32; % number of points on x axis
N_y=N_x;
Delta_x=L_x/(N_x-1) % x axis spacing
Delta_y=Delta_x;
prefactor = 1/(2*n*k_0*Delta_x^2);
courant_number= 0.1;

h=5*(Delta_x+Delta_x); % propagation step

% h=courant_number*(Delta_x^2)
z_range1=10;
z_range2=10;
N_z=round(z_range1/h) % number of propagation steps
% N_z=2*round(z_range2/h); % number of propagation steps

% plotting=zeros(N_x,N_z); % storage for plotting
x=linspace(-0.5*L_x,0.5*L_x,N_x); % coordinates along x-axis
[X,Y]=meshgrid(x,x);
r=sqrt(X.^2 + Y.^2);
E=exp(-(r/(w_0)).^2); % initial Gaussian field
%E=LaguerreGaussianE([5,0,0,,r );


%  E=exp(-(r/(w_0)).^2).*besselj(0,beta0*r);
 %E=exp(-(r/(w_0)).^2).*sinc(X*12);
% E=ee;
%
amax=max(abs(E(:).^2));

z = 0;
z_plot = zeros(N_z);

plotting=zeros(N_x,N_y,N_z);
save('plotSplit.mat','plotting','-v7.3');
clear plotting;

m=matfile('plotSplit.mat','writable',true);
 for nn=1:N_z % BPM stepping
     nn
     sprintf('%d steps remaining',N_z-nn)
 z = z + h;
 z_plot(nn) = z + h;
 m.plotting(:,:,nn)=abs(E).^2;
 tic;
 for jj=1:N_y
    v=step(Delta_x,k_0,h,n,E(jj,:)');
    E(jj,:)=v';
 end
 
 for ii=1:N_x
     E(:,ii)=step(Delta_x,k_0,h,n,E(:,ii));
 end
 toc;
%  E=propTF(E,L_x,lambda,h);
 end;
 
 for ll = 1:(z/h)
 % choosing 2D plots every step
 if ll==1
         surf(X,Y,squeeze(m.plotting(:,:,ll)));
         shading interp;
         colormap jet;
         xlim([-5 5])
         ylim([-5 5])
         zlim([0 1])
 %         shading interp
 %         colorbar
         caxis manual
         caxis([0,amax])
         set(gca,'nextplot','replacechildren','visible','off')
         f = getframe(1);
         [im,map] = rgb2ind(f.cdata,256,'nodither');
         im(1,1,1,40) = 0;
     else
         surf(X,Y,squeeze(m.plotting(:,:,ll)));
         shading interp;
         colormap jet;
 xlim([-5 5])
         ylim([-5 5])
          shading interp
          colormap jet
 %         colorbar
         zlim([0 1])
         caxis manual
         caxis([0,amax])
         f = getframe(1);
         im(:,:,1,ll) = rgb2ind(f.cdata,map,'nodither');
     end
 % size of tick marks on both axes
 % pause(0.5)
 % pause
 end
 imwrite(im,map,'gaussian_split_64.gif','DelayTime',0.1,'LoopCount',inf)
 close all
 
 
 
 
 % %
 % for k = 1:N_z/10:N_z % choosing 2D plots every 10-th step
 % plot(plotting(:,k),'LineWidth',1.5)
 % set(gca,'FontSize',14); % size of tick marks on both axes
 % hold on
 % end
 % pause
 % close all
 % %
 for k = 1:ceil((z/h)/25):(z/h) % choosing 3D plots every 10-th step
 y = z_plot(k)*ones(size(x)); % spread out along y-axis
 plot3(x,y,squeeze(m.plotting(round(N_x/2),:,k)),'LineWidth',1.5)
 hold on
 end
 grid on
 xlabel('x','FontSize',14)
 ylabel('z','FontSize',14) % along propagation direction
 set(gca,'FontSize',14); % size of tick marks on both axes
 print('gauss_split_128.png','-dpng')
 % pause
 close all
 
 toc;