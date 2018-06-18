% File name: bpm_tbc.m
% Illustrates propagation of Gaussian pulse in a free space
% using BPM with transparent boundary conditions
% Operator P is determined in a separate function
close all;
clear all; clc;
L_x=10.0; % transversal dimension (along x-axis)
w_0=1.0; % width of input Gaussian pulse

deL=1e-4;
kPu=8.05;
kPr=kPu+deL; kCo=kPu-deL;
% larmbdaPu = 0.78; % wavelength
beta0=12; %bg number

theta=1*pi/180;

%%%-20Mhz two-photon detuning
Xi_pp=0.05;
Xi_cc=0.01;
Xi_pc=1;%0.02;
Xi_cp=1;%0.02;
EPS=6.5e-6;




n=1.0; % refractive index of the medium
npu=n-EPS;
npr=sqrt(1+Xi_pp);
% k_0=2*pi/lambda; % wavenumber
N_x=512; % number of points on x axis
N_y=N_x;
Delta_x=L_x/(N_x-1) % x axis spacing
Delta_y=Delta_x;
% prefactor = 1/(2*n*k_0*Delta_x^2);
% courant_number= 0.1;

h=0.2*(Delta_x+Delta_x); % propagation step

% h=courant_number*(Delta_x^2)
z_range1=0.2;
z_range2=1.5;
N_z=round(z_range1/h) % number of propagation steps
 %N_z=2*round(z_range2/h); % number of propagation steps

% plotting=zeros(N_x,N_z); % storage for plotting
x=linspace(-0.5*L_x,0.5*L_x,N_x); % coordinates along x-axis
[X,Y]=meshgrid(x,x);
r=sqrt(X.^2 + Y.^2);
%E=exp(-(r/(w_0)).^2); % initial Gaussian field

 Epr=exp(-(r/(0.2*w_0)).^2);%.*besselj(0,beta0*r);
 Econj=eps.*ones(size(X));%0.1.*Epr;
% E=tilt(E,L_x,lambda,0,theta);% E=ee;
%
% Ipu=abs(exp(-(r/(w_0)).^2).*sinc(Y*5));
E=rect(Y*2).*exp(-(r/(0.4*w_0)).^2);
% E=propTF(E,L_x,0.6,0.01);
E=propTF(E,L_x,0.6,0.01+4*0.2035);

pause;
Ipu=(abs(E).^2);
% phiPr=normalize_this(-(r.^2)/2); %phiPr=phiGa-min(phiGa(:));
% [kPrx,kPry,kPrz]=surfnorm(X,Y,phiPr);
% 
% kPrz=kPrz*cos(theta)-kPry*sin(theta);
% phiPu=phiPr; phiPu(Ipu<0.05)=0;[kPux,kPuy,kPuz]=surfnorm(X,Y,phiPu);
% kPux(phiPu==0)=0;kPuy(phiPu==0)=0;kPuz(phiPu==0)=0;
theta=-1*pi/180;

% deltaKz=DeltaKz(X,Y,theta,20,Ipu,kPu,kPr,kCo,npu,npr);
deltaKz=DeltaKz(X,Y,theta,20,Ipu,kPu,kPr,kCo,npu,npr);

surfMe(deltaKz);colorbar;

deltaKz_=deltaKz;deltaKz_(deltaKz>1)=1e5;deltaKz_(deltaKz<=1)=0;
% pause;
% IdeltaKz=deltaKz;IdeltaKz(real(exp(1i*deltaKz))<0.9)=0;
% pause;
%abs(2*kPu*npu*kPuz-kPr*npr*kPrz-kCo*cos(theta));
% deltaKz=abs(exp(-(r/(w_0)).^2).*sinc(Y*5));
% pause;
clear Ipu;

amax=max(abs(Epr(:).^2));
% f=@(z,A)10*[i*kPr*Xi_pc*exp(1i*deltaKz*z).*conj(A(end/2+1:end,:))/2;...
%     i*kCo*Xi_cp*exp(1i*deltaKz*z).*conj(A(1:end/2,:))/2];
% f=@(z,A)3i*[kPr*Xi_pc*exp(-1*deltaKz_*z).*conj(A(end/2+1:end,:))/2;...
%     kCo*Xi_cp*exp(-1*deltaKz_*z).*conj(A(1:end/2,:))/2];


z = 0; z_plot = zeros(N_z); 
% plotting{1}=abs(E).^2;
% E=E(:);

% E=E(:);
% 'line 44 E size is'
% size(E)

plotting=zeros(N_x,N_y,N_z);
plottingConj=plotting;
save('plotFFT.mat','plotting','plottingConj','-v7.3');
clear plotting;clear plottingConj;

m=matfile('plotFFT.mat','writable',true);



for nn=1:N_z
    nn
    sprintf('%d more steps',N_z-nn)
z = z + h;
% pause;


z_plot(nn) = z + h;
m.plotting(:,:,nn)=abs(Epr).^2;
m.plottingConj(:,:,nn)=abs(Econj).^2;
% 'line 50 E size'
% size(E)
tic;
% E = nstep(Delta_x,k_0,h,n,N_x,N_y,E);
E=propTF(E,L_x,0.6,4*h);
% Ipu=abs(E).^2;
% clear Ipu;
f=@(z,A)10*[i*kPr*Xi_pc*exp(1i*deltaKz*z).*conj(A(end/2+1:end,:))/2;...
    i*kCo*Xi_cp*exp(1i*deltaKz*z).*conj(A(1:end/2,:))/2];


Epr=fftBpmDisp(Epr,L_x,2*pi/(kPr),Xi_pp,h/2);
Econj=fftBpmDisp(Epr,L_x,2*pi/(kCo),Xi_cc,h/2);
[Epr,Econj]=splitThat(RK47step(f,z,[Epr;Econj],h));

Epr=fftBpmDisp(Epr,L_x,2*pi/(kPr),Xi_pp,h/2);
Econj=fftBpmDisp(Epr,L_x,2*pi/(kCo),Xi_cc,h/2);
deltaKz=DeltaKz(X,Y,theta,20,abs(E).^2,kPu,kPr,kCo,npu,npr);


% 'line 54 E size'
% size(E)
% E=v2v(E,N_x,N_y);

% E=v2vonur(E,N_x,N_y);
% 'line 57 E size'
% size(E)

% % E=reshape(E',1,N_x*N_y);
% E = nstep(Delta_y,k_0,h,n,N_y,N_x,E);
% 'line 62 E size'
% size(E)
% % E=reshape(reshape(E,N_x,N_y)',1,N_x*N_y);
% E=iv2v(E,N_x,N_y);

% E=iv2vonur(E,N_x,N_y);
toc;
% E=iv2vonur(E,N_x,N_y);
% %'line 66 E size'
% size(E)
end;
toc;

save('plotFFT.mat','X','Y','z','h','-append');
% clear all;
% for ll = 1:(z/h)%N_z;%ceil(N_z/20):N_z
% % choosing 2D plots every step
% if ll==1
%         surf(X,Y,squeeze(m.plotting(:,:,ll)));
%         xlim([-5 5])
%         ylim([-5 5])
%         zlim([0 1])
%          shading interp
%          colormap jet
%         caxis manual
%         caxis([0,amax])
%         set(gca,'nextplot','replacechildren','visible','off')
%         f = getframe(1);
%         [im,map] = rgb2ind(f.cdata,256,'nodither');
%         im(1,1,1,40) = 0;
%     else
%         surf(X,Y,squeeze(m.plotting(:,:,ll)));
% xlim([-5 5])
%         ylim([-5 5])
%          shading interp
%          colormap jet
% %         colorbar
%         zlim([0 1])
%         caxis manual
%         caxis([0,amax])
%         f = getframe(1);
%         im(:,:,1,ll) = rgb2ind(f.cdata,map,'nodither');
%     end
% % size of tick marks on both axes
% % pause(0.5)
% % pause
% end
% imwrite(im,map,'gaussian_nosplit_35.gif','DelayTime',0.1,'LoopCount',inf)
% close all
% 
% 
% 
% 
% % %
% % for k = 1:N_z/10:N_z % choosing 2D plots every 10-th step
% % plot(plotting(:,k),'LineWidth',1.5)
% % set(gca,'FontSize',14); % size of tick marks on both axes
% % hold on
% % end
% % pause
% % close all
% % %
% % for k = 1:round(N_z/20):N_z % choosing 3D plots every 10-th step
% for k = 1:ceil((z/h)/20):(z/h) % choosing 3D plots every 10-th step
% 
% y = z_plot(k)*ones(size(x)); % spread out along y-axis
% plot3(x,y,squeeze(m.plotting(round(N_x/2),:,k)),'LineWidth',1.5)
% hold on
% end
% grid on
% xlabel('x','FontSize',14)
% ylabel('z','FontSize',14) % along propagation direction
% set(gca,'FontSize',14); % size of tick marks on both axes
% print('gaussian_nosplit_35_plot3.png','-dpng')
% % pause
% % close all
% 
% 
