% File name: bpm_tbc.m
% Illustrates propagation of Gaussian pulse in a free space
% using BPM with transparent boundary conditions
% Operator P is determined in a separate function
close all;
clear all;
L_x=10.0; % transversal dimension (along x-axis)
w_0=1.0; % width of input Gaussian pulse
lambda = 0.6; % wavelength
beta0=12; %bg number

zR=pi*(w_0^2)/pi;



n=1.0; % refractive index of the medium
k_0=2*pi/lambda; % wavenumber
N_x=512; % number of points on x axis
N_y=N_x;
Delta_x=L_x/(N_x-1) % x axis spacing
Delta_y=Delta_x;
prefactor = 1/(2*n*k_0*Delta_x^2);
courant_number= 0.1;

h=0.5*(Delta_x+Delta_x); % propagation step





% h=courant_number*(Delta_x^2)
z_range1=17;
z_range2=1.5;
% N_z=round(z_range1/h); % number of propagation steps
N_z=2*round(z_range2/h); % number of propagation steps
N_z=round(N_z/9);
% plotting=zeros(N_x,N_z); % storage for plotting
x=linspace(-0.5*L_x,0.5*L_x,N_x); % coordinates along x-axis
[X,Y]=meshgrid(x,x);
r=sqrt(X.^2 + Y.^2);
% E=exp(-(r/(w_0)).^2); % initial Gaussian field



w=@(z)w_0*sqrt(1+(z/zR).^2); %beam waist
R=@(z)(z+((zR^2)./z)); %beam curvature
beta=@(z)beta0./(1+1i*(z/zR));
expon=@(z)exp(-((beta0*z/k_0)./w(z)).^2);

gau=@(r,z)exp(-(r/w(z)).^2);
bgau=@(r,z)(w_0./w(z))*gau(r,z).*expon(z).*besselj(0,beta(z)*r);%,N_x,N_y);

E=exp(-(r/(w_0)).^2).*besselj(0,beta0*r);
% E=ee;
%
% E=bgau(r,10);
% E(r>4.5)=0;
amax=max(abs(E(:).^2));

z = 0;
z_plot = zeros(N_z);

for nn=1:N_z % BPM stepping
    ARG{nn}=beta(z)*r;
    plo=bgau(r,z);plo(abs(plo)>1)=0;plo(r>2.5)=0;
%     plo=(w_0/w(z))*gau(r,z);
    plotting{nn}=abs(plo).^2;

z = z + h;
z_plot(nn) = z + h;

% for jj=1:N_y
%    v=step(Delta_x,k_0,h,n,E(jj,:)');
%    E(jj,:)=v';
% end
% 
% for ii=1:N_x
%     E(:,ii)=step(Delta_x,k_0,h,n,E(:,ii));
% end
end;
% 
for ll = 1:N_z
% choosing 2D plots every step
if ll==1
        surf(X,Y,plotting{ll});
        xlim([-5 5])
        ylim([-5 5])
        zlim([0 1])
        shading interp
%         colorbar
colormap jet;
        caxis manual
        caxis([0,1])
        set(gca,'nextplot','replacechildren','visible','off')
        f = getframe(1);
        [im,map] = rgb2ind(f.cdata,256,'nodither');
        im(1,1,1,40) = 0;
    else
        surf(X,Y,plotting{ll});
xlim([-5 5])
        ylim([-5 5])
        shading interp
        colormap jet;
%         colorbar
        zlim([0 1])
        caxis manual
        caxis([0,1])
        f = getframe(1);
        im(:,:,1,ll) = rgb2ind(f.cdata,map,'nodither');
    end
% size of tick marks on both axes
% pause(0.5)
% pause
end
imwrite(im,map,'analy.gif','DelayTime',0.001,'LoopCount',inf)
close all
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
for k = 1:ceil(N_z/25):N_z % choosing 3D plots every 10-th step
y = z_plot(k)*ones(size(x)); % spread out along y-axis
plot3(x,y,plotting{k}(round(N_x/2),:),'LineWidth',1.5)
hold on
end
grid on
xlabel('x','FontSize',14)
ylabel('z','FontSize',14) % along propagation direction
set(gca,'FontSize',14); % size of tick marks on both axes
print('fd_bpm_plot3.png','-dpng')
% % pause
% % close all






% arg1=ARG{1};arg2=ARG{2};
% arg3=ARG{3};arg4=ARG{4};
% arg5=ARG{5};arg6=ARG{6};
% arg7=ARG{7};arg8=ARG{8};
% arg9=ARG{9};arg10=ARG{10};
% arg11=ARG{11};arg12=ARG{12};
% arg13=ARG{13};arg14=ARG{14};
% arg15=ARG{15};arg16=ARG{16};
% arg17=ARG{17};arg18=ARG{18};
% cd 'C:\Users\Onur\Desktop\Temp';
% save('arg1.mat','arg1');save('arg2.mat','arg2');
% save('arg3.mat','arg3');save('arg4.mat','arg4');
% save('arg5.mat','arg5');save('arg6.mat','arg6');
% save('arg7.mat','arg7');save('arg8.mat','arg8');
% save('arg9.mat','arg9');save('arg10.mat','arg10');
% save('arg11.mat','arg11');save('arg12.mat','arg12');
% save('arg13.mat','arg13');save('arg14.mat','arg14');
% save('arg15.mat','arg15');save('arg16.mat','arg16');
% save('arg17.mat','arg17');save('arg18.mat','arg18');

% ahmet{1}=bess1;ahmet{2}=bess2;ahmet{3}=bess3;ahmet{4}=bess4;
% ahmet{5}=bess5;ahmet{6}=bess6; save('ahmet.mat','ahmet')
% load('c:\Users\Onur\Desktop\Temp\ahmet.mat');
% close all;
% 
% for kk=1:6;ahmet{kk}(ahmet{kk}>1)=0;%ahmet{kk}=normalize_this(ahmet{kk});end;
% end
% for ll = 1:6
% % choosing 2D plots every step
% if ll==1
%         surf(X,Y,ahmet{ll}.*plotting{ll});
%         xlim([-5 5])
%         ylim([-5 5])
%         zlim([0 0.4])
%         shading interp
% %         colorbar
% colormap jet;
%         caxis manual
%         caxis([0,1])
%         set(gca,'nextplot','replacechildren','visible','off')
%         f = getframe(1);
%         [im,map] = rgb2ind(f.cdata,256,'nodither');
%         im(1,1,1,40) = 0;
%     else
%         surf(X,Y,ahmet{ll}.*plotting{ll});
% xlim([-5 5])
%         ylim([-5 5])
%         shading interp
%         colormap jet;
% %         colorbar
%         zlim([0 0.4])
%         caxis manual
%         caxis([0,1])
%         f = getframe(1);
%         im(:,:,1,ll) = rgb2ind(f.cdata,map,'nodither');
%     end
% % size of tick marks on both axes
% % pause(0.5)
% pause
% end
% imwrite(im,map,'emek.gif','DelayTime',0.1,'LoopCount',inf)
