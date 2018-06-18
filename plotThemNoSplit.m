% ft_bpm2D;
clear all;
close all;
% m=matfile('plotNoSplitC.mat','writable',true);
% m=matfile('plotNoSplit.mat','writable',true);
m=matfile('plotFFT2nd.mat','writable',true);



z=m.z; h=m.h;X=m.X;Y=m.Y;amax=max(max(m.plotting(:,:,1)));
z_plot=0:h:z;x=X(1,:);y=Y(:,1);N_x=size(x); N_y=N_x;
% 
% for ll = 1:(z/h)%N_z;%ceil(N_z/20):N_z
% % choosing 2D plots every step
% ll
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
%                 caxis([0,max(max(max(squeeze(m.plotting))))])
%         f = getframe(1);
%         im(:,:,1,ll) = rgb2ind(f.cdata,map,'nodither');
%     end
% % size of tick marks on both axes
% % pause(0.5)
% % pause
% end
% imwrite(im,map,'probe.gif','DelayTime',0.1,'LoopCount',inf)
% close all

for ll = 1:(z/h)%N_z;%ceil(N_z/20):N_z
% choosing 2D plots every step
ll
if ll==1
        surf(zoom_(X,0.3,0.7),zoom_(Y,0.3,0.7),zoom_(TCfilter(squeeze(m.plotting(:,:,ll)),X,Y,3,0.15),0.3,0.7));
%         xlim([-5 5])
%         ylim([-5 5])
%         zlim([0 1])
        view(2);
         shading interp
         colormap jet
        caxis manual
                caxis([0,max(max(max(squeeze(m.plotting))))])

%         caxis([0,amax])
        axis equal tight;
        set(gca,'nextplot','replacechildren','visible','off')
        f = getframe(1);
        im = rgb2ind(f.cdata,jet,'nodither');
%         im(1,1,1,40) = 0;
    else
        surf(zoom_(X,0.3,0.7),zoom_(Y,0.3,0.7),zoom_(TCfilter(squeeze(m.plotting(:,:,ll)),X,Y,3,0.15),0.3,0.7));
% xlim([-5 5])
%         ylim([-5 5])
         shading interp
         colormap jet
%         colorbar
%         zlim([0 1])
        view(2);
        caxis manual
                        caxis([0,max(max(max(squeeze(m.plotting))))])

%         caxis([0,amax])
        axis equal tight;
        f = getframe(gcf);
        im(:,:,1,ll) = rgb2ind(f.cdata,jet,'nodither');
    end
% size of tick marks on both axes
% pause(0.5)
% pause
end
imwrite(im,jet,'probe_2d.gif','DelayTime',0.1,'LoopCount',inf)
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
% for k = 1:round(N_z/20):N_z % choosing 3D plots every 10-th step
% for k = 1:ceil((z/h)/20):(z/h) % choosing 3D plots every 10-th step
% 
% y = z_plot(k)*ones(size(X(:,1))); % spread out along y-axis
% plot3(x,y,squeeze(m.plotting(round(N_x/2),:,k)),'LineWidth',1.5)
% hold on
% end
% grid on
% xlabel('x','FontSize',14)
% ylabel('z','FontSize',14) % along propagation direction
% set(gca,'FontSize',14); % size of tick marks on both axes
% print('probe_plot3X.png','-dpng')
% % pause
% close all
% 
% for k = 1:ceil((z/h)/20):(z/h) % choosing 3D plots every 10-th step
% 
% y = z_plot(k)*ones(size(X(:,1))); % spread out along y-axis
% plot3(x,y,squeeze(m.plotting(:,round(N_x/2),k)),'LineWidth',1.5)
% hold on
% end
% grid on
% xlabel('x','FontSize',14)
% ylabel('z','FontSize',14) % along propagation direction
% set(gca,'FontSize',14); % size of tick marks on both axes
% print('probe_plot3Y.png','-dpng')
% close all;
% 
% 
% for ll = 1:(z/h)%N_z;%ceil(N_z/20):N_z
% % choosing 2D plots every step
% ll
% if ll==1
%         surf(X,Y,squeeze(m.plottingConj(:,:,ll)));
% %         surf(zoom_(X,0.3,0.7),zoom_(Y,0.3,0.7),zoom_(squeeze(m.plottingConj(:,:,ll)),0.3,0.7))
%         xlim([-5 5])
%         ylim([-5 5])
%         zlim([0 1])
% %         title(sprintf('%d',ll))
%          shading interp
%          colormap jet
% %         caxis manual
% %         caxis([0,max(max(max(squeeze(m.plottingConj))))])
%         set(gca,'nextplot','replacechildren','visible','off')
%         f = getframe(1);
%         [im,map] = rgb2ind(f.cdata,256,'nodither');
%         im(1,1,1,40) = 0;
%     else
%         surf(X,Y,squeeze(m.plottingConj(:,:,ll)));
% %         surf(zoom_(X,0.3,0.7),zoom_(Y,0.3,0.7),zoom_(squeeze(m.plottingConj(:,:,ll)),0.3,0.7))
% 
% xlim([-5 5])
%         ylim([-5 5])
% %         title(sprintf('%d',ll))
%          shading interp
%          colormap jet
% %         colorbar
%         zlim([0 1])
% %         caxis manual
% %         caxis([0,max(max(max(squeeze(m.plottingConj))))])
%         f = getframe(1);
%         im(:,:,1,ll) = rgb2ind(f.cdata,map,'nodither');
%     end
% % size of tick marks on both axes
% % pause(0.5)
% % pause
% end
% imwrite(im,map,'conj.gif','DelayTime',0.1,'LoopCount',inf)
% close all

for ll = 1:(z/h)%N_z;%ceil(N_z/20):N_z
% choosing 2D plots every step
ll
if ll==1
        surf(zoom_(X,0.3,0.7),zoom_(Y,0.3,0.7),zoom_(fliplr(TCfilter(squeeze(m.plottingConj(:,:,ll)),X,Y,3,0.15)),0.3,0.7));
%         xlim([-5 5])
%         ylim([-5 5])
%         zlim([0 1])
        view(2);
%         title(sprintf('%d',ll))
         shading interp
         colormap jet
        caxis manual
        caxis([0,max(max(max(squeeze(m.plottingConj))))])
        axis equal tight;
        set(gca,'nextplot','replacechildren','visible','off')
        f = getframe(1);
%         [im,map] = rgb2ind(f.cdata,256,'nodither');
                [im,map] = rgb2ind(f.cdata,256,'nodither');
%  imwrite(imind,cm,outfile,'gif','DelayTime',0,'loopcount',inf);

imwrite(im,jet,'conj_2d.gif','gif','DelayTime',0.1,'LoopCount',inf)
%         im(1,1,1,ll) = 0;
    else
        surf(zoom_(X,0.3,0.7),zoom_(Y,0.3,0.7),zoom_(fliplr(TCfilter(squeeze(m.plottingConj(:,:,ll)),X,Y,3,0.15)),0.3,0.7));
% xlim([-5 5])
%         ylim([-5 5])
         shading interp
         colormap jet
%         colorbar
%         zlim([0 1])
        view(2);
        caxis manual
                caxis([0,max(max(max(squeeze(m.plottingConj))))])

%         caxis([0,amax])
        axis equal tight;
        f = getframe(1);
%         im(:,:,1,ll) = rgb2ind(f.cdata,map,'nodither');
%                 im(:,:,1,ll) = rgb2ind(frame2im(f),map,'nodither');
                 im= rgb2ind(f.cdata,jet,'nodither');

    imwrite(im,jet,'conj_2d.gif','gif','DelayTime',0.1,'writemode','append');

                
    end
% size of tick marks on both axes
% pause(0.5)
% pause
end
% imwrite(imind,cm,outfile,'gif','DelayTime',0,'writemode','append');
% imwrite(im,map,'conj_2d.gif','gif','DelayTime',0.1,'writemode','append');

% imwrite(im,map,'conj_2d.gif','DelayTime',0.5,'LoopCount',inf)
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
% for k = 1:round(N_z/20):N_z % choosing 3D plots every 10-th step
for k = 1:ceil((z/h)/20):(z/h) % choosing 3D plots every 10-th step

y = z_plot(k)*ones(size(X(:,1))); % spread out along y-axis
plot3(x,y,squeeze(m.plottingConj(round(N_x/2),:,k)),'LineWidth',1.5)
hold on
end
grid on
xlabel('x','FontSize',14)
ylabel('z','FontSize',14) % along propagation direction
set(gca,'FontSize',14); % size of tick marks on both axes
print('conj_plot3X.png','-dpng')
% pause
close all

for k = 1:ceil((z/h)/20):(z/h) % choosing 3D plots every 10-th step

y = z_plot(k)*ones(size(X(:,1))); % spread out along y-axis
plot3(x,y,squeeze(m.plottingConj(:,round(N_x/2),k)),'LineWidth',1.5)
hold on
end
grid on
xlabel('x','FontSize',14)
ylabel('z','FontSize',14) % along propagation direction
set(gca,'FontSize',14); % size of tick marks on both axes
print('conj_plot3Y.png','-dpng')
close all;

BB=zoom_(TCfilter(squeeze(m.plotting(:,:,ll)),X,Y,3,0.2),0.3,0.7);
AA=zoom_(fliplr(TCfilter(squeeze(m.plottingConj(:,:,ll)),X,Y,3,0.2)),0.3,0.7);
A=[AA,BB];
cScale=[0,max(max(max(squeeze(m.plotting))))];
figure;surf(A);shading interp;colormap jet;view(2);caxis manual;caxis(cScale);axis equal tight;
set(gca,'xtick',[])
set(gca,'ytick',[])
% print(gcf,'2nd','-depsc','-tiff')
% print(gcf,'1st','-dpng')
print(gcf,'2nd','-dpng')
