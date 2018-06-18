% File name: bpm_tbc.m
% Illustrates propagation of Gaussian pulse in a free space
% using BPM with transparent boundary conditions
% Operator P is determined in a separate function
clear all;
L_x=10.0; % transversal dimension (along x-axis)
w_0=1.0; % width of input Gaussian pulse
lambda = 0.6; % wavelength
n=1.0; % refractive index of the medium
k_0=2*pi/lambda; % wavenumber
N_x=128; % number of points on x axis
Delta_x=L_x/(N_x-1); % x axis spacing
h=5*Delta_x; % propagation step
N_z=100; % number of propagation steps
plotting=zeros(N_x,N_z); % storage for plotting
x=linspace(-0.5*L_x,0.5*L_x,N_x); % coordinates along x-axis
x = x';
E=exp(-(x/w_0).^2); % initial Gaussian field
%
z = 0;
z_plot = zeros(N_z);
for r=1:N_z % BPM stepping
z = z + h;
z_plot(r) = z + h;
plotting(:,r)=abs(E).^2;
E = step(Delta_x,k_0,h,n,E); % Propagates pulse over one step
end;
%
for k = 1:N_z/10:N_z % choosing 2D plots every 10-th step
plot(plotting(:,k),'LineWidth',1.5)
set(gca,'FontSize',14); % size of tick marks on both axes
hold on
end
pause
close all
%
for k = 1:N_z/10:N_z % choosing 3D plots every 10-th step
y = z_plot(k)*ones(size(x)); % spread out along y-axis
plot3(x,y,plotting(:,k),'LineWidth',1.5)
hold on
end
grid on
xlabel('x','FontSize',14)
ylabel('z','FontSize',14) % along propagation direction
set(gca,'FontSize',14); % size of tick marks on both axes
pause
% close all