% File name: fd_bpm_free.m
% Propagation of Gaussian pulse in a free space by Crank-Nicholson method
% No boundary conditions are introduced
clear all
L_x=10.0; % transversal dimension (along x-axis)
w_0=1.0; % width of input Gaussian pulse
lambda = 0.6; % wavelength
n=1.0; % refractive index of the medium
k_0=2*pi/lambda; % wavenumber
N_x=25; % points on x axis
Delta_x=L_x/(N_x-1); % x axis spacing
h=5*Delta_x; % propagation step along z-axis
N_z=100; % number of propagation steps
plotting=zeros(N_x,N_z); % storage for plotting
x=linspace(-0.5*L_x,0.5*L_x,N_x); % coordinates along x-axis
x = x';
E=exp(-(x/w_0).^2); % initial Gaussian field
%
% beta = n*k_0. With this choice, last term in propagator vanishes
prefactor = 1/(2*n*k_0*Delta_x^2); main = ones(N_x,1); above =...
ones(N_x-1,1); below = above;
P = prefactor*(diag(above,-1)-2*diag(main,0)+diag(below,1)); % matrix P
%
step_plus = eye(N_x) + 0.5i*h*P; % step forward
step_minus =eye(N_x)-0.5i*h*P; % step backward
%
z = 0; z_plot = zeros(N_z); for r=1:N_z
    z = z + h;
z_plot(r) = z + h;
plotting(:,r)=abs(E).^2;
E=step_plus\step_minus*E;
end;
%
for k = 1:N_z/10:N_z % choosing 2D plots every 10-th step
plot(plotting(:,k),'LineWidth',1.5)
set(gca,'FontSize',14); % size of tick marks on both axes
hold on
end; pause; close all;
%
for k = 1:N_z/10:N_z % choosing 3D plots every 10-th step
y = z_plot(k)*ones(size(x)); % spread out along y-axis
plot3(x,y,plotting(:,k),'LineWidth',1.5)
hold on
end; grid on; xlabel('x (mm)','FontSize',14)
ylabel('z (mm)','FontSize',14) % along propagation direction
set(gca,'FontSize',14); % size of tick marks on both axes
pause; close all;