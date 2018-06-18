function E_new = nstep(Delta_x,k_0,h,n,N_x,N_y,E_old)
% Function propagates BPM solution along one step
% determine size of the system
%--- Defines operator P outside of a boundary
prefactor = 1/(2*n*k_0*Delta_x^2); main = ones(N_x*N_y,1); 
above =ones(N_x*N_y-1,1); 
P = prefactor*sparse(diag(above,-1)-2*diag(main,0)+diag(above,1)); % matrix P
below=ones(N_x*N_y-N_x,1);
Q = prefactor*sparse(diag(below,-N_x)-2*diag(main,0)+diag(below,N_x));
clear main; clear above; clear below;
L_plus = speye(N_x*N_y) + 0.25i*h*P;
% step forward
L_minus =speye(N_x*N_y)-0.25i*h*Q;

clear P; clear Q;
% step backward
%
%---- Implementation of boundary conditions
%
pref = 0.25i*h/(2*k_0*Delta_x^2);
for j=1:N_y
k1=1i/Delta_x*log(E_old((j-1)*N_x+2)/E_old((j-1)*N_x+1));
k2=1i/Delta_x*log(E_old(N_x+j)/E_old(j));
if real(k1)<0
k1=1i*imag(k1);
end;
if real(k2)<0
    k2=1i*imag(k2);
end
left = pref*exp(1i*k1*Delta_x);
% left correction for next step
bottom=pref*exp(1i*k2*Delta_x);
% bottom correction for next step
L_plus((j-1)*N_x+1,(j-1)*N_x+1) = L_plus((j-1)*N_x+1,(j-1)*N_x+1)+left;
L_minus((j-1)*N_x+1,(j-1)*N_x+1) = L_minus((j-1)*N_x+1,(j-1)*N_x+1)-left;
L_plus(j,j)=L_plus(j,j)+bottom;
L_minus(j,j)=L_minus(j,j)-bottom;
%
k1=-1i/Delta_x*log(E_old(j*N_x)/E_old(j*N_x-1));
k2=-1i/Delta_x*log(E_old((N_y-1)*N_x+j)/E_old((N_y-2)*N_x+j));
if real(k1)<0
k1=1i*imag(k1);
end;
if real(k2)<0
    k2=1i*imag(k2);
end
right = pref*exp(1i*k1*Delta_x);
% right correction for next step
top=pref*exp(1i*k2*Delta_x);
% top correction for next step
L_plus(j*N_x,j*N_x) = L_plus(j*N_x,j*N_x) + right;
L_minus(j*N_x,j*N_x) = L_minus(j*N_x,j*N_x) - right;
L_plus((N_y-1)*N_x+j,(N_y-1)*N_x+j)=L_plus((N_y-1)*N_x+j,(N_y-1)*N_x+j)+top;
L_minus((N_y-1)*N_x+j,(N_y-1)*N_x+j)=L_minus((N_y-1)*N_x+j,(N_y-1)*N_x+j)-top;
% L_minus(N_x*(N_y-1)+j) = L_minus(N_x*(N_y-1)+j) - right;
end
%
% 'Lminus size is' 
% size(L_minus)
% 'Lplus size is'
% size(L_plus)
% 'Eold sizeis'
% size(E_old)

E_new = L_minus\L_plus*E_old;
% determine new solution
