% File name: step.m
function E_new = step(Delta_x,k_0,h,n,E_old)
% Function propagates BPM solution along one step
N_x =size(E_old,1); % determine size of the system
%--- Defines operator P outside of a boundary
prefactor = 1/(2*n*k_0*Delta_x^2);
main = ones(N_x,1);
above = ones(N_x-1,1);
below = above;
P = prefactor*(diag(above,-1)-2*diag(main,0)+diag(below,1)); % matrix P
%
L_plus = sparse(eye(N_x) + 0.5i*h*P); % step forward
L_minus = sparse(eye(N_x)-0.5i*h*P); % step backward
%
%---- Implementation of boundary conditions
%
pref = 0.5i*h/(2*k_0*Delta_x^2);
k=1i/Delta_x*log(E_old(2)/E_old(1));
if real(k)<0
k=1i*imag(k);
end;
left = pref*exp(1i*k*Delta_x); % left correction for next step
L_plus(1) = L_plus(1)+left;
L_minus(1) = L_minus(1)-left;
%
k=-1i/Delta_x*log(E_old(N_x)/E_old(N_x-1));
if real(k)<0
k=1i*imag(k);
end;
right = pref*exp(1i*k*Delta_x); % right correction for nest step
L_plus(N_x) = L_plus(N_x) + right;
L_minus(N_x) = L_minus(N_x) - right;
%
E_new = L_minus\L_plus*E_old; % determine new solution
