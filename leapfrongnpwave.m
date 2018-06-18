function y=leapfrognpwave(x,t,c,func,prev)
N=length(t);
T=max(t); k=abs(t(2)-t(1));
M=length(x);h=abs(x(2)-x(1));
F=func(x);
% An=analy(x,-k);
% Sltn_store=sparse(N,M);
alpha=c*k/(h);

% D=sparse(1:M,1:M,(1-alpha)*ones(1,M),M,M);E=sparse(1:M-1,2:M,alpha*ones(1,M-1),M,M);

D=sparse(1:M,1:M,ones(1,M),M,M);E=sparse(1:M-1,2:M,ones(1,M-1),M,M);
S=sparse(M,M);

D=2*(1-alpha^2)*D;
% uw=E-D; %uw(M,1)=1;
lf=(alpha^2)*(E+E'); 
lf(M,M-1)=0; 

A=D+lf;

% A(1,1)=(1-alpha/2);
% A(1,2)=alpha/2;
% A(M,M)=(1-alpha/2);
% A(M,M-1)=0;
% 
% 'onur cond number is '
% 
% cond(full(A))
% 
% 'onur norm is '
% norm(full(A))


% 
% above=ones(M-1,1);
% Alin=diag(above,1)+diag(-above,-1);
% Alin=alpha*Alin;
% Alin(1,1)=1-alpha/2; Alin(1,2)=alpha/2;
% Alin(M,M-1)=0; Alin(M,M)=1-alpha/2;

% 'lin cond number is '
% cond(Alin)
% 
% 'lin norm is '
% norm(Alin)
% 
% 'norm diff'
% norm(Alin-full(A))
% 
% 'max diff'
% max(max(Alin-full(A)))
% 
% imwrite(Alin-full(A),'matrix.png')
% A=D+E;


% A(M,1)=alpha;
% 
% sprintf('Condition Number of A=%.d',cond(A))
% sprintf('2 norm of A=%.d',normest(A))



Un=F';

[mm,nn]=size(prev);
if mm<nn
prev=prev';
end
prev(1)=0;
prev(end)=0;
% Sltn_store[1,:]=Un;
% 
% 'size of prev is '
% size(prev)
% 'size of A is '
% size(A)
% 'size of Un is '
% size(Un)
plotting{1}=real(Un);%.^2;

for n=2:N
    Unp1=A*Un - prev;
    prev=Un;
   prev(1)=0;
 prev(end)=0;
    Un=Unp1;
    plotting{n}=real(Un);%abs(Un).^2;

%     Sltn_store[n,:]=Un;
end

y=plotting;