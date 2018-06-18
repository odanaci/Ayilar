function [y,TC,CC]=DeltaKzA(X,Y,Ipu,a,w,nR,tol)

Ipu=normalize_this(Ipu);
Ipu=nthroot(Ipu,nR);
% Ipu=log(Ipu);
r=sqrt(X.^2+Y.^2);
RR=sqrt((X+a).^2+Y.^2);
TC=exp(-((RR-a)/w).^2); %tolerance cone
CC=ones(size(RR));%cutoff cone
CC(RR>a+1.5*w)=0;
CC(RR<a-1.5*w)=0;

y=TC.*Ipu;
% y(y<1e-2)=0;

 y=(1-y);
 y(y>1-tol)=1e3;
%  y(y==1)=1e2;
%  y(CC==0)=0;
%  y(CC==0)=1e5;

%  y(y<0.3)=0;
%  y=(1-y)*1e3;
%  y(CC==1)=1e5;
% y(TC==0)=1e5;


return;