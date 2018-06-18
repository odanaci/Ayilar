function y=TCfilter(E,X,Y,a,w)

r=sqrt(X.^2+Y.^2);
RR=sqrt((X+a).^2+Y.^2);
TC=exp(-((RR-a)/w).^2); %tolerance cone

y=E.*TC;

return;