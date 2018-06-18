function y=DeltaKz(X,Y,theta,R,Ipu,kPu,kPr,kCo,npu,npr)

r=sqrt(X.^2+Y.^2);
phiPr=normalize_this(-(r.^2)/R); %phiPr=phiPr-min(phiPr(:));
[kPrx,kPry,kPrz]=surfnorm(X,Y,phiPr);
% kPrz=phiPr;
kPrz=kPrz*cos(theta)-kPrx*sin(theta);
phiPu=phiPr; 
% phiPu(Ipu<0.05)=0;
phiPu(Ipu<0.01)=0;

% phiPu=phiPu.*Ipu;
[kPux,kPuy,kPuz]=surfnorm(X,Y,phiPu);
kPux(phiPu==0)=0;kPuy(phiPu==0)=0;kPuz(phiPu==0)=0;

y=abs(2*kPu*npu*kPuz-kPr*npr*kPrz-kCo*kPrz); %cos(theta));

return;