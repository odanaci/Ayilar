function y=RK47(f,t,w0)

y(:,1)=w0;
h=abs(t(2)-t(1));
N=(max(t)-min(t))/h;

for jj=1:N
k1=f(t(jj),y(:,jj));
k2=f(t(jj) + (h/2),y(:,jj) + (h/2)*k1);
k3=f(t(jj)+(h/2),y(:,jj)+(h/2)*k2);
k4=f(t(jj)+h,y(:,jj)+h*k3);
y(:,jj +1 ) = y(:,jj) + (h/6)*(k1 + 2*k2 + 2*k3 +k4);
end

[m,n]= size(y);
if m==1
    y=y';
end

return