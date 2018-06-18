function Anp1=RK47step(f,z,An,h)


k1=f(z,An);
k2=f(z + (h/2),An + (h/2)*k1);
k3=f(z+(h/2),An+(h/2)*k2);
k4=f(z+h,An+h*k3);
Anp1 = An + (h/6)*(k1 + 2*k2 + 2*k3 +k4);


% [m,n]= size(y);
% if m==1
%     y=y';
% end

return