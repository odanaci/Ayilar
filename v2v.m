function [new_E]=v2v(E,m,n);
% transform E order
new_E=E;
for i=1:m
    for j=1:n
        new_E((i-1)*n+j)=E((j-1)*m+i);
    end
end

