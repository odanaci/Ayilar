function [new_E]=iv2v(E,m,n);
% inverse of transform E order
new_E=E;
for i=1:m
    for j=1:n
        new_E((j-1)*m+i)=E((i-1)*n+j);
    end
end

