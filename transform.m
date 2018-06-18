function [vE]=transform(mE);
% transform matrix E to vector E
[m,n]=size(mE);
vE=zeros(n*m,1);
for i=1:n
    vE(m*(i-1)+1:i*m)=mE(:,i);
end

