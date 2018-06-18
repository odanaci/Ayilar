
function [mE]=v2m(vE,m,n);
% transform matrix E to vector E
mE=zeros(m,n);
for i=1:n
    mE(:,i)=vE(m*(i-1)+1:i*m);
end

