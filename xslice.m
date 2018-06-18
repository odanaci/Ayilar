function out=xslice(A)

[M,N]=size(A);
out=A(:,M/2+1);

return