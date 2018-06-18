function[out]=rect(x);
%
% rectangle function
%
% evaluates rect(x) 
% note: returns odd number of samples for full width
%
out=abs(x)<=1/2;
out=double(out);
end 