function out=fourt2(A)

out=ifftshift(fft2(fftshift(A)));

return