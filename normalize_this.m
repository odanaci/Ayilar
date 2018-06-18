function out=normalize_this(F)

out=(F-min(F(:)))/(max(F(:))-min(F(:)));

return