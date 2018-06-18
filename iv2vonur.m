function x=iv2vonur(y,N_x,N_y)

YY=reshape(y,N_x,N_y);

x=reshape(YY',N_x*N_y,1);

return