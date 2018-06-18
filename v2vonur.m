function y=v2vonur(x,N_x,N_y)

XX=reshape(x,N_x,N_y);

y=reshape(XX',N_x*N_y,1);

return