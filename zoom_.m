function Y=zoom_(X,zoom1,zoom2)

[M,N]=size(X);

if M==1|N == 1
    
    if M>1
    Y=X(round(zoom1*M):round(zoom2*M),1);
    else
     Y=X(1,round(zoom1*N):round(zoom2*N));
    end

else

Y=X(round(zoom1*M):round(zoom2*M),round(zoom1*M):round(zoom2*M));
        
    end


end