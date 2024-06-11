function SS=TransformData_inv(Z,phi,ZNBC) 

if ZNBC==0      %Indentity transformation
    SS=Z;
%     SS(Z<=0)=0;
    
elseif ZNBC==1  %Box-Cox transformation
    if phi==0
        SS=exp(Z);
    elseif phi>0
        SS=(Z*phi+1).^(1/phi);
        SS(Z<-1/phi)=0;
    elseif phi<0
        SS=(Z*phi+1).^(1/phi);
        SS(Z>=-1/phi)=Inf;
    end
    
elseif ZNBC==2 %Squared root transformation
    SS=Z.^2;
%     SS(Z<=0)=0;
end

end