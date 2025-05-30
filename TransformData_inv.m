function SS=TransformData_inv(Z,phi,ID_BC_or_SR) 

if ID_BC_or_SR==0 %Identity transformation
    SS=Z;
    
elseif ID_BC_or_SR==1 %Box-Cox transformation
    if phi==0
        SS=exp(Z);
    elseif phi>0
        SS=(Z*phi+1).^(1/phi);
        SS(Z<-1/phi)=0;
    elseif phi<0
        SS=(Z*phi+1).^(1/phi);
        SS(Z>=-1/phi)=Inf;
    end
    
elseif ID_BC_or_SR==2 %Square root transformation
    SS=Z.^2;

end

end