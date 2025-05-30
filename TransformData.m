function [Z,LogAbsJacobian]=TransformData(SS,phi,ID_BC_or_SR) 

if any(SS<=0,'all')
    return
end

if ID_BC_or_SR==0 %Identity transformation
    Z=SS;   
    dZ=1;    
elseif ID_BC_or_SR==1 %Box-Cox transformation
    if phi~=0
        Z=(SS.^phi-1)/phi;
        dZ=SS.^(phi-1);
    elseif phi==0
        Z=log(SS);    
        dZ=1./SS;
    end
elseif ID_BC_or_SR==2 %Square root transformation
    Z=SS.^0.5;    
    dZ=1;
end
LogAbsJacobian=sum(log(dZ)); %log of absolute Jacobian = log of Jacobian

end