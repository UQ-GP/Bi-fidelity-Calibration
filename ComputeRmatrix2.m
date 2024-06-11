function R=ComputeRmatrix2(Xs,Thetas,nuggetValue) 
if nargin==2
    nugget=1e-6;
elseif nargin==3
    nugget=nuggetValue;
end
[nX,Dim]=size(Xs);


R=ones(nX,nX);
for id=1:Dim
    x=Xs(:,id);    
    y=Xs(:,id)';
    r=abs(x-y) .*Thetas(id);
    R=R.* [ exp(-r) .* (1+r) ] ;       
end

R=R+ speye(nX)*nugget;

end