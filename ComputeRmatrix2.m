function R=ComputeRmatrix2(Xs,theta,nuggetValue) 
if nargin==2
    nugget=1e-6;
elseif nargin==3
    nugget=nuggetValue;
end
[nXs,d]=size(Xs);

R=ones(nXs,nXs);
for i=1:d
    x=Xs(:,i);    
    y=Xs(:,i)';
    r=abs(x-y).*theta(i);
    R=R.*(exp(-r).*(1+r));       
end

R=R+speye(nXs)*nugget;

end