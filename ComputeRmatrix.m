function R=ComputeRmatrix(Xs,Ys,theta) 
[nXs,d]=size(Xs);
[nYs,~]=size(Ys);

R=ones(nXs,nYs);
for i=1:d
    x=Xs(:,i);    
    y=Ys(:,i)';
    r=abs(x-y).*theta(i);
    R=R.*(exp(-r).*(1+r));      
end

end