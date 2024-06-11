function R=ComputeRmatrix(Xs,Ys,Thetas) 
[nX,Dim]=size(Xs);
[nY,~]=size(Ys);

R=ones(nX,nY);
for id=1:Dim
    x=Xs(:,id);    
    y=Ys(:,id)';
    r=abs(x-y) .*Thetas(id);
    R=R.* [ exp(-r) .* (1+r) ];      
end

end