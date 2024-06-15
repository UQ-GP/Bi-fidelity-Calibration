function [Dl,Dh]=GenerateNestedLHD(nl,nh,Dim,nsim) 

t=nl/nh;
if(round(t)~=t)
    disp('error')
    return
end
BestDist=0; BestDesign=[ ];

for idx=1:nsim
    D=zeros(nl,Dim);

    for idxDim=1:Dim
        V=randperm(nh)';
        Vl=(V-1)*t;
        Select=unifrnd(0,t,nh,1);
        Select=ceil(Select);
        tau=Vl+Select;
        remain=setdiff(1:nl,tau);
        rho=remain(randperm(nl-nh));
        D(:,idxDim)=[tau(:);rho(:)];
    end
    Design=(D-unifrnd(0,1,nl,Dim))/nl;
    Dist=min(pdist(Design));
    if(Dist>BestDist)
        BestDist=Dist;
        BestDesign=Design;
    end
end
Dl=BestDesign;
Dh=BestDesign(1:nh,:);
end%GenerateNestedLHD

