function [xstarML,ShxstarML,exitflag]=Example2FindTrueMLE(w,Shxstar,xstar)
Case=2;
d=3; nlevels=21; 
GridPoints=(fullfact(nlevels*ones(1,d))-1)/(nlevels-1);
NGP=nlevels^d;
AllYh=zeros(NGP,1210);
parfor i=1:NGP
    AllYh(i,:)=Simulator(GridPoints(i,:),2,Case);
end
AllSh=sum((AllYh-w).^2,2);  

AllSh=[AllSh;Shxstar];
GridPoints=[GridPoints;xstar];

[~,sortidx]=sort(AllSh);

lb=0*ones(1,d); ub=1*ones(1,d);
options=optimoptions('fmincon','Display','off');

ShFun=@(x) sum((Simulator(x,2,Case)-w).^2);
xstarMLTry=zeros(20,d); ShFunBestValTry=zeros(20,1); exitflagTry=zeros(20,1); SortedGridPoints=GridPoints(sortidx,:);
parfor id=1:20
    id
    [xstarMLTry(id,:),ShFunBestValTry(id,:),exitflagTry(id)]=fmincon(ShFun,SortedGridPoints(id,:),[],[],[],[],lb,ub,[],options);
end

[~,minidx]=min(ShFunBestValTry);
xstarML0=xstarMLTry(minidx,:);

[xstarML,ShxstarML,exitflag]=fmincon(ShFun,xstarML0,[],[],[],[],lb,ub,[],options);

save Example2FindTrueMLEData.mat
save('Example2TrueMLE.mat','xstarML','ShxstarML','exitflag')