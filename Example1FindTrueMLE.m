clear all
load Example1GridData26.mat
AllSh=sum((AllYh-w).^2,2);
[~,sortidx]=sort(AllSh);
%%
lb=0*ones(1,d); ub=1*ones(1,d);
options=optimoptions('patternsearch','MaxIterations',500,'Display','iter');

ShFun=@(x) sum((Simulator(x,2,Case)-w).^2);
xstarMLTry=zeros(20,d); ShFunBestValTry=zeros(20,1); exitflagTry=zeros(20,1); 
for id=1:20
    id
    [xstarMLTry(id,:),ShFunBestValTry(id,:),exitflagTry(id)]=patternsearch(ShFun,GridPoints(sortidx(id),:),[],[],[],[],lb,ub,[],options);
    save Example1TrueMLE.mat
end

[ShxstarML,minidx]=min(ShFunBestValTry);
xstarML=xstarMLTry(minidx,:);

save Example1TrueMLE.mat