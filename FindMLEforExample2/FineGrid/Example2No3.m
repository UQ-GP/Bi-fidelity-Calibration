clear,clc,format compact 
load Example2WrongMLE.mat PhysData SSE_XTrue XTrue
Case=2;
Dim=3; NL=21^Dim;
Points=(fullfact(21*ones(1,Dim))-1)/20;
AllYh=zeros(NL,1210); AllYl=zeros(NL,1210);

parfor i=1:NL
    AllYh(i,:)=Simulator(Points(i,:),2,Case);
    AllYl(i,:)=Simulator(Points(i,:),1,Case);
end
AllSh=sum((AllYh-PhysData).^2,2); %21^3

AllSh=[AllSh;SSE_XTrue];
Points=[Points;XTrue];

[sortAllSh,sortidx]=sort(AllSh);
save Example2XMLEFineGridOkOptTop20.mat
%%
clear all
load Example2XMLEFineGridOkOptTop20.mat

lb = 0*ones(1,Dim);ub = 1*ones(1,Dim);
options=optimoptions('patternsearch','MaxIterations',10^6,'MeshTolerance',10^-3,'TolFun',10^-6,'MaxFunEvals',10^8);

SSHFun=@(x) sum((Simulator(x,2,Case)-PhysData).^2);
parfor id=1:20
    id
    [XMLETry(id,:),fval(id,:),exitflag(id)] = patternsearch(SSHFun,Points(sortidx(id),:),[],[],[],[],lb,ub,[],options)  ;
end

[SSE_XMLE,minidx]=min(fval);
XMLE=XMLETry(minidx,:);

save Example2XMLEFineGridOkOptTop20.mat
 
SSE_XMLE2=SSE_XMLE;
XMLE2=XMLE;
clearvars -except SSE_XMLE2 XMLE2
save Example2TrueMLE.mat