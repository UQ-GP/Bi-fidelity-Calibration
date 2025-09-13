clear all
tic 
 Case=2;
 GridPoints=(fullfact([41 41 41])-1)/40;
 parfor jd=1:size(GridPoints,1)
     jd
     AllYh(jd,:)=Simulator(GridPoints(jd,:),2,Case);
     AllYl(jd,:)=Simulator(GridPoints(jd,:),1,Case);
 end
 RunTime=toc; 
 save Example2GridData.mat AllYh AllYl GridPoints RunTime
%% 
clear all
tic
load Example2.mat
load Example2GridData.mat AllYh GridPoints RunTime
format long g
t=0.99; d=3;
NGP=size(GridPoints,1);
lb=0*ones(1,d); ub=1*ones(1,d);
ApproxShsGrid=zeros(NGP,100);
AllSh=sum((AllYh-SingleDataInput(1).w).^2,2); 
for id=1:1:100
    id
    DataInput=SingleDataInput(id);
    Dh=DataInput.Dh;
    Yh=(DataInput.Yh)';
    
    xstar=DataInput.xstar;
    w=(DataInput.w)'; 
    Case=DataInput.Case;
    ybarh=mean(Yh,2);
    Yhc=Yh-ybarh;
    [U0,S0,~]=svd(Yhc);
    S0=diag(S0).^2;
    PossibleFractions=cumsum(S0)/sum(S0);
    p=find(PossibleFractions>=t,1);
    U=U0(:,1:p);
    
    ApproxShFun=@(x) sum((U*U'*(Simulator(x,2,Case)'-ybarh)+ybarh-w).^2);
    TrueSh=@(x) sum((Simulator(x,2,Case)'-w).^2);
    
    xhatstarMLEnd=T_SVD{id}.xhatstarMLs(end,:);
    xhatstarMLsEnd(id,:)=xhatstarMLEnd;
    TrueShxhatstarMLsEnd(id,1)=T_SVD{id}.ShxhatstarMLs(end,:);
    ApproxShxhatstarMLsEnd(id,1)=ApproxShFun(xhatstarMLEnd);
        
    ApproxShFun_y=@(y) sum((U*U'*(y'-ybarh)+ybarh-w).^2);
    parfor jd=1:NGP
        ApproxShsGrid(jd,id)=ApproxShFun_y(AllYh(jd,:));
    end
    [SortedApproxShs0,SortedApproxShsIndices0]=sort(ApproxShsGrid(:,id));
    PrctileApproxShs(id,:)=[SortedApproxShs0(1:100)' SortedApproxShs0(ceil(41^3/2)) SortedApproxShs0(end)];
    SmallestApproxShsPoints(:,:,id)=GridPoints(SortedApproxShsIndices0(1:100),:);
    
    CandidatePoints=[GridPoints;xstar;xstarML;xhatstarMLsEnd(id,:)];
    ApproxShxstar=ApproxShFun(xstar);
    ApproxShxstarML=ApproxShFun(xstarML);
    ApproxShs=[ApproxShsGrid(:,id);ApproxShxstar;ApproxShxstarML;ApproxShxhatstarMLsEnd(id,1)];
    [SortedApproxShs,SortedApproxShsIndices]=sort(ApproxShs); 
    options=optimoptions('fmincon','Display','off'); TiesInd(id)=SortedApproxShs(21)-SortedApproxShs(20);
    SortedCandidatePoints=CandidatePoints(SortedApproxShsIndices,:); StartPoints(:,:,id)=SortedCandidatePoints(1:20,:);
    parfor jd=1:20
        jd
        StartPoint=StartPoints(jd,:,id);
        [ApproxxstarMLTry(jd,:,id),ApproxShBestValTry(jd,id),exitflag(jd,id),~,~,~,hessian(:,:,jd,id)]=fmincon(ApproxShFun,StartPoint,[],[],[],[],lb,ub,[],options);
    end
    
    [ApproxShApproxxstarML,minidx]=min(ApproxShBestValTry(:,id)); 
    ApproxxstarML=ApproxxstarMLTry(minidx,:,id);
    store1(id,:)=[xstar Shxstar ApproxShxstar]
    store2(id,:)=[xstarML ShxstarML ApproxShxstarML]
    TrueShApproxxstarML=TrueSh(ApproxxstarML);
    store3(id,:)=[ApproxxstarML TrueShApproxxstarML ApproxShApproxxstarML]
    store4(id,:)=[xhatstarMLsEnd(id,:) TrueShxhatstarMLsEnd(id,:) ApproxShxhatstarMLsEnd(id,:)] 
    
    xhatstarMLs=T_SVD{id}.xhatstarMLs(nh0:end,:);
    TrueShxhatstarMLs=T_SVD{id}.ShxhatstarMLs(nh0:end,:)';
    NoxhatstarMLs=size(xhatstarMLs,1);
    for jd=1:NoxhatstarMLs
        ApproxShxhatstarMLs(jd)=ApproxShFun(xhatstarMLs(jd,:)); 
        DisttoApproxxstarML(jd)=norm(xhatstarMLs(jd,:)-store3(id,1:d));
        DisttoxstarML(jd)=norm(xhatstarMLs(jd,:)-xstarML);
        NumberGridPointsBetterApproxSh(jd)=sum(SortedApproxShs0<ApproxShxhatstarMLs(jd));     
        NumberGridPointsBetterSh(jd)=sum(AllSh<TrueShxhatstarMLs(jd));         
    end
    store5(id,:)=ApproxShxhatstarMLs;     
    store6(id,:)=TrueShxhatstarMLs;
    store7(id,:)=DisttoApproxxstarML;
    store8(id,:)=DisttoxstarML; 
    store9(id,:)=NumberGridPointsBetterApproxSh;  
    store10(id,:)=NumberGridPointsBetterSh;    
end
RunTime=RunTime+toc;
clearvars -except RunTime store1 store2 store3 store4 store5 store6 store7 store8 store9 store10 PrctileApproxShs SmallestApproxShsPoints AllSh TiesInd StartPoints ApproxxstarMLTry ApproxShBestValTry exitflag hessian
save DataforAppendixJ_SVD.mat