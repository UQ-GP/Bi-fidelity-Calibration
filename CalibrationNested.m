function [RecordTable,RecordData,RunTime]=CalibrationNested(DataInput,AccuracyLevel)
% Implements a bi-fidelity calibration method, i.e., the Nested method.
tic
nugget=1e-6;

Dl=DataInput.Dl;
Yl=DataInput.Yl;
Dh=DataInput.Dh;
Yh=DataInput.Yh;
xstar=DataInput.xstar;
w=DataInput.w;
CostRatio=DataInput.CostRatio;
Budget=DataInput.Budget;
Case=DataInput.Case;

[nl,d]=size(Dl);
nh=size(Dh,1);
n=nl+nh;
D=[Dl;Dh];
Level=[ones(nl,1);2*ones(nh,1)];
ZhVec=(sum((Yh-w).^2,2)).^0.5;
ZlVec=(sum((Yl-w).^2,2)).^0.5;

ZVec=[ZlVec;ZhVec];
Budget=Budget-(nl*1+nh*CostRatio);
MinBudget=CostRatio+1;

DlnotDh=Dl(nh+1:nl,:);
YlnotDh=Yl(nh+1:nl,:);
ZlVecnotDh=ZlVec(nh+1:nl,:);
if d==2
    if(AccuracyLevel==1)
    nlevel=2501;
    else
    nlevel=3001;    
    end
elseif d==3
    if(AccuracyLevel==1)
    nlevel=201;
    else
    nlevel=226;
    end
end
GridPoints=(fullfact(nlevel*ones(1,d))-1)/(nlevel-1);

lb=0*ones(1,d); ub=1*ones(1,d);
ZFit=1;

MaxAFVals(n,1)=0;
if(AccuracyLevel==1)
NoS=90;
else
NoS=100;    
end
while (1)
    %Fits the GP model and finds xhat^*_{ML} (estimate of the MLE of the calibration parameter vector).
    if ZFit==1
        [thetal,mul,sigma2l,invRl,invRlRes,condRl,MinM2LogLikelihoodl]=GPFitLF(Dl,ZlVec,nugget,AccuracyLevel);
        [thetah,rho,muh,sigma2h,invRh,invRhRes,condRh,MinM2LogLikelihoodh]=GPFitHF(Dh,ZhVec,ZlVec,nugget,AccuracyLevel);
        ZFit=0;
        sigma2ls(n-1:n,:)=[sigma2l;sigma2l];
        thetals(n-1:n,:)=[thetal;thetal];
        rhos(n-1:n,:)=[rho;rho];
        sigma2hs(n-1:n,:)=[sigma2h;sigma2h];
        thetahs(n-1:n,:)=[thetah;thetah];
        mus(n-1:n,:)=[mul muh; mul muh];
        
        MinM2LogLikelihoodhs(n-1:n,:)=MinM2LogLikelihoodh;
        MinM2LogLikelihoodls(n-1:n,:)=MinM2LogLikelihoodl;
        MinM2LogLikelihoods(n-1:n,:)=MinM2LogLikelihoodh+MinM2LogLikelihoodl;

    else
        [invRl,invRlRes,condRl]=M2LogLikelihoodLFSame(Dl,ZlVec,thetal,mul,nugget);

        [invRh,invRhRes,condRh]=M2LogLikelihoodHFSame(Dh,ZhVec,ZlVec,thetah,rho,muh,nugget);
        
    end

    [ZhhatxhatstarML,minidx]=min(ZhVec);
    xhatstarML=Dh(minidx,:);
    
    xhatstarMLs(n-1:n,:)=[xhatstarML;xhatstarML];
    Zhminhats(n-1:n,:)=[ZhhatxhatstarML;ZhhatxhatstarML];
    
    %Evaluates the true HF SSE at xhat^*_{ML}, i.e., S_h(xhat^*_{ML}).
    yhxhatstarML=Yh(minidx,:);
    ShxhatstarML=sum((yhxhatstarML-w).^2);
    yhxhatstarMLs(n-1:n,:)=[yhxhatstarML;yhxhatstarML];
    ShxhatstarMLs(n-1:n,:)=[ShxhatstarML;ShxhatstarML];
    condRlRhs(n-1:n,:)=[condRl,condRh;condRl,condRh];

    disp(['Based on ' num2str(n) ' design points, the estimate xhat^*_{ML} of S_h''s minimizer and S_h(xhat^*_{ML}) are [' num2str(xhatstarMLs(n,:),' %.3f') '] and ' num2str(ShxhatstarMLs(n,:)) ' respectively.'])
    
    if Budget<MinBudget
        break
    end
 
    MinusEIVals=CompMinusEI(GridPoints,Dh,Dl,thetal,mul,sigma2l,invRl,invRlRes,thetah,rho,muh,sigma2h,invRh,invRhRes,ZhhatxhatstarML,nugget);
    %Adds the follow-up design points by maximizing the EI AF for the HF response.
    [~,SortedMinusEIValsIndices]=sort(MinusEIVals);
    options=optimoptions('patternsearch','Display','off');
    
    MinusEIFun=@(x) CompMinusEI(x,Dh,Dl,thetal,mul,sigma2l,invRl,invRlRes,thetah,rho,muh,sigma2h,invRh,invRhRes,ZhhatxhatstarML,nugget);
    XBestTry=zeros(NoS,d);
    fBestTry=zeros(NoS,1); SortedGridPoints=GridPoints(SortedMinusEIValsIndices,:);     
    parfor id=1:NoS
        [XBestTry(id,:),fBestTry(id,1)]=patternsearch(MinusEIFun,SortedGridPoints(id,:),[],[],[],[],lb,ub,[],options);
    end
    [fBest,minidx]=min(fBestTry);
    NextPoint=XBestTry(minidx,:);
    MaxAFVal=-fBest;
    MaxAFVals([n+1 n+2],:)=MaxAFVal;
    disp(['Remaining budget=' num2str(Budget) '. Next, the ' num2str(n+1) '-th and ' num2str(n+2) '-th runs will be made at the point [' num2str(NextPoint,' %1.3f') '].'])
    
    Dh(nh+1,:)=NextPoint;
    Yh(nh+1,:)=Simulator(NextPoint,2,Case);
    ZhVec(nh+1,:)=sum((Yh(nh+1,:)-w).^2)^0.5;
    
    Dl=[Dh;DlnotDh];
    Yl_new=Simulator(NextPoint,1,Case);
    Yl=[Yl(1:nh,:);Yl_new;YlnotDh];
    Zl_new=sum((Yl_new-w).^2)^0.5;
    ZlVec=[ZlVec(1:nh,:);Zl_new;ZlVecnotDh];
    
    D=[D;NextPoint;NextPoint];
    Level=[Level;2;1];
    ZVec=[ZVec;ZhVec(nh+1,:);Zl_new];
    
    Budget=Budget-CostRatio-1;
    nh=nh+1;
    nl=nl+1;
    n=n+2;
end
sigma2ls(n-1:n,:)=[sigma2l;sigma2l];
thetals(n-1:n,:)=[thetal;thetal];
rhos(n-1:n,:)=[rho;rho];
sigma2hs(n-1:n,:)=[sigma2h;sigma2h];
thetahs(n-1:n,:)=[thetah;thetah];
mus(n-1:n,:)=[mul muh; mul muh];
MinM2LogLikelihoodhs(n-1:n,:)=MinM2LogLikelihoodh;
MinM2LogLikelihoodls(n-1:n,:)=MinM2LogLikelihoodl;
MinM2LogLikelihoods(n-1:n,:)=MinM2LogLikelihoodh+MinM2LogLikelihoodl;

%Stores design points, the simulator outputs at design points, and other information.
RecordData.Dl=Dl;
RecordData.Yl=Yl;
RecordData.ZlVec=ZlVec;
RecordData.Dh=Dh;
RecordData.Yh=Yh;
RecordData.ZhVec=ZhVec;
RecordData.xstar=xstar;
RecordData.w=w;
RecordData.CostRatio=CostRatio;
RecordData.Budget=Budget;
RecordData.yhxhatstarMLs=yhxhatstarMLs;

%Stores GP emulator parameter estimates, values of xhat^*_{ML} (rows of xhatstarMLs), values of S_h(xhat^*_{ML}) (elements of ShxhatstarMLs), and other important information with a table.
RecordTable=table(D,Level,ZVec,MaxAFVals,xhatstarMLs,Zhminhats,ShxhatstarMLs,sigma2ls,thetals,rhos,sigma2hs,thetahs,mus,condRlRhs,MinM2LogLikelihoodhs,MinM2LogLikelihoodls,MinM2LogLikelihoods);
RunTime=toc;
end
%%
function [thetal,mul,sigma2l,invRl,invRlRes,condRl,MinM2LogLikelihoodl]=GPFitLF(Dl,ZlVec,nugget,AccuracyLevel)

d=size(Dl,2);
lb=(0.25)*ones(1,d);
ub=(15)*ones(1,d);
Ranges=ub-lb;
npar=numel(lb);
M2LogLlFun=@(Parl) M2LogLikelihoodLF(Dl,ZlVec,Parl,nugget);
if(AccuracyLevel==1)
HNoS=90;
else
HNoS=100;
end

Sobolset=sobolset(npar,'Skip',1e3,'Leap',1e2);
if(AccuracyLevel==1)
StandardPoints=net(Sobolset,7000*npar);
else
StandardPoints=net(Sobolset,8000*npar);  
end
CandidatePoints=lb+Ranges.*StandardPoints;

NoCandidatePoints=size(CandidatePoints,1);
M2LogLlVals=zeros(NoCandidatePoints,1);
parfor id=1:NoCandidatePoints
    M2LogLlVals(id,1)=M2LogLlFun(CandidatePoints(id,:));
end
[~,SortedM2LogLlValsIndices]=sort(M2LogLlVals);
options=optimoptions('patternsearch','Display','off');

Selected_StandardPoints=StandardPoints(SortedM2LogLlValsIndices(1:HNoS),:);
Remaining_StandardPoints=StandardPoints(SortedM2LogLlValsIndices((HNoS+1):end),:);
Count=0;
for kd=1:size(Remaining_StandardPoints,1)
    if min(pdist2(Remaining_StandardPoints(kd,:),Selected_StandardPoints),[],2)>sqrt(npar*0.1^2)
        Selected_StandardPoints=[Selected_StandardPoints; Remaining_StandardPoints(kd,:)];
        Count=Count+1;
        if Count==HNoS
            break
        end
    end
end

StartingPoints=lb+Ranges.*Selected_StandardPoints; 

NoStartingPoints=size(StartingPoints,1);
fBestTry=zeros(NoStartingPoints,1);
XBestTry=zeros(NoStartingPoints,npar);
parfor id=1:NoStartingPoints
    [XBestTry(id,:),fBestTry(id,1)]=patternsearch(M2LogLlFun,StartingPoints(id,:),[],[],[],[],lb,ub,[],options); 
end
[~,minidx]=min(fBestTry);
OptParl=XBestTry(minidx,:);
[MinM2LogLikelihoodl,mul,sigma2l,thetal,invRl,invRlRes,condRl]=M2LogLlFun(OptParl);

end
%% 
function [M2LogLikelihoodlVal,mul,sigma2l,thetal,invRl,invRlRes,condRl]=M2LogLikelihoodLF(Dl,ZlVec,thetal,nugget)
nl=size(ZlVec,1);
FlT=ones(1,nl);
Rl=ComputeRmatrix2(Dl,thetal,nugget);
[invRl,logdetRl,condRl]=invandlogdet(Rl);
FlTinvRl=FlT*invRl;

mul=(FlTinvRl*ZlVec)/sum(invRl,'all');
Res=ZlVec-mul;
invRlRes=invRl*Res;
sigma2l=Res'*invRlRes/nl;
M2LogLikelihoodlVal=nl*log(sigma2l)+logdetRl;
if ~isfinite(mul) || sigma2l<=0 || ~isfinite(sigma2l)
    M2LogLikelihoodlVal=Inf;mul=[];sigma2l=[];thetal=[];invRl=[];invRlRes=[];condRl=[];
    return
end
end
%%
function [invRl,invRlRes,condRl]=M2LogLikelihoodLFSame(Dl,ZlVec,thetal,mul,nugget)
Rl=ComputeRmatrix2(Dl,thetal,nugget);
[invRl,~,condRl]=invandlogdet(Rl);
Res=ZlVec-mul;
invRlRes=invRl*Res;

if any(~isfinite(invRlRes),'all') 
    invRl=[];invRlRes=[];condRl=[];
    return
end

end
%%
function [thetah,rho,muh,sigma2h,invRh,invRhRes,condRh,MinM2LogLikelihoodh]=GPFitHF(Dh,ZhVec,ZlVec,nugget,AccuracyLevel)

d=size(Dh,2);
lb=(0.25)*ones(1,d);
ub=(15)*ones(1,d);
Ranges=ub-lb;
npar=numel(lb);
if(AccuracyLevel==1)
    HNoS=90;
else
    HNoS=100;
end

M2LogLhFun=@(Parh) M2LogLikelihoodHF(Dh,ZhVec,ZlVec,Parh,nugget);

Sobolset=sobolset(npar,'Skip',1e3,'Leap',1e2);
if(AccuracyLevel==1)
StandardPoints=net(Sobolset,7000*npar);
else
StandardPoints=net(Sobolset,8000*npar);   
end
CandidatePoints=lb+Ranges.*StandardPoints;

NoCandidatePoints=size(CandidatePoints,1);
M2LogLhVals=zeros(NoCandidatePoints,1);
parfor id=1:NoCandidatePoints
    M2LogLhVals(id,1)=M2LogLhFun(CandidatePoints(id,:));
end
[~,SortedM2LogLhValsIndices]=sort(M2LogLhVals);
options=optimoptions('patternsearch','Display','off');

Selected_StandardPoints=StandardPoints(SortedM2LogLhValsIndices(1:HNoS),:);
Remaining_StandardPoints=StandardPoints(SortedM2LogLhValsIndices((HNoS+1):end),:);
Count=0;
for kd=1:size(Remaining_StandardPoints,1)
    if min(pdist2(Remaining_StandardPoints(kd,:),Selected_StandardPoints),[],2)>sqrt(npar*0.1^2)
        Selected_StandardPoints=[Selected_StandardPoints; Remaining_StandardPoints(kd,:)];
        Count=Count+1;
        if Count==HNoS
            break
        end
    end
end

StartingPoints=lb+Ranges.*Selected_StandardPoints; 

NoStartingPoints=size(StartingPoints,1);
fBestTry=zeros(NoStartingPoints,1);
XBestTry=zeros(NoStartingPoints,npar);
parfor id=1:NoStartingPoints
    [XBestTry(id,:),fBestTry(id,1)]=patternsearch(M2LogLhFun,StartingPoints(id,:),[],[],[],[],lb,ub,[],options);
end
[~,minidx]=min(fBestTry);
OptParh=XBestTry(minidx,:);
[MinM2LogLikelihoodh,rho,muh,sigma2h,thetah,invRh,invRhRes,condRh]=M2LogLhFun(OptParh);

end
%%
function [M2LogLikelihoodhVal,rho,muh,sigma2h,thetah,invRh,invRhRes,condRh]=M2LogLikelihoodHF(Dh,ZhVec,ZlVec,thetah,nugget)
nh=size(ZhVec,1);
Fh=[ZlVec(1:nh,:) ones(nh,1)];
Rh=ComputeRmatrix2(Dh,thetah,nugget);
[invRh,logdetRh,condRh]=invandlogdet(Rh);
FhTinvRh=Fh'*invRh;

betah=(FhTinvRh*Fh)\(FhTinvRh*ZhVec);
rho=betah(1);
muh=betah(2);
Res=ZhVec-Fh*betah;
invRhRes=invRh*Res;
sigma2h=Res'*invRhRes/nh;
M2LogLikelihoodhVal=nh*log(sigma2h)+logdetRh;
if any(~isfinite(betah),'all') || sigma2h<=0 || ~isfinite(sigma2h)
    M2LogLikelihoodhVal=Inf;rho=[];muh=[];sigma2h=[];thetah=[];invRh=[];invRhRes=[];condRh=[];
    return
end
end
%%
function [invRh,invRhRes,condRh]=M2LogLikelihoodHFSame(Dh,ZhVec,ZlVec,thetah,rho,muh,nugget)
nh=size(ZhVec,1);
Fh=[ZlVec(1:nh,:) ones(nh,1)];
Rh=ComputeRmatrix2(Dh,thetah,nugget);
[invRh,~,condRh]=invandlogdet(Rh);
betah=[rho,muh]';
Res=ZhVec-Fh*betah;
invRhRes=invRh*Res;

if any(~isfinite(invRhRes),'all') 
    invRh=[];invRhRes=[];condRh=[];
    return
end

end
%%
function MinusEIVals=CompMinusEI(xs,Dh,Dl,thetal,mul,sigma2l,invRl,invRlRes,thetah,rho,muh,sigma2h,invRh,invRhRes,ZhhatxhatstarML,nugget)

Noxs=size(xs,1);
rlT=ComputeRmatrix(xs,Dl,thetal);
Zlhats=mul+rlT*invRlRes;

rhT=ComputeRmatrix(xs,Dh,thetah);
fh=[Zlhats ones(Noxs,1)];
betah=[rho;muh];
Zhhats=fh*betah+rhT*invRhRes;

ZlVars=sigma2l*(1+nugget-sum((rlT*invRl).*rlT,2));
ZlVars=max(ZlVars,0);
ZhVars=rho^2*ZlVars+sigma2h*(1+nugget-sum((rhT*invRh).*rhT,2));
ZhVars=max(ZhVars,0);

sVals=ZhVars.^0.5;
DeltaVals=ZhhatxhatstarML-Zhhats;
DeltaValsoversVals=DeltaVals./sVals;
EIVals=DeltaVals.*normcdf(DeltaValsoversVals)+sVals.*normpdf(DeltaValsoversVals);
EIVals( (ZhVars==0) | (min(pdist2(xs,Dl),[],2) < 10^(-3)) )=0;
MinusEIVals=-EIVals;
end