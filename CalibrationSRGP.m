function [RecordTable,RecordData,RunTime]=CalibrationSRGP(DataInput,AccuracyLevel)
% Implements a single-fidelity calibration method, i.e., the SR-GP method.
tic
nugget=1e-6;

Dh=DataInput.Dh;
Yh=DataInput.Yh;
xstar=DataInput.xstar;
w=DataInput.w;
CostRatio=DataInput.CostRatio;
Budget=DataInput.Budget;
Case=DataInput.Case;
[n,d]=size(Dh);
Level=2*ones(n,1);

ZhVec=(mean((Yh-w).^2,2)).^0.5;
Budget=Budget-(n*CostRatio);

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

HistoryxhatstarMLs=[];
MaxAFVals(n,1)=0;
if(AccuracyLevel==1)
NoS=90;
else
NoS=100;
end
while (1)
    %Fits the GP model and finds xhat^*_{ML} (estimate of the MLE of the calibration parameter vector).
    if ZFit==1
        [sigma2,theta,mu,invR,invRRes,condR,MinM2LogLikelihood]=GPFit(Dh,ZhVec,nugget,AccuracyLevel);
        ZFit=0;
        sigma2s(n,:)=sigma2;                         thetas(n,:)=theta; mus(n,:)=mu;
        MinM2LogLikelihoods(n,:)=MinM2LogLikelihood;

        [Zhhats1,ZhVars1,rT]=CompPosMeanVar(GridPoints,Dh,theta,mu,sigma2,invR,invRRes,nugget);
        
    else
        [invR,invRRes,condR]=M2LogLikelihoodSame(Dh,ZhVec,theta,mu,nugget);

        rT=[rT,ComputeRmatrix(GridPoints,Dh(n,:),theta)];
        [Zhhats1,ZhVars1]=CompPosMeanVarGrid(mu,sigma2,invR,invRRes,nugget,rT);
    end
    
    xhatstarMLCandidatePoints=[GridPoints;Dh;HistoryxhatstarMLs];
    
    [Zhhats2,ZhVars2]=CompPosMeanVar(Dh,Dh,theta,mu,sigma2,invR,invRRes,nugget);
    [Zhhats3,ZhVars3]=CompPosMeanVar(HistoryxhatstarMLs,Dh,theta,mu,sigma2,invR,invRRes,nugget);
    Zhhats=[Zhhats1;Zhhats2;Zhhats3]; ZhVars=[ZhVars1;ZhVars2;ZhVars3];
    ZhQuantileVals=Zhhats+norminv(0.9)*ZhVars.^0.5;
    [~,SortedZhQuantileValsIndices]=sort(ZhQuantileVals);
    options=optimoptions('patternsearch','Display','off');  
    
    ZhQuantileFun=@(x) CompZhQuantile(x,Dh,theta,mu,sigma2,invR,invRRes,nugget);
    XBestTry=zeros(NoS,d);
    fBestTry=zeros(NoS,1); SortedxhatstarMLCandidatePoints=xhatstarMLCandidatePoints(SortedZhQuantileValsIndices,:);   
    parfor id=1:NoS
        [XBestTry(id,:),fBestTry(id,1)]=patternsearch(ZhQuantileFun,SortedxhatstarMLCandidatePoints(id,:),[],[],[],[],lb,ub,[],options);
    end
    [Zhminhats(n,:),minidx]=min(fBestTry);
    xhatstarML=XBestTry(minidx,:);
    xhatstarMLs(n,:)=xhatstarML;    
    condRs(n,:)=condR;
    %Evaluates the true HF SSE at xhat^*_{ML}, i.e., S_h(xhat^*_{ML}).
    yhxhatstarMLs(n,:)=Simulator(xhatstarMLs(n,:),2,Case);
    ShxhatstarMLs(n,:)=sum((yhxhatstarMLs(n,:)-w).^2);
    
    disp(['Based on ' num2str(n) ' design points, the estimate xhat^*_{ML} of S_h''s minimizer and S_h(xhat^*_{ML}) are [' num2str(xhatstarMLs(n,:),' %.3f') '] and ' num2str(ShxhatstarMLs(n,:)) ' respectively.'])
    
    if Budget<CostRatio
        break
    end
    HistoryxhatstarMLs=[HistoryxhatstarMLs;xhatstarML];
    MaxAFCandidatePoints=[GridPoints;HistoryxhatstarMLs];

    [ZhhatxhatstarML,ZhVars4]=CompPosMeanVar(xhatstarML,Dh,theta,mu,sigma2,invR,invRRes,nugget);
    ZhhatsAF=[Zhhats1;Zhhats3;ZhhatxhatstarML]; ZhVarsAF=[ZhVars1;ZhVars3;ZhVars4]; 

    %Adds a follow-up design point by maximizing the EI AF.
    MinusEIVals=CompMinusEIGrid(MaxAFCandidatePoints,Dh,ZhhatxhatstarML,ZhhatsAF,ZhVarsAF);
    [~,SortedMinusEIValsIndices]=sort(MinusEIVals);
    options=optimoptions('patternsearch','Display','off');

    MinusEIFun=@(xs) CompMinusEI(xs,Dh,ZhhatxhatstarML,theta,mu,sigma2,invR,invRRes,nugget);    
    XBestTry=zeros(NoS,d);
    fBestTry=zeros(NoS,1); SortedMaxAFCandidatePoints=MaxAFCandidatePoints(SortedMinusEIValsIndices,:);    
    parfor id=1:NoS
        [XBestTry(id,:),fBestTry(id,1)]=patternsearch(MinusEIFun,SortedMaxAFCandidatePoints(id,:),[],[],[],[],lb,ub,[],options);
    end
    [fBest,minidx]=min(fBestTry);
    NextPoint=XBestTry(minidx,:);
    MaxAFVal=-fBest; MaxAFVals(n+1,:)=MaxAFVal;
    disp(['Remaining budget=' num2str(Budget) '. Next, the ' num2str(n+1) '-th run will be made at the point [' num2str(NextPoint,' %1.3f') '].'])
    n=n+1;
    Dh(n,:)=NextPoint;
    Yh(n,:)=Simulator(NextPoint,2,Case);
    Level(n,:)=2;
    
    ZhVec(n,:)=(mean((Yh(n,:)-w).^2))^0.5;
    Budget=Budget-CostRatio;
end

sigma2s(n,:)=sigma2;                         thetas(n,:)=theta; mus(n,:)=mu;
MinM2LogLikelihoods(n,:)=MinM2LogLikelihood;

%Stores design points, the simulator outputs at design points, and other information.
D=Dh;
RecordData.Dl=[];
RecordData.Yl=[];
RecordData.Dh=Dh;
RecordData.Yh=Yh;
RecordData.ZhVec=ZhVec;
RecordData.xstar=xstar;
RecordData.w=w;
RecordData.CostRatio=CostRatio;
RecordData.Budget=Budget;
RecordData.yhxhatstarMLs=yhxhatstarMLs;

%Stores GP emulator parameter estimates, values of xhat^*_{ML} (rows of xhatstarMLs), values of S_h(xhat^*_{ML}) (elements of ShxhatstarMLs), and other important information with a table.
RecordTable=table(D,Level,ZhVec,MaxAFVals,xhatstarMLs,Zhminhats,ShxhatstarMLs,sigma2s,thetas,mus,condRs,MinM2LogLikelihoods);
RunTime=toc;
end
%%
function [sigma2,theta,mu,invR,invRRes,condR,MinM2LogLikelihood]=GPFit(Dh,ZhVec,nugget,AccuracyLevel)
[n,d]=size(Dh);

lb=(0.25)*ones(1,d);
ub=(15)*ones(1,d);
Ranges=ub-lb;
M2LogLFun=@(Par) M2LogLikelihood(Dh,ZhVec,Par,n,nugget);
npar=numel(ub);
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
M2LogLVals=zeros(NoCandidatePoints,1);
parfor id=1:NoCandidatePoints
    M2LogLVals(id,1)=M2LogLFun(CandidatePoints(id,:));
end
[~,SortedM2LogLValsIndices]=sort(M2LogLVals);
options=optimoptions('patternsearch','Display','off');

Selected_StandardPoints=StandardPoints(SortedM2LogLValsIndices(1:HNoS),:);
Remaining_StandardPoints=StandardPoints(SortedM2LogLValsIndices((HNoS+1):end),:);
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
    [XBestTry(id,:),fBestTry(id,1)]=patternsearch(M2LogLFun,StartingPoints(id,:),[],[],[],[],lb,ub,[],options);
end
[~,minidx]=min(fBestTry);
OptPar=XBestTry(minidx,:);
[MinM2LogLikelihood,mu,sigma2,theta,invR,invRRes,condR]=M2LogLFun(OptPar);

end
%%
function [M2LogLikelihoodVal,mu,sigma2,theta,invR,invRRes,condR]=M2LogLikelihood(Dh,ZhVec,theta,n,nugget)

FT=ones(1,n);
R=ComputeRmatrix2(Dh,theta,nugget);
[invR,logdetR,condR]=invandlogdet(R);
FTinvR=FT*invR;
mu=(FTinvR*ZhVec)/sum(invR,'all');

Res=ZhVec-mu;
invRRes=invR*Res;
sigma2=Res'*invRRes/n;
M2LogLikelihoodVal=n*log(sigma2)+logdetR;

if ~isfinite(mu) || sigma2<=0 || ~isfinite(sigma2)
    M2LogLikelihoodVal=Inf;mu=[];sigma2=[];theta=[];invR=[];invRRes=[];condR=[];
    return
end
end
%%
function [invR,invRRes,condR]=M2LogLikelihoodSame(Dh,ZhVec,theta,mu,nugget)

R=ComputeRmatrix2(Dh,theta,nugget);
[invR,~,condR]=invandlogdet(R);
Res=ZhVec-mu;
invRRes=invR*Res;

if any(~isfinite(invRRes),'all') 
    invR=[];invRRes=[];condR=[];
    return
end

end
%%
function ZhQuantiles=CompZhQuantile(xs,Dh,theta,mu,sigma2,invR,invRRes,nugget)
[Zhhats,ZhVars]=CompPosMeanVar(xs,Dh,theta,mu,sigma2,invR,invRRes,nugget);
ZhQuantiles=Zhhats+norminv(0.9)*ZhVars.^0.5;

end
%%
function [Zhhats,ZhVars,rT]=CompPosMeanVar(xs,Dh,theta,mu,sigma2,invR,invRRes,nugget)
rT=ComputeRmatrix(xs,Dh,theta);
rTinvR=rT*invR;
Zhhats=mu+rT*invRRes;
ZhVars=sigma2*(1+nugget-sum(rTinvR.*rT,2));
ZhVars=max(ZhVars,0);
end
%% 
function [Zhhats,ZhVars]=CompPosMeanVarGrid(mu,sigma2,invR,invRRes,nugget,rT)
rTinvR=rT*invR;
Zhhats=mu+rT*invRRes;
ZhVars=sigma2*(1+nugget-sum(rTinvR.*rT,2));
ZhVars=max(ZhVars,0);
end
%%
function MinusEIVals=CompMinusEI(xs,Dh,ZhhatxhatstarML,theta,mu,sigma2,invR,invRRes,nugget)
[Zhhats,ZhVars]=CompPosMeanVar(xs,Dh,theta,mu,sigma2,invR,invRRes,nugget);
sVals=ZhVars.^0.5;
DeltaVals=ZhhatxhatstarML-Zhhats;
DeltaValsoversVals=DeltaVals./sVals;
EIVals=DeltaVals.*normcdf(DeltaValsoversVals)+sVals.*normpdf(DeltaValsoversVals);
EIVals( (ZhVars==0) | (min(pdist2(xs,Dh),[],2) < 10^(-3)) )=0;
MinusEIVals=-EIVals;
end
%%
function MinusEIVals=CompMinusEIGrid(xs,Dh,ZhhatxhatstarML,Zhhats,ZhVars)
sVals=ZhVars.^0.5;
DeltaVals=ZhhatxhatstarML-Zhhats;
DeltaValsoversVals=DeltaVals./sVals;
EIVals=DeltaVals.*normcdf(DeltaValsoversVals)+sVals.*normpdf(DeltaValsoversVals);
EIVals( (ZhVars==0) | (min(pdist2(xs,Dh),[],2) < 10^(-3)) )=0;
MinusEIVals=-EIVals;
end