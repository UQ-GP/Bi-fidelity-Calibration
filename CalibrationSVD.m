function [RecordTable,RecordData,RunTime]=CalibrationSVD(DataInput,AccuracyLevel,t)
% Implements a single-fidelity calibration method, i.e., the SVD method.
tic
nugget=1e-6;
Dh=DataInput.Dh;
Yh=(DataInput.Yh)';%HF output at each design point is a column vector.
N=size(Yh,1);
xstar=DataInput.xstar;
w=(DataInput.w)';%Field data is a column vector.
CostRatio=DataInput.CostRatio;
Budget=DataInput.Budget;
Case=DataInput.Case;
[n,d]=size(Dh);
Level=2*ones(n,1);

ybarh=mean(Yh,2);
w2=w-ybarh;
Yhc=Yh-ybarh;
[U0,S0,~]=svd(Yhc);
S0=diag(S0).^2;
PossibleFractions=cumsum(S0)/sum(S0);
p=find(PossibleFractions>=t,1)
U=U0(:,1:p); U2=U.^2;
Ur=U0(:,(p+1):end);

omegamatrix=Yhc'*U;

epsilon=Yhc-U*(omegamatrix');
sigma2epsilon=var(epsilon(:),1);

cvec=((-w2)')*U;
dvec=((-w2)')*Ur;
dvecdvecT=dvec*(dvec');

S=(Yh-w).^2;
Smin=min(S,[],2); srSmin=Smin.^0.5;
SSEs=sum(S,1);
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
        [theta,sigma2,invR,invRRes,condR,MinM2LogLikelihood]=GPFitomega(Dh,omegamatrix,nugget,AccuracyLevel);
        thetas(n,:)=(theta(:))'; sigma2s(n,:)=sigma2;
        ps(n,:)=p;               MinM2LogLikelihoods(n,:)=MinM2LogLikelihood;
        ZFit=0; 
        
        [EShVals1,omegaMeans1,omegaVars1,rhT]=CompESh(GridPoints,Dh,omegamatrix,theta,sigma2,invR,invRRes,nugget,cvec,dvecdvecT,sigma2epsilon,N);

    else
        [invR,invRRes,condR]=M2LogLikelihoodSame(theta,Dh,omegamatrix,nugget);

        for jd=1:p
            rhT{jd}=[rhT{jd},ComputeRmatrix(GridPoints,Dh(n,:),theta(jd,:))];
        end
        [EShVals1,omegaMeans1,omegaVars1]=CompEShGrid(GridPoints,omegamatrix,sigma2,invR,invRRes,nugget,cvec,dvecdvecT,sigma2epsilon,N,rhT);
        Smin=min([Smin (Yh(:,n)-w).^2],[],2);
        srSmin=Smin.^0.5;
    end
    
    xhatstarMLCandidatePoints=[GridPoints;Dh;HistoryxhatstarMLs];
        
    EShFun=@(x) CompESh(x,Dh,omegamatrix,theta,sigma2,invR,invRRes,nugget,cvec,dvecdvecT,sigma2epsilon,N);
    EShVals2=EShFun(Dh);
    [EShVals3,omegaMeans3,omegaVars3]=EShFun(HistoryxhatstarMLs);
    EShVals=[EShVals1;EShVals2;EShVals3];
    [~,SortedEShValsIndices]=sort(EShVals);
    options=optimoptions('patternsearch','Display','off');
    
    XBestTry=zeros(NoS,d);
    fBestTry=zeros(NoS,1); SortedxhatstarMLCandidatePoints=xhatstarMLCandidatePoints(SortedEShValsIndices,:);
    parfor id=1:NoS
        [XBestTry(id,:),fBestTry(id,1)]=patternsearch(EShFun,SortedxhatstarMLCandidatePoints(id,:),[],[],[],[],lb,ub,[],options);
    end
    [Shminhats(n,:),minidx]=min(fBestTry);
    xhatstarML=XBestTry(minidx,:);
    
    %Evaluates the true HF SSE at xhat^*_{ML}, i.e., S_h(xhat^*_{ML}).
    xhatstarMLs(n,:)=xhatstarML;
    condRs(n,:)=condR;
    yhxhatstarMLs(n,:)=Simulator(xhatstarML,2,Case);
    ShxhatstarMLs(n,:)=sum(((yhxhatstarMLs(n,:)')-w).^2);

    disp(['Based on ' num2str(n) ' design points, the estimate xhat^*_{ML} of S_h''s minimizer and S_h(xhat^*_{ML}) are [' num2str(xhatstarMLs(n,:),' %.3f') '] and ' num2str(ShxhatstarMLs(n,:)) ' respectively.'])
    
    if Budget<CostRatio
        break
    end
    HistoryxhatstarMLs=[HistoryxhatstarMLs;xhatstarML];    
    MaxAFCandidatePoints=[GridPoints;HistoryxhatstarMLs];
    
    %Adds a follow-up design point by maximizing the AF.
    MinusSumEIsFun=@(x) CompMinusSumEIs(x,Dh,omegamatrix,theta,sigma2,invR,invRRes,nugget,w2,U,U2,Smin,srSmin,sigma2epsilon);
        
    MinusSumEIsVals1=CompMinusSumEIsGrid(MaxAFCandidatePoints(1:(end-1),:),Dh,w2,U,U2,Smin,srSmin,sigma2epsilon,[omegaMeans1;omegaMeans3],[omegaVars1;omegaVars3]);
    MinusSumEIsVals2=MinusSumEIsFun(xhatstarML);
    MinusSumEIsVals=[MinusSumEIsVals1;MinusSumEIsVals2];
    [~,SortedMinusSumEIsValsIndices]=sort(MinusSumEIsVals);
    options=optimoptions('patternsearch','Display','off');
    
    XBestTry=zeros(NoS,d);
    fBestTry=zeros(NoS,1); SortedMaxAFCandidatePoints=MaxAFCandidatePoints(SortedMinusSumEIsValsIndices,:);
    parfor id=1:NoS
        [XBestTry(id,:),fBestTry(id,1)]=patternsearch(MinusSumEIsFun,SortedMaxAFCandidatePoints(id,:),[],[],[],[],lb,ub,[],options);
    end
    [fBest,minidx]=min(fBestTry);
    NextPoint=XBestTry(minidx,:);
    MaxAFVal=-fBest; MaxAFVals(n+1,:)=MaxAFVal;
    disp(['Remaining budget=' num2str(Budget) '. Next, the ' num2str(n+1) '-th run will be made at the point [' num2str(NextPoint,' %1.3f') '].'])
    n=n+1;
    Dh(n,:)=NextPoint;
    Yh(:,n)=Simulator(NextPoint,2,Case)';
    Yhc(:,n)=Yh(:,n)-ybarh;
    
    omegamatrix(n,:)=(Yhc(:,n)')*U;
    
    SSEs(n)=sum((Yh(:,n)-w).^2);
    Level(n,:)=2;
    
    Budget=Budget-CostRatio;
end

thetas(n,:)=(theta(:))'; sigma2s(n,:)=sigma2; 
ps(n,:)=p;               MinM2LogLikelihoods(n,:)=MinM2LogLikelihood;

%Stores design points, the simulator outputs at design points, and other information.
SSEs=SSEs';
D=Dh;
RecordData.Dl=[];
RecordData.Yl=[];
RecordData.Dh=Dh;
RecordData.Yh=Yh;
RecordData.xstar=xstar;
RecordData.w=w;
RecordData.CostRatio=CostRatio;
RecordData.Budget=Budget;
RecordData.yhxhatstarMLs=yhxhatstarMLs;

%Stores GP emulator parameter estimates, values of xhat^*_{ML} (rows of xhatstarMLs), values of S_h(xhat^*_{ML}) (elements of ShxhatstarMLs), and other important information with a table.
RecordTable=table(D,Level,SSEs,MaxAFVals,xhatstarMLs,Shminhats,ShxhatstarMLs,ps,sigma2s,thetas,condRs,MinM2LogLikelihoods);
RunTime=toc;
end
%%
function [theta,sigma2,invR,invRRes,maxcondR,MinM2LogLikelihood]=GPFitomega(Dh,omegamatrix,nugget,AccuracyLevel)

d=size(Dh,2);
p=size(omegamatrix,2);
lb=(0.25)*ones(1,d);
ub=(15)*ones(1,d);
Ranges=ub-lb;
npar=numel(lb);
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
theta=zeros(p,d);
sigma2=zeros(1,p);
invR=cell(1,p);
invRRes=cell(1,p);
condR=zeros(1,p); 
MinM2LogLikelihood=0;
options=optimoptions('patternsearch','Display','off'); 
for jd=1:p
    omega0=omegamatrix(:,jd);
    M2LogLFun=@(Par) M2LogLikelihood0(Par,Dh,omega0,nugget);
    
    parfor id=1:NoCandidatePoints
        M2LogLVals(id,1)=M2LogLFun(CandidatePoints(id,:));
    end
    [~,SortedM2LogLValsIndices]=sort(M2LogLVals);
    
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
    [MinM2LogLikelihood0,sigma20,theta0,invR0,invRRes0,condR0]=M2LogLFun(OptPar);
    
    theta(jd,:)=theta0;
    sigma2(jd)=sigma20;
    invR{jd}=invR0;
    invRRes{jd}=invRRes0;
    condR(jd)=condR0;
    MinM2LogLikelihood=MinM2LogLikelihood+MinM2LogLikelihood0;
end
maxcondR=max(condR);
end
%%
function [M2LogLikelihood0Val,sigma20,theta0,invR0,invRRes0,condR0]=M2LogLikelihood0(theta0,Dh,omega0,nugget)
nh=size(omega0,1);
R0=ComputeRmatrix2(Dh,theta0,nugget);
[invR0,logdetR0,condR0]=invandlogdet(R0);
Res0=omega0;
invRRes0=invR0*Res0;
sigma20=(Res0')*invRRes0/nh;
M2LogLikelihood0Val=nh*log(sigma20)+logdetR0;
if sigma20<=0 || ~isfinite(sigma20)
    M2LogLikelihood0Val=Inf;sigma20=[];theta0=[];invR0=[];invRRes0=[];condR0=[];
    return
end
end
%%
function [invR,invRRes,maxcondR]=M2LogLikelihoodSame(theta,Dh,omegamatrix,nugget)
p=size(omegamatrix,2);
invR=cell(1,p);
invRRes=cell(1,p);
condR=zeros(1,p);
Flag=zeros(1,p)~=zeros(1,p);
parfor jd=1:p
    omega0=omegamatrix(:,jd);
    theta0=theta(jd,:);
    R0=ComputeRmatrix2(Dh,theta0,nugget);
    [invR0,~,condR0]=invandlogdet(R0);
    Res0=omega0;
    invRRes0=invR0*Res0;
    if any(~isfinite(invRRes0),'all')
        Flag(jd)=true;
    end    
    invRRes{jd}=invRRes0;
    invR{jd}=invR0;
    condR(jd)=condR0;
end
if(any(Flag))
    invR=[];invRRes=[];maxcondR=[];
    return
end
maxcondR=max(condR);
end
%%
function [EShs,omegaMeans,omegaVars,rhT]=CompESh(xs,Dh,omegamatrix,theta,sigma2,invR,invRRes,nugget,cvec,dvecdvecT,sigma2epsilon,N)
[omegaMeans,omegaVars,rhT]=CompPosMeanVar(xs,Dh,omegamatrix,theta,sigma2,invR,invRRes,nugget);
Part1=sum((omegaMeans+cvec).^2,2);
Part2=sum(omegaVars,2)+N*sigma2epsilon;
EShs=Part1+Part2+dvecdvecT;
end
%%
function [EShs,omegaMeans,omegaVars]=CompEShGrid(xs,omegamatrix,sigma2,invR,invRRes,nugget,cvec,dvecdvecT,sigma2epsilon,N,rhT)
[omegaMeans,omegaVars]=CompPosMeanVarGrid(xs,omegamatrix,sigma2,invR,invRRes,nugget,rhT);
Part1=sum((omegaMeans+cvec).^2,2);
Part2=sum(omegaVars,2)+N*sigma2epsilon;
EShs=Part1+Part2+dvecdvecT;
end
%%
function [omegaMeans,omegaVars,rhT]=CompPosMeanVar(xs,Dh,omegamatrix,theta,sigma2,invR,invRRes,nugget)
p=size(omegamatrix,2);
Noxs=size(xs,1);
omegaMeans=zeros(Noxs,p);
omegaVars=zeros(Noxs,p); 
rhT=cell(1,p);
parfor jd=1:p
    theta0=theta(jd,:);
    invR0=invR{jd};    
    invRRes0=invRRes{jd};
    rhT{jd}=ComputeRmatrix(xs,Dh,theta0);
    omegaMeans(:,jd)=rhT{jd}*invRRes0;
    omegaVars0=sigma2(jd)*(1+nugget-sum((rhT{jd}*invR0).*rhT{jd},2));
    omegaVars(:,jd)=max(omegaVars0,0);
end
end
%%
function [omegaMeans,omegaVars]=CompPosMeanVarGrid(xs,omegamatrix,sigma2,invR,invRRes,nugget,rhT)
p=size(omegamatrix,2);
Noxs=size(xs,1);
omegaMeans=zeros(Noxs,p);
omegaVars=zeros(Noxs,p);

parfor jd=1:p
    invR0=invR{jd};
    invRRes0=invRRes{jd};    
    omegaMeans(:,jd)=rhT{jd}*invRRes0;
    omegaVars0=sigma2(jd)*(1+nugget-sum((rhT{jd}*invR0).*rhT{jd},2));
    omegaVars(:,jd)=max(omegaVars0,0);
end
end
%%
function [MinusSumEIsVals]=CompMinusSumEIs(xs,Dh,omegamatrix,theta,sigma2,invR,invRRes,nugget,w2,U,U2,Smin,srSmin,sigma2epsilon)
[omegaMeans,omegaVars]=CompPosMeanVar(xs,Dh,omegamatrix,theta,sigma2,invR,invRRes,nugget);

Noxs=size(xs,1);
SumEIsVals=zeros(Noxs,1);
parfor id=1:Noxs
    yhMean_centered=U*(omegaMeans(id,:)');
    yhVar=U2*(omegaVars(id,:)')+sigma2epsilon;
    
    sryhVar=sqrt(yhVar);
    Res=w2-yhMean_centered;
    Qplus=(Res+srSmin)./sryhVar;
    Qminus=(Res-srSmin)./sryhVar;
    EIPart1=(Smin-Res.^2-yhVar).*(normcdf(Qplus)-normcdf(Qminus));
    EIPart2=(srSmin-Res).*normpdf(Qplus)+(srSmin+Res).*normpdf(Qminus);
    EIs=EIPart1+sryhVar.*EIPart2;
    I0=find(yhVar==0);
    EIs(I0)=max(Smin(I0)-Res(I0).^2,0);    
    SumEIsVals(id,:)=sum(EIs);
end
SumEIsVals( min(pdist2(xs,Dh),[],2) < 10^(-3) )=0;
MinusSumEIsVals=-SumEIsVals;
end
%%
function [MinusSumEIsVals]=CompMinusSumEIsGrid(xs,Dh,w2,U,U2,Smin,srSmin,sigma2epsilon,omegaMeans,omegaVars)

Noxs=size(xs,1);
SumEIsVals=zeros(Noxs,1);
parfor id=1:Noxs
    yhMean_centered=U*(omegaMeans(id,:)');    
    yhVar=U2*(omegaVars(id,:)')+sigma2epsilon;
    
    sryhVar=sqrt(yhVar);
    Res=w2-yhMean_centered;
    Qplus=(Res+srSmin)./sryhVar;
    Qminus=(Res-srSmin)./sryhVar;
    EIPart1=(Smin-Res.^2-yhVar).*(normcdf(Qplus)-normcdf(Qminus));
    EIPart2=(srSmin-Res).*normpdf(Qplus)+(srSmin+Res).*normpdf(Qminus);
    EIs=EIPart1+sryhVar.*EIPart2;
    I0=find(yhVar==0);
    EIs(I0)=max(Smin(I0)-Res(I0).^2,0);
    SumEIsVals(id,:)=sum(EIs);
end
SumEIsVals( min(pdist2(xs,Dh),[],2) < 10^(-3) )=0;
MinusSumEIsVals=-SumEIsVals;
end