function [RecordTable,RecordData,RunTime]=CalibrationSVDAGP(DataInput,AccuracyLevel,t)
% Implements a bi-fidelity calibration method, i.e., the SVD-AGP method.
tic
nugget=1e-6;
Dl=DataInput.Dl;
Yl=(DataInput.Yl)';%LF output at each design point is a column vector.
Dh=DataInput.Dh;
Yh=(DataInput.Yh)';%HF output at each design point is a column vector.
nh=size(Dh,1); [nl,d]=size(Dl);
xstar=DataInput.xstar;
w=(DataInput.w)';%Field data is a column vector.
CostRatio=DataInput.CostRatio;
Budget=DataInput.Budget;
Case=DataInput.Case;

ybarl=mean(Yl,2);
Ylc=Yl-ybarl;
[Ul0,Sl0,~]=svd(Ylc);
Sl0=diag(Sl0).^2;
PossibleFractions=cumsum(Sl0)/sum(Sl0);
pl=find(PossibleFractions>=t,1)
Ul=Ul0(:,1:pl); Ul2=Ul.^2;

Deltah=Yh-Yl(:,1:nh);
Deltabarh=mean(Deltah,2);
Deltahc=Deltah-Deltabarh;
[Uh0,Sh0,~]=svd(Deltahc);
Sh0=diag(Sh0).^2;
PossibleFractions=cumsum(Sh0)/sum(Sh0);
ph=find(PossibleFractions>=t,1)
Uh=Uh0(:,1:ph); Uh2=Uh.^2;

omegalmatrix=Ylc'*Ul;
omegahmatrix=Deltahc'*Uh;

epsilonl=Ylc-Ul*(omegalmatrix');
sigma2epsilonl=var(epsilonl(:),1);

epsilonh=Deltahc-Uh*(omegahmatrix');
sigma2epsilonh=var(epsilonh(:),1);

n=nl+nh;
D=[Dl;Dh];
Level=[ones(nl,1);2*ones(nh,1)];

S=(Yh-w).^2;
Smin=min(S,[],2); srSmin=Smin.^0.5;
SSEs=sum([(Yl-w).^2 S],1);
Budget=Budget-(nl*1+nh*CostRatio);
MinBudget=CostRatio+1;

DlnotDh=Dl(nh+1:end,:);
omegalnotDh=omegalmatrix(nh+1:end,:);
YlnotDh=Yl(:,nh+1:end);
YlcnotDh=Ylc(:,nh+1:end);

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
        [thetal,sigma2l,invRl,invRlRes,condRl,MinM2LogLikelihoodl]=GPFitomega(Dl,omegalmatrix,nugget,AccuracyLevel);
        [thetah,sigma2h,invRh,invRhRes,condRh,MinM2LogLikelihoodh]=GPFitomega(Dh,omegahmatrix,nugget,AccuracyLevel);
        thetals(n-1:n,:)=repmat((thetal(:))',2,1); sigma2ls(n-1:n,:)=[sigma2l;sigma2l];
        thetahs(n-1:n,:)=repmat((thetah(:))',2,1); sigma2hs(n-1:n,:)=[sigma2h;sigma2h];
        
        pls(n-1:n,:)=[pl;pl]; MinM2LogLikelihoodls(n-1:n,:)=[MinM2LogLikelihoodl;MinM2LogLikelihoodl];
        phs(n-1:n,:)=[ph;ph]; MinM2LogLikelihoodhs(n-1:n,:)=[MinM2LogLikelihoodh;MinM2LogLikelihoodh];
        ZFit=0;
    else
        [invRl,invRlRes,condRl]=M2LogLikelihoodSame(thetal,Dl,omegalmatrix,nugget);
        [invRh,invRhRes,condRh]=M2LogLikelihoodSame(thetah,Dh,omegahmatrix,nugget);
        Smin=min([Smin (Yh(:,nh)-w).^2],[],2);
        srSmin=Smin.^0.5;
    end
    
    xhatstarMLCandidatePoints=[GridPoints;Dl;HistoryxhatstarMLs]; %All rows of Dh are rows of Dl.
    
    EShFun=@(x) CompESh(x,Dh,omegahmatrix,thetah,sigma2h,invRh,invRhRes,Dl,omegalmatrix,thetal,sigma2l,invRl,invRlRes,nugget,Ul,Ul2,Uh,Uh2,ybarl,Deltabarh,w,sigma2epsilonl,sigma2epsilonh);
    [EShVals1,omegalMeans1,omegalVars1,omegahMeans1,omegahVars1]=EShFun(GridPoints);
    EShVals2=EShFun(Dl);
    [EShVals3,omegalMeans3,omegalVars3,omegahMeans3,omegahVars3]=EShFun(HistoryxhatstarMLs);
    
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
    condRlRhs(n-1:n,:)=[condRl,condRh;condRl,condRh];
    %Evaluates the true HF SSE at xhat^*_{ML}, i.e., S_h(xhat^*_{ML}).
    xhatstarMLs(n-1:n,:)=[xhatstarML;xhatstarML];
    yhxhatstarML=Simulator(xhatstarML,2,Case);
    ShxhatstarML=sum(((yhxhatstarML')-w).^2);
    yhxhatstarMLs(n-1:n,:)=[yhxhatstarML;yhxhatstarML];
    ShxhatstarMLs(n-1:n,:)=[ShxhatstarML;ShxhatstarML];
    
    disp(['Based on ' num2str(n) ' design points, the estimate xhat^*_{ML} of S_h''s minimizer and S_h(xhat^*_{ML}) are [' num2str(xhatstarMLs(n,:),' %.3f') '] and ' num2str(ShxhatstarMLs(n,:)) ' respectively.'])
 
    if Budget<MinBudget
        break
    end
    HistoryxhatstarMLs=[HistoryxhatstarMLs;xhatstarML];
    MaxAFCandidatePoints=[GridPoints;HistoryxhatstarMLs];
    
    %Adds a follow-up design point by maximizing the AF.
    MinusSumEIsFun=@(x) CompMinusSumEIs(x,Dh,omegahmatrix,thetah,sigma2h,invRh,invRhRes,Dl,omegalmatrix,thetal,sigma2l,invRl,invRlRes,nugget,Ul,Ul2,Uh,Uh2,ybarl,Deltabarh,Smin,srSmin,w,sigma2epsilonl,sigma2epsilonh);
    
    MinusSumEIsVals1=CompMinusSumEIsGrid(GridPoints,Dl,Ul,Ul2,Uh,Uh2,ybarl,Deltabarh,Smin,srSmin,w,sigma2epsilonl,sigma2epsilonh,omegalMeans1,omegalVars1,omegahMeans1,omegahVars1);
    MinusSumEIsVals3=CompMinusSumEIsGrid(HistoryxhatstarMLs(1:(end-1),:),Dl,Ul,Ul2,Uh,Uh2,ybarl,Deltabarh,Smin,srSmin,w,sigma2epsilonl,sigma2epsilonh,omegalMeans3,omegalVars3,omegahMeans3,omegahVars3);
    MinusSumEIsVals4=MinusSumEIsFun(xhatstarML);
    MinusSumEIsVals=[MinusSumEIsVals1;MinusSumEIsVals3;MinusSumEIsVals4];
    [~,SortedMinusSumEIsValsIndices]=sort(MinusSumEIsVals);
    options=optimoptions('patternsearch','Display','off');
    
    XBestTry=zeros(NoS,d);
    fBestTry=zeros(NoS,1); SortedMaxAFCandidatePoints=MaxAFCandidatePoints(SortedMinusSumEIsValsIndices,:);
    parfor id=1:NoS
        [XBestTry(id,:),fBestTry(id,1)]=patternsearch(MinusSumEIsFun,SortedMaxAFCandidatePoints(id,:),[],[],[],[],lb,ub,[],options);
    end
    [fBest,minidx]=min(fBestTry);
    NextPoint=XBestTry(minidx,:);
    MaxAFVal=-fBest; MaxAFVals((n+1):(n+2),:)=[MaxAFVal;MaxAFVal];   
    disp(['Remaining budget=' num2str(Budget) '. Next, the ' num2str(n+1) '-th and ' num2str(n+2) '-th runs will be made at the point [' num2str(NextPoint,' %1.3f') '].'])
    Dh(nh+1,:)=NextPoint;
    Yh(:,nh+1)=(Simulator(NextPoint,2,Case))';
    
    Dl=[Dh;DlnotDh];
    Yl=[Yl(:,1:nh),(Simulator(NextPoint,1,Case))',YlnotDh];
    Ylc=[Ylc(:,1:nh),Yl(:,nh+1)-ybarl,YlcnotDh]; 
    omegal_new=(Ylc(:,nh+1)')*Ul;
    omegalmatrix=[omegalmatrix(1:nh,:);omegal_new;omegalnotDh];
    
    Deltah(:,nh+1)=Yh(:,nh+1)-Yl(:,nh+1);
    Deltahc(:,nh+1)=Deltah(:,nh+1)-Deltabarh;
    omegahmatrix(nh+1,:)=(Deltahc(:,nh+1)')*Uh;
    
    D=[D;NextPoint;NextPoint];
    Level=[Level;2;1];
    
    SSEs=[SSEs sum((Yh(:,end)-w).^2) sum((Yl(:,nh+1)-w).^2)];
    
    Budget=Budget-CostRatio-1;
    nh=nh+1;
    nl=nl+1;
    n=n+2;
end
thetals(n-1:n,:)=repmat((thetal(:))',2,1); sigma2ls(n-1:n,:)=[sigma2l;sigma2l];
thetahs(n-1:n,:)=repmat((thetah(:))',2,1); sigma2hs(n-1:n,:)=[sigma2h;sigma2h];

pls(n-1:n,:)=[pl;pl]; MinM2LogLikelihoodls(n-1:n,:)=[MinM2LogLikelihoodl;MinM2LogLikelihoodl];
phs(n-1:n,:)=[ph;ph]; MinM2LogLikelihoodhs(n-1:n,:)=[MinM2LogLikelihoodh;MinM2LogLikelihoodh];

%Stores design points, the simulator outputs at design points, and other information.
SSEs=SSEs';
RecordData.Dl=Dl;
RecordData.Yl=Yl;
RecordData.Dh=Dh;
RecordData.Yh=Yh;
RecordData.xstar=xstar;
RecordData.w=w;
RecordData.CostRatio=CostRatio;
RecordData.Budget=Budget;
RecordData.yhxhatstarMLs=yhxhatstarMLs;

%Stores GP emulator parameter estimates, values of xhat^*_{ML} (rows of xhatstarMLs), values of S_h(xhat^*_{ML}) (elements of ShxhatstarMLs), and other important information with a table.
RecordTable=table(D,Level,SSEs,MaxAFVals,xhatstarMLs,Shminhats,ShxhatstarMLs,pls,phs,thetals,sigma2ls,thetahs,sigma2hs,condRlRhs,MinM2LogLikelihoodhs,MinM2LogLikelihoodls);
RunTime=toc;
end
%%
function [theta,sigma2,invR,invRRes,maxcondR,MinM2LogLikelihood]=GPFitomega(D,omegamatrix,nugget,AccuracyLevel)

d=size(D,2);
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
    M2LogLFun=@(Par) M2LogLikelihood0(Par,D,omega0,nugget);
    
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
function [M2LogLikelihood0Val,sigma20,theta0,invR0,invRRes0,condR0]=M2LogLikelihood0(theta0,D,omega0,nugget)
n=size(omega0,1);
R0=ComputeRmatrix2(D,theta0,nugget); 
[invR0,logdetR0,condR0]=invandlogdet(R0);
Res0=omega0;
invRRes0=invR0*Res0;
sigma20=(Res0')*invRRes0/n;
M2LogLikelihood0Val=n*log(sigma20)+logdetR0;
if sigma20<=0 || ~isfinite(sigma20)
    M2LogLikelihood0Val=Inf;sigma20=[];theta0=[];invR0=[];invRRes0=[];condR0=[];
    return
end
end
%%
function [invR,invRRes,maxcondR]=M2LogLikelihoodSame(theta,D,omegamatrix,nugget)
p=size(omegamatrix,2);
invR=cell(1,p);
invRRes=cell(1,p);
condR=zeros(1,p);
Flag=zeros(1,p)~=zeros(1,p);
parfor jd=1:p
    omega0=omegamatrix(:,jd);
    theta0=theta(jd,:);
    R0=ComputeRmatrix2(D,theta0,nugget); 
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
function [EShs,omegalMeans,omegalVars,omegahMeans,omegahVars,rlT,rhT]=CompESh(xs,Dh,omegahmatrix,thetah,sigma2h,invRh,invRhRes,Dl,omegalmatrix,thetal,sigma2l,invRl,invRlRes,nugget,Ul,Ul2,Uh,Uh2,ybarl,Deltabarh,w,sigma2epsilonl,sigma2epsilonh)
[omegalMeans,omegalVars,rlT]=CompPosMeanVar(xs,Dl,omegalmatrix,thetal,sigma2l,invRl,invRlRes,nugget);
[omegahMeans,omegahVars,rhT]=CompPosMeanVar(xs,Dh,omegahmatrix,thetah,sigma2h,invRh,invRhRes,nugget);
Noxs=size(xs,1); EShs=zeros(Noxs,1); 
parfor id=1:Noxs
    ylMean=Ul*(omegalMeans(id,:)')+ybarl;
    DeltahMean=Uh*(omegahMeans(id,:)')+Deltabarh;
    yhMean=ylMean+DeltahMean;
    yhVar=Ul2*(omegalVars(id,:)')+sigma2epsilonl+Uh2*(omegahVars(id,:)')+sigma2epsilonh;
    
    Part1=sum((yhMean-w).^2,1);
    Part2=sum(yhVar,1);
    EShs(id)=Part1+Part2;
end
end
%%
function [omegaMeans,omegaVars,rT]=CompPosMeanVar(xs,D,omega,theta,sigma2,invR,invRRes,nugget)
p=size(omega,2);
Noxs=size(xs,1);
omegaMeans=zeros(Noxs,p);
omegaVars=zeros(Noxs,p); 
rT=cell(1,p);
parfor jd=1:p
    theta0=theta(jd,:);
    invR0=invR{jd};    
    invRRes0=invRRes{jd};
    rT{jd}=ComputeRmatrix(xs,D,theta0);
    omegaMeans(:,jd)=rT{jd}*invRRes0;
    omegaVars0=sigma2(jd)*(1+nugget-sum((rT{jd}*invR0).*rT{jd},2));
    omegaVars(:,jd)=max(omegaVars0,0);
end
end
%%
function MinusSumEIsVals=CompMinusSumEIs(xs,Dh,omegahmatrix,thetah,sigma2h,invRh,invRhRes,Dl,omegalmatrix,thetal,sigma2l,invRl,invRlRes,nugget,Ul,Ul2,Uh,Uh2,ybarl,Deltabarh,Smin,srSmin,w,sigma2epsilonl,sigma2epsilonh)
[omegalMeans,omegalVars]=CompPosMeanVar(xs,Dl,omegalmatrix,thetal,sigma2l,invRl,invRlRes,nugget); 
[omegahMeans,omegahVars]=CompPosMeanVar(xs,Dh,omegahmatrix,thetah,sigma2h,invRh,invRhRes,nugget);

Noxs=size(xs,1);
SumEIsVals=zeros(Noxs,1);
parfor id=1:Noxs
    
    ylMean=Ul*(omegalMeans(id,:)')+ybarl;
    DeltahMean=Uh*(omegahMeans(id,:)')+Deltabarh;
    yhMean=ylMean+DeltahMean;
    yhVar=(Ul2)*(omegalVars(id,:)')+sigma2epsilonl+(Uh2)*(omegahVars(id,:)')+sigma2epsilonh;
    
    sryhVar=sqrt(yhVar);
    Res=w-yhMean;
    Qplus=(Res+srSmin)./sryhVar;
    Qminus=(Res-srSmin)./sryhVar;
    EIPart1=(Smin-Res.^2-yhVar).*(normcdf(Qplus)-normcdf(Qminus));
    EIPart2=(srSmin-Res).*normpdf(Qplus)+(srSmin+Res).*normpdf(Qminus);
    EIs=EIPart1+sryhVar.*EIPart2;
    I0=find(yhVar==0);
    EIs(I0)=max(Smin(I0)-Res(I0).^2,0);      
    SumEIsVals(id,:)=sum(EIs);
end
SumEIsVals( min(pdist2(xs,Dl),[],2) < 10^(-3) )=0;
MinusSumEIsVals=-SumEIsVals;
end
%%
function MinusSumEIsVals=CompMinusSumEIsGrid(xs,Dl,Ul,Ul2,Uh,Uh2,ybarl,Deltabarh,Smin,srSmin,w,sigma2epsilonl,sigma2epsilonh,omegalMeans,omegalVars,omegahMeans,omegahVars)

Noxs=size(xs,1);
SumEIsVals=zeros(Noxs,1);
parfor id=1:Noxs
    
    ylMean=Ul*(omegalMeans(id,:)')+ybarl;
    DeltahMean=Uh*(omegahMeans(id,:)')+Deltabarh;
    yhMean=ylMean+DeltahMean;
    yhVar=(Ul2)*(omegalVars(id,:)')+sigma2epsilonl+(Uh2)*(omegahVars(id,:)')+sigma2epsilonh;
    
    sryhVar=sqrt(yhVar);
    Res=w-yhMean;
    Qplus=(Res+srSmin)./sryhVar;
    Qminus=(Res-srSmin)./sryhVar;
    EIPart1=(Smin-Res.^2-yhVar).*(normcdf(Qplus)-normcdf(Qminus));
    EIPart2=(srSmin-Res).*normpdf(Qplus)+(srSmin+Res).*normpdf(Qminus);
    EIs=EIPart1+sryhVar.*EIPart2;
    I0=find(yhVar==0);
    EIs(I0)=max(Smin(I0)-Res(I0).^2,0);       
    SumEIsVals(id,:)=sum(EIs);
end
SumEIsVals( min(pdist2(xs,Dl),[],2) < 10^(-3) )=0;
MinusSumEIsVals=-SumEIsVals;
end