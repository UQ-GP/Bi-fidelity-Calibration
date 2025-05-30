function [RecordTable,RecordData,RunTime]=CalibrationAGP(DataInput,ID_BC_or_SR,ZMLFSSEorZLFSSE,AccuracyLevel)
% Implements four bi-fidelity calibration methods, i.e., the MBC-AGP (activated by specifying ID_BC_or_SR=1 and ZMLFSSEorZLFSSE=1 as inputs), BC-AGP (activated by specifying ID_BC_or_SR=1 and ZMLFSSEorZLFSSE=0 as inputs), MID-AGP (activated by specifying ID_BC_or_SR=0 and ZMLFSSEorZLFSSE=1 as inputs), and SR-AGP (activated by specifying ID_BC_or_SR=2 and ZMLFSSEorZLFSSE=0 as inputs) methods.
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
if ZMLFSSEorZLFSSE==1
    N=numel(w);
    [~,idxinDl,idxinDh]=intersect(Dl,Dh,'rows','stable');
    SameYl=Yl(idxinDl,:);
    SameYh=Yh(idxinDh,:);
    OnesVec=ones(numel(idxinDh),1);
    ahati_bhati=zeros(2,N);
    for kd=1:N
        Ylkd=SameYl(:,kd);
        Yhkd=SameYh(:,kd);
        ModelMatrix=[OnesVec,Ylkd];
        lastwarn('');
        ahati_bhati(:,kd)=regress(Yhkd,ModelMatrix);
        [warnMsg,~]=lastwarn;
        if contains(warnMsg,'X is rank deficient to within machine precision.')
            return
        end
    end
    YlModified=ahati_bhati(1,:)+Yl.*ahati_bhati(2,:);
    SlorSlpVec=sum((YlModified-w).^2,2);
    ShVec=sum((Yh-w).^2,2);
elseif ZMLFSSEorZLFSSE==0
    ShVec=sum((Yh-w).^2,2);
    SlorSlpVec=sum((Yl-w).^2,2);
end

SVec=[SlorSlpVec;ShVec];
Budget=Budget-(nl*1+nh*CostRatio);

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
        [sigma2l,thetal,rho,gamma2,sigma2h,thetah,phi,mu,invVprime,invVprimeRes,condVprime,MinM2LogLikelihood]=AGPFit(Dl,SlorSlpVec,Dh,ShVec,ID_BC_or_SR,nugget,AccuracyLevel);
        ZFit=0;
        MinM2LogLikelihoods(n,:)=MinM2LogLikelihood; sigma2ls(n,:)=sigma2l; thetals(n,:)=thetal; rhos(n,:)=rho;
        gamma2s(n,:)=gamma2;                         sigma2hs(n,:)=sigma2h; thetahs(n,:)=thetah;
        phis(n,:)=phi;                               mus(n,:)=mu;
    else
        [invVprime,invVprimeRes,condVprime]=M2LogLikelihoodSame(Dl,SlorSlpVec,Dh,ShVec,ID_BC_or_SR,thetal,thetah,rho,gamma2,phi,mu,nl,nh,nugget);
        
    end
    Dunique=unique([Dl;Dh],'rows','stable');
    [Zhhats1,ZhVars1,~,Corrs1]=CompPosMeanVarCorr(GridPoints,1,Dl,Dh,sigma2l,thetal,rho,gamma2,thetah,mu,invVprime,invVprimeRes,nugget);
    [Zhhats2,ZhVars2]=CompPosMeanVarCorr(Dunique,2,Dl,Dh,sigma2l,thetal,rho,gamma2,thetah,mu,invVprime,invVprimeRes,nugget);
    [Zhhats3,ZhVars3,~,Corrs3]=CompPosMeanVarCorr(HistoryxhatstarMLs,1,Dl,Dh,sigma2l,thetal,rho,gamma2,thetah,mu,invVprime,invVprimeRes,nugget);
    Zhhats=[Zhhats1;Zhhats2;Zhhats3]; ZhVars=[ZhVars1;ZhVars2;ZhVars3];
    ZhQuantileVals=Zhhats+norminv(0.9)*ZhVars.^0.5;
    [~,SortedZhQuantileValsIndices]=sort(ZhQuantileVals);
    options=optimoptions('patternsearch','Display','off'); 

    xhatstarMLCandidatePoints=[GridPoints;Dunique;HistoryxhatstarMLs];
    ZhQuantileFun=@(x) CompZhQuantile(x,2,Dl,Dh,sigma2l,thetal,rho,gamma2,thetah,mu,invVprime,invVprimeRes,nugget);
    XBestTry=zeros(NoS,d);
    fBestTry=zeros(NoS,1); SortedxhatstarMLCandidatePoints=xhatstarMLCandidatePoints(SortedZhQuantileValsIndices,:);
    parfor id=1:NoS
        [XBestTry(id,:),fBestTry(id,1)]=patternsearch(ZhQuantileFun,SortedxhatstarMLCandidatePoints(id,:),[],[],[],[],lb,ub,[],options);
    end
    [MinZhQuantile,minidx]=min(fBestTry);
    xhatstarML=XBestTry(minidx,:);
    
    xhatstarMLs(n,:)=xhatstarML;
    Shminhats(n,:)=TransformData_inv(MinZhQuantile,phi,ID_BC_or_SR);
    condVprimes(n,:)=condVprime;
    %Evaluates the true HF SSE at xhat^*_{ML}, i.e., S_h(xhat^*_{ML}).
    yhxhatstarMLs(n,:)=Simulator(xhatstarML,2,Case);
    ShxhatstarMLs(n,:)=sum((yhxhatstarMLs(n,:)-w).^2);
    
    disp(['Based on ' num2str(n) ' design points, the estimate xhat^*_{ML} of S_h''s minimizer and S_h(xhat^*_{ML}) are [' num2str(xhatstarMLs(n,:),' %.3f') '] and ' num2str(ShxhatstarMLs(n,:)) ' respectively.'])
    
    if Budget<1
        break
    end
    HistoryxhatstarMLs=[HistoryxhatstarMLs;xhatstarML];
    [ZhhatxhatstarML,ZhVars4,~,Corrs4]=CompPosMeanVarCorr(xhatstarML,1,Dl,Dh,sigma2l,thetal,rho,gamma2,thetah,mu,invVprime,invVprimeRes,nugget);
    ZhhatsAF=[Zhhats1;Zhhats3;ZhhatxhatstarML]; ZhVarsAF=[ZhVars1;ZhVars3;ZhVars4]; CorrsAF=[Corrs1;Corrs3;Corrs4];

    MaxAFCandidatePoints=[GridPoints;HistoryxhatstarMLs];
    %Adds a follow-up design point by maximizing the AEI AF.
    [NextPoint,NextLevel,MaxAFVal]=FindNextDesignPoint(Budget,CostRatio,Dl,Dh,sigma2l,thetal,rho,gamma2,thetah,mu,invVprime,invVprimeRes,MaxAFCandidatePoints,nugget,ZhhatxhatstarML,ZhhatsAF,ZhVarsAF,CorrsAF,AccuracyLevel);
    MaxAFVals(n+1,:)=MaxAFVal;
    disp(['Remaining budget=' num2str(Budget) '. Next, the ' num2str(n+1) '-th run will be made at the point [' num2str(NextPoint,' %1.3f') '] at fidelity level ' num2str(NextLevel) '.'])
    n=n+1;
    if NextLevel==1
        nl=nl+1;
        Dl(nl,:)=NextPoint;
        Yl(nl,:)=Simulator(NextPoint,1,Case);
        if ZMLFSSEorZLFSSE==1
            YlModified=[YlModified; ahati_bhati(1,:)+Yl(nl,:).*ahati_bhati(2,:)];
            SlorSlpVec=[SlorSlpVec; sum((YlModified(end,:)-w).^2,2)];
        elseif ZMLFSSEorZLFSSE==0
            SlorSlpVec=[SlorSlpVec; sum((Yl(end,:)-w).^2,2)];
        end
        Budget=Budget-1;
        SVec(n,1)=SlorSlpVec(nl,:);
        
    else
        nh=nh+1;
        Dh(nh,:)=NextPoint;
        Yh(nh,:)=Simulator(NextPoint,2,Case);
        ShVec=[ShVec; sum((Yh(end,:)-w).^2,2)];
        Budget=Budget-CostRatio;
        SVec(n,1)=ShVec(nh,:);
    end
    
    D(n,:)=NextPoint;
    Level(n,:)=NextLevel;
end
MinM2LogLikelihoods(n,:)=MinM2LogLikelihood; sigma2ls(n,:)=sigma2l; thetals(n,:)=thetal; rhos(n,:)=rho;
gamma2s(n,:)=gamma2;                         sigma2hs(n,:)=sigma2h; thetahs(n,:)=thetah;
phis(n,:)=phi;                               mus(n,:)=mu;    

%Stores design points, the simulator outputs at design points, and other information.
RecordData.Dl=Dl;
RecordData.Yl=Yl;
RecordData.SlorSlpVec=SlorSlpVec;
RecordData.Dh=Dh;
RecordData.Yh=Yh;
RecordData.ShVec=ShVec;
RecordData.xstar=xstar;
RecordData.w=w;
RecordData.CostRatio=CostRatio;
RecordData.Budget=Budget;
RecordData.yhxhatstarMLs=yhxhatstarMLs;

%Stores GP emulator parameter estimates, values of xhat^*_{ML} (rows of xhatstarMLs), values of S_h(xhat^*_{ML}) (elements of ShxhatstarMLs), and other important information with a table.
RecordTable=table(D,Level,SVec,MaxAFVals,xhatstarMLs,Shminhats,ShxhatstarMLs,sigma2ls,thetals,rhos,gamma2s,sigma2hs,thetahs,phis,mus,condVprimes,MinM2LogLikelihoods);
RunTime=toc;
end
%%
function [sigma2l,thetal,rho,gamma2,sigma2h,thetah,phi,mu,invVprime,invVprimeRes,condVprime,MinM2LogLikelihood]=AGPFit(Dl,SlorSlpVec,Dh,ShVec,ID_BC_or_SR,nugget,AccuracyLevel)
[nl,d]=size(Dl);
nh=size(Dh,1);
if(AccuracyLevel==1)
HNoS=90;
else
HNoS=100;
end
if ID_BC_or_SR==1
    lb=[(0.25)*ones(1,d) (0.25)*ones(1,d) 0.5 -6 -0.5];
    ub=[(15)*ones(1,d) (15)*ones(1,d) 1.5 0 1];
    Ranges=ub-lb;
    M2LogLFun=@(Par) M2LogLikelihood(Dl,SlorSlpVec,Dh,ShVec,ID_BC_or_SR,Par(1:d),Par((d+1):(2*d)),Par(2*d+1),10^Par(2*d+2),Par(2*d+3),nl,nh,nugget);
    npar=numel(lb);
    Sobolset=sobolset(npar,'Skip',1e3,'Leap',1e2);
    if(AccuracyLevel==1)
    StandardPoints=net(Sobolset,7000*npar);
    else
    StandardPoints=net(Sobolset,8000*npar);    
    end
    CandidatePoints=lb+Ranges.*StandardPoints;
    
elseif ID_BC_or_SR==0 || ID_BC_or_SR==2
    lb=[(0.25)*ones(1,d) (0.25)*ones(1,d) 0.5 -6];
    ub=[(15)*ones(1,d) (15)*ones(1,d) 1.5 0];
    Ranges=ub-lb;    
    M2LogLFun=@(Par) M2LogLikelihood(Dl,SlorSlpVec,Dh,ShVec,ID_BC_or_SR,Par(1:d),Par((d+1):(2*d)),Par(2*d+1),10^Par(2*d+2),0,nl,nh,nugget);
    npar=numel(lb);
    Sobolset=sobolset(npar,'Skip',1e3,'Leap',1e2);   
    if(AccuracyLevel==1)
    StandardPoints=net(Sobolset,7000*npar);
    else
    StandardPoints=net(Sobolset,8000*npar);    
    end
    CandidatePoints=lb+Ranges.*StandardPoints; 
    
end

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

[MinM2LogLikelihood,sigma2l,thetal,rho,gamma2,sigma2h,thetah,phi,mu,invVprime,invVprimeRes,condVprime]=M2LogLFun(OptPar);

end
%%
function [M2LogLikelihoodVal,sigma2l,thetal,rho,gamma2,sigma2h,thetah,phi,mu,invVprime,invVprimeRes,condVprime]=M2LogLikelihood(Dl,SlorSlpVec,Dh,ShVec,ID_BC_or_SR,thetal,thetah,rho,gamma2,phi,nl,nh,nugget)

[Zphi,LogAbsJacobian]=TransformData([SlorSlpVec;ShVec],phi,ID_BC_or_SR);
VprimeDl_Dh=rho*ComputeRmatrix(Dl,Dh,thetal);
Vprime=[ ComputeRmatrix2(Dl,thetal,nugget) , VprimeDl_Dh;
         VprimeDl_Dh' ,                      rho^2*ComputeRmatrix(Dh,Dh,thetal)+gamma2*ComputeRmatrix2(Dh,thetah,nugget) ];

F=[ ones(nl,1) ,     zeros(nl,1);
    rho*ones(nh,1) , ones(nh,1) ];
[invVprime,logdetVprime,condVprime]=invandlogdet(Vprime);

FTinvVprime=F'*invVprime;
invFTinvVprimeF=invandlogdet(FTinvVprime*F);

mu=invFTinvVprimeF*FTinvVprime*Zphi;
Res=Zphi-F*mu;

invVprimeRes=invVprime*Res;
sigma2l=(Res'*invVprimeRes)/(nl+nh);

M2LogLikelihoodVal=logdetVprime+(nl+nh)*log(sigma2l)-2*LogAbsJacobian;

if any(~isfinite(mu),'all') || sigma2l<=0 || ~isfinite(sigma2l)
    M2LogLikelihoodVal=Inf;sigma2l=[];thetal=[];rho=[];gamma2=[];sigma2h=[];thetah=[];phi=[];mu=[];invVprime=[];invVprimeRes=[];condVprime=[];
    return
end
sigma2h=gamma2*sigma2l;

end
%%
function [invVprime,invVprimeRes,condVprime]=M2LogLikelihoodSame(Dl,SlorSlpVec,Dh,ShVec,ID_BC_or_SR,thetal,thetah,rho,gamma2,phi,mu,nl,nh,nugget)

[Zphi]=TransformData([SlorSlpVec;ShVec],phi,ID_BC_or_SR);
VprimeDl_Dh=rho*ComputeRmatrix(Dl,Dh,thetal);
Vprime=[ ComputeRmatrix2(Dl,thetal,nugget) , VprimeDl_Dh;
         VprimeDl_Dh' ,                      rho^2*ComputeRmatrix(Dh,Dh,thetal)+gamma2*ComputeRmatrix2(Dh,thetah,nugget) ];

F=[ ones(nl,1) ,     zeros(nl,1);
    rho*ones(nh,1) , ones(nh,1) ];
[invVprime,~,condVprime]=invandlogdet(Vprime);
Res=Zphi-F*mu;
invVprimeRes=invVprime*Res;

if any(~isfinite(invVprimeRes),'all')
    invVprime=[];invVprimeRes=[];condVprime=[];
    return
end

end
%%
function ZhQuantiles=CompZhQuantile(xs,Choice,Dl,Dh,sigma2l,thetal,rho,gamma2,thetah,mu,invVprime,invVprimeRes,nugget)
[Zhhats,ZhVars]=CompPosMeanVarCorr(xs,Choice,Dl,Dh,sigma2l,thetal,rho,gamma2,thetah,mu,invVprime,invVprimeRes,nugget);
ZhQuantiles=Zhhats+norminv(0.9)*ZhVars.^0.5;
end

%%
function [Zhhats,ZhVars,ZlVars,Corrs,CrossCovs]=CompPosMeanVarCorr(xs,Choice,Dl,Dh,sigma2l,thetal,rho,gamma2,thetah,mu,invVprime,invVprimeRes,nugget)

RlT1=ComputeRmatrix(xs,Dl,thetal);
RlT2=ComputeRmatrix(xs,Dh,thetal);
RhT=ComputeRmatrix(xs,Dh,thetah);
    
if Choice==1
    rlToverSigma2l=[RlT1,rho*RlT2];
    rlToverSigma2linvVprime=rlToverSigma2l*invVprime;
    ZlVars=sigma2l*(1+nugget-sum(rlToverSigma2linvVprime.*rlToverSigma2l,2));
    ZlVars=max(ZlVars,0);
    
    fh=[rho,1];
    rhToverSigma2l=[rho*RlT1,rho^2*RlT2+gamma2*RhT];
    rhToverSigma2linvVprime=rhToverSigma2l*invVprime;
    Zhhats=fh*mu+rhToverSigma2l*invVprimeRes;
    ZhVars=sigma2l*(rho^2+gamma2*(1+nugget)-sum(rhToverSigma2linvVprime.*rhToverSigma2l,2));
    ZhVars=max(ZhVars,0);
        
    CrossCovs=sigma2l*(rho-sum(rlToverSigma2linvVprime.*rhToverSigma2l,2));
    Corrs=( CrossCovs ./ ((ZlVars).^0.5) ) ./ ((ZhVars).^0.5);
    Corrs( (ZhVars==0) | (ZlVars==0) )=0;
    Corrs=max(min(Corrs,1),-1);
    
elseif Choice==2
    
    fh=[rho,1];
    rhToverSigma2l=[rho*RlT1,rho^2*RlT2+gamma2*RhT];
    rhToverSigma2linvVprime=rhToverSigma2l*invVprime;
    Zhhats=fh*mu+rhToverSigma2l*invVprimeRes;
    ZhVars=sigma2l*(rho^2+gamma2*(1+nugget)-sum(rhToverSigma2linvVprime.*rhToverSigma2l,2));     
    ZhVars=max(ZhVars,0);
end

end
%%
function [NextPoint,NextLevel,MaxAFVal]=FindNextDesignPoint(Budget,CostRatio,Dl,Dh,sigma2l,thetal,rho,gamma2,thetah,mu,invVprime,invVprimeRes,MaxAFCandidatePoints,nugget,ZhhatxhatstarML,ZhhatsAF,ZhVarsAF,CorrsAF,AccuracyLevel)
d=size(Dl,2);
lb=0*ones(1,d); ub=1*ones(1,d);
if(AccuracyLevel==1) 
NoS=90;
else
NoS=100;
end
if Budget<CostRatio
    
    NextLevel=1;
    
    MinusAEIVals=CompMinusAEIGrid(MaxAFCandidatePoints,ZhhatxhatstarML,CostRatio,Dl,Dh,ZhhatsAF,ZhVarsAF,CorrsAF);
    [~,SortedMinusAEIValsIndices]=sort(MinusAEIVals);
    options=optimoptions('patternsearch','Display','off');
    MinusAEIFun=@(x) CompMinusAEI(x,ZhhatxhatstarML,CostRatio,NextLevel,Dl,Dh,sigma2l,thetal,rho,gamma2,thetah,mu,invVprime,invVprimeRes,nugget);
    XBestTry=zeros(NoS,d);
    fBestTry=zeros(NoS,1); SortedMaxAFCandidatePoints=MaxAFCandidatePoints(SortedMinusAEIValsIndices,:);
    parfor id=1:NoS
        [XBestTry(id,:),fBestTry(id,1)]=patternsearch(MinusAEIFun,SortedMaxAFCandidatePoints(id,:),[],[],[],[],lb,ub,[],options);
    end
    [fBest,minidx]=min(fBestTry);
    NextPoint=XBestTry(minidx,:);
    MaxAFVal=-fBest;
    
else
    
    Level1=1;

    [MinusAEIVals,MinusEIVals]=CompMinusAEIGrid(MaxAFCandidatePoints,ZhhatxhatstarML,CostRatio,Dl,Dh,ZhhatsAF,ZhVarsAF,CorrsAF);
    [~,SortedMinusAEIValsIndices]=sort(MinusAEIVals);
    options=optimoptions('patternsearch','Display','off');   
    MinusAEI1Fun=@(x) CompMinusAEI(x,ZhhatxhatstarML,CostRatio,Level1,Dl,Dh,sigma2l,thetal,rho,gamma2,thetah,mu,invVprime,invVprimeRes,nugget);
    XBestTry=zeros(NoS,d);
    fBestTry=zeros(NoS,1); SortedMaxAFCandidatePoints=MaxAFCandidatePoints(SortedMinusAEIValsIndices,:);
    parfor id=1:NoS
        [XBestTry(id,:),fBestTry(id,1)]=patternsearch(MinusAEI1Fun,SortedMaxAFCandidatePoints(id,:),[],[],[],[],lb,ub,[],options);
    end
    [MinMinusAEI1,minidx]=min(fBestTry);
    NextPoint1=XBestTry(minidx,:);

    Level2=2;

    [~,SortedMinusEIValsIndices]=sort(MinusEIVals);
    options=optimoptions('patternsearch','Display','off');
    MinusAEI2Fun=@(x) CompMinusAEI(x,ZhhatxhatstarML,CostRatio,Level2,Dl,Dh,sigma2l,thetal,rho,gamma2,thetah,mu,invVprime,invVprimeRes,nugget);
    XBestTry=zeros(NoS,d);
    fBestTry=zeros(NoS,1); SortedMaxAFCandidatePoints=MaxAFCandidatePoints(SortedMinusEIValsIndices,:);   
    parfor id=1:NoS
        [XBestTry(id,:),fBestTry(id,1)]=patternsearch(MinusAEI2Fun,SortedMaxAFCandidatePoints(id,:),[],[],[],[],lb,ub,[],options);
    end
    [MinMinusAEI2,minidx]=min(fBestTry);
    NextPoint2=XBestTry(minidx,:);
    
    if MinMinusAEI1<=MinMinusAEI2
        NextPoint=NextPoint1; NextLevel=1; MaxAFVal=-MinMinusAEI1;
    else
        NextPoint=NextPoint2; NextLevel=2; MaxAFVal=-MinMinusAEI2;
    end
end

end
%%
function [MinusAEIVals]=CompMinusAEI(xs,ZhhatxhatstarML,CostRatio,AEILevel,Dl,Dh,sigma2l,thetal,rho,gamma2,thetah,mu,invVprime,invVprimeRes,nugget)

if AEILevel==2
    Choice2=2;
    [Zhhats,ZhVars]=CompPosMeanVarCorr(xs,Choice2,Dl,Dh,sigma2l,thetal,rho,gamma2,thetah,mu,invVprime,invVprimeRes,nugget);
    DeltaVals=ZhhatxhatstarML-Zhhats;
    sVals=ZhVars.^0.5;
    DeltaValsoversVals=DeltaVals./sVals;
    EIs=DeltaVals.*normcdf(DeltaValsoversVals)+sVals.*normpdf(DeltaValsoversVals);
    AEIVals=EIs;
    AEIVals( (ZhVars==0) | (min(pdist2(xs,Dh),[],2) < 10^(-3)) ) =0;
elseif AEILevel==1
    Choice1=1;
    [Zhhats,ZhVars,~,Corrs]=CompPosMeanVarCorr(xs,Choice1,Dl,Dh,sigma2l,thetal,rho,gamma2,thetah,mu,invVprime,invVprimeRes,nugget);
    DeltaVals=ZhhatxhatstarML-Zhhats;
    sVals=ZhVars.^0.5;
    DeltaValsoversVals=DeltaVals./sVals;
    EIs=DeltaVals.*normcdf(DeltaValsoversVals)+sVals.*normpdf(DeltaValsoversVals);
    AEIVals=EIs.*Corrs*CostRatio;
    AEIVals( (ZhVars==0) |(min(pdist2(xs,Dl),[],2) < 10^(-3)) | (min(pdist2(xs,Dh),[],2) < 10^(-3)) | (Corrs < 0) )=0; 
end
MinusAEIVals=-AEIVals;
end
%%
function [MinusAEIVals,MinusEIVals]=CompMinusAEIGrid(xs,ZhhatxhatstarML,CostRatio,Dl,Dh,Zhhats,ZhVars,Corrs)

DeltaVals=ZhhatxhatstarML-Zhhats;
sVals=ZhVars.^0.5;
DeltaValsoversVals=DeltaVals./sVals;
EIs=DeltaVals.*normcdf(DeltaValsoversVals)+sVals.*normpdf(DeltaValsoversVals);
AEIs=EIs.*Corrs*CostRatio;
Ind=(ZhVars==0) | (min(pdist2(xs,Dh),[],2) < 10^(-3));
AEIs( Ind | (min(pdist2(xs,Dl),[],2) < 10^(-3)) | (Corrs < 0) )=0;
MinusAEIVals=-AEIs;
EIs( Ind )=0;
MinusEIVals=-EIs;
end