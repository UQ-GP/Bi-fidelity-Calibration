function [RecordTable,RecordData]=CalibrationSRGP(DataInput,Val)
% Single fidelity calibration method : SR-GP
%           Resph: is a vector of HF response modeled as a GP
%           Xhats: is MLE of the calibration parameter vector at all iterations
%           Resphminhats: is posterior means/medians of the Resph at Xhats at all iterations
%           SSETrue_Xhats: is the true values of SSE evaluated at Xhats at all iterations
% Referennces:
% Ranjan, P., Thomas, M., Teismann, H., & Mukhoti, S. (2016). Inverse problem for a time-series valued computer simulator via scalarization.
% Open Journal of Statistics, 6(3), 528-544.

nugget=1e-6;
Dh=DataInput.Dh;
Yh=DataInput.Yh;
XTrue=DataInput.XTrue;
PhysData=DataInput.PhysData;
RatioCost=DataInput.RatioCost;
Budget=DataInput.Budget;
Case=DataInput.Case;
[n,Dim]=size(Dh);
Level=2*ones(n,1) ;


Resph=(mean((Yh-PhysData).^2,2)).^0.5;
Budget =Budget-(n*RatioCost);

if Dim==2
    if(Val==1)
    nlevel=2501;
    else
    nlevel=3001;    
    end
elseif Dim==3
    if(Val==1)
    nlevel=201;
    else
    nlevel=226;
    end
end
AFPoints=(fullfact(nlevel*ones(1,Dim))-1)/(nlevel-1);

lb=0*ones(1,Dim);ub=1*ones(1,Dim);
options=optimoptions('patternsearch','disp','off');

HistoryXhats=[];
AFs(n,:)=0;
%Bayesian optimization
ZFit=1;
if(Val==1)
NoS=90;
else
NoS=100;    
end
while (1)
    if ZFit==1
        [Mu,Sigma,Theta,invR,invRRes,condR,~,PDF]=GPFit(Dh,Resph,nugget,Val);
        ZFit=0;
        Sigmas(n,:)=Sigma;    Thetas(n,:)=Theta;    Mus(n,:)=Mu;
        PDFs(n,:)=PDF;

        [ZhPreds1,ZhCovs1,rT]=Fun_GPPrediction(AFPoints,Dh,Resph,Theta,Mu,Sigma,invR,invRRes,nugget);
        
    else
        [invR,invRRes,condR,FT,R]=Fun_NegLogLikelihood_Same(Dh,Resph,Theta,n,nugget,Mu,Sigma);
        PDF=mvnpdf(Resph,FT'*Mu,Sigma*R);

        rT=[rT, ComputeRmatrix(AFPoints,Dh(n,:),Theta)];
        [ZhPreds1,ZhCovs1]=Fun_GPPredictionGrid(AFPoints,Dh,Resph,Theta,Mu,Sigma,invR,invRRes,nugget,rT);
    end
    NewXhatPoints=[AFPoints;Dh;HistoryXhats];
    
    [ZhPreds2,ZhCovs2]=Fun_GPPrediction(Dh,Dh,Resph,Theta,Mu,Sigma,invR,invRRes,nugget);
    [ZhPreds3,ZhCovs3]=Fun_GPPrediction(HistoryXhats,Dh,Resph,Theta,Mu,Sigma,invR,invRRes,nugget);
    ZhPreds=[ZhPreds1;ZhPreds2;ZhPreds3]; ZhCovs=[ZhCovs1;ZhCovs2;ZhCovs3];
    Fval_ZhQuantile=ZhPreds+norminv(0.9)*ZhCovs.^0.5;
    [~,Sortidx]=sort(Fval_ZhQuantile);
    
    ObjZhQuantilePreds=@(x) Fun_GPQuantilePrediction(x,Dh,Resph,Theta,Mu,Sigma,invR,invRRes,nugget);
    XBestTry=zeros(NoS,Dim);
    fBestTry=zeros(NoS,1);   
    parfor id=1:NoS
        [XBestTry(id,:),fBestTry(id,1)]= patternsearch(ObjZhQuantilePreds,NewXhatPoints(Sortidx(id),:),[],[],[],[],lb,ub,[],options);
    end
    [Zhminhats(n,:),minidx]=min(fBestTry);
    Xhats(n,:)=XBestTry(minidx,:) ;
    Xhat_new=Xhats(n,:);
    %Stores GP model parameters and the MLE of the calibration parameter vector
    condRs(n,:)=condR;
    Yh_Xhats(n,:)=Simulator(Xhats(n,:),2,Case);
    SSETrue_Xhats(n,:)=sum( (Yh_Xhats(n,:)-PhysData).^2);
    
    disp(['Xhat_new and SSETrue_Xhats at  ' num2str(n) ' -iter '  num2str(Xhats(n,:))  '    ' num2str(SSETrue_Xhats(n,:))  ])
    
    if Budget<RatioCost
        break
    end
    HistoryXhats=[HistoryXhats; Xhat_new];
    NewAFPoints=[AFPoints;HistoryXhats];

    [MinZhRefEI,ZhCovs4]=Fun_GPPrediction(Xhat_new,Dh,Resph,Theta,Mu,Sigma,invR,invRRes,nugget);
    ZhPredsAF=[ZhPreds1;ZhPreds3;MinZhRefEI]; ZhCovsAF=[ZhCovs1;ZhCovs3;ZhCovs4]; 

    %%%%%%%Maximizes the EI AF and adds a follow-up design point
    fvals=Fun_MinusEIGrid(NewAFPoints,Dh,Resph,Theta,Mu,Sigma,invR,invRRes,MinZhRefEI,nugget,ZhPredsAF,ZhCovsAF);
    [~,Sortidx]=sort(fvals);
    
    MinusEIObj=@(TeD) Fun_MinusEI(TeD,Dh,Resph,Theta,Mu,Sigma,invR,invRRes,MinZhRefEI,nugget);    
    XBestTry=zeros(NoS,Dim);
    fBestTry=zeros(NoS,1);    
    parfor id=1:NoS
        [XBestTry(id,:),fBestTry(id,1)]= patternsearch(MinusEIObj,NewAFPoints(Sortidx(id),:),[],[],[],[],lb,ub,[],options);
    end
    [fBest,minidx]=min(fBestTry);
    NextPoint=XBestTry(minidx,:) ;
    AF=-fBest;AFs(n+1,:)=AF;
    disp(['Current budget=' num2str(Budget) '. The ' num2str(n+1) '-th run will be at point ' num2str(NextPoint,' %1.3f ')   ])
    n=n+1;
    Dh(n,:)=NextPoint;
    Yh(n,:)=Simulator(NextPoint,2,Case);
    Level(n,:)=2;
    
    Resph(n,:)=(mean(  (Yh(n,:)-PhysData).^2))^0.5;
    Budget=Budget-RatioCost ;
end

Sigmas(n,:)=Sigma;    Thetas(n,:)=Theta;    Mus(n,:)=Mu;condRs(n,:)=condR;
PDFs(n,:)=PDF;

%Stores design points and their corresponding simulator output
Resplh=Resph;
D=Dh;
RecordData.Dl=[];
RecordData.Yl=[];
RecordData.Respl=[];
RecordData.Dh=Dh;
RecordData.Yh=Yh;
RecordData.Resph=Resph;
RecordData.XTrue=XTrue;
RecordData.PhysData=PhysData;
RecordData.RatioCost=RatioCost;
RecordData.Budget=Budget;
RecordData.Yh_Xhats=Yh_Xhats;

%Stores GP model parameters and the MLE of the calibration parameter vector at all iterations with a table
RecordTable=table(D,Level,Resplh,Sigmas,Thetas,Mus,Xhats,Zhminhats,SSETrue_Xhats,condRs,PDFs,AFs);
end
%%
function [Mu,Sigma,Theta,invR,invRRes,condR,Objfval,PDF]=GPFit(Dh,Resph,nugget,Val)
[n,Dim]=size(Dh);

lb=[ (0.25)*ones(1,Dim) ];
ub=[ (15)*ones(1,Dim) ];
ObjS=@(Par) Fun_NegLogLikelihood(Dh,Resph,Par,n,nugget);
nvar=numel(ub);
options=optimoptions('patternsearch','disp','off');
if(Val==1)
HNoS=90;
else
HNoS=100;
end

Sobolset=sobolset(nvar,'Skip',1e3,'Leap',1e2);
if(Val==1)
StandardPoints=[net(Sobolset,7000*nvar)];
else
StandardPoints=[net(Sobolset,8000*nvar)];    
end

SIZE=size(StandardPoints,1);
fvals=zeros(SIZE,1);
parfor id=1:SIZE
    Point=lb + (ub-lb).*StandardPoints(id,:);
    fvals(id,1)=ObjS(Point);
end
[~,Sortfvalsidx]=sort(fvals);

OPt_StandardPoints=[ StandardPoints(Sortfvalsidx(1:HNoS),:)];
Remain_StandardPoints=StandardPoints(Sortfvalsidx((HNoS+1):end),:);
Count=0;
for kd=1:size(Remain_StandardPoints,1)
    if min(pdist2(Remain_StandardPoints(kd,:),OPt_StandardPoints),[],2)>sqrt(nvar*0.1^2)
        OPt_StandardPoints=[OPt_StandardPoints; Remain_StandardPoints(kd,:)];
        Count=Count+1;
        if Count==HNoS
            break
        end
    end
end

LengthOPt_StandardPoints=size(OPt_StandardPoints,1);
fBestTry=zeros(LengthOPt_StandardPoints,1);
XBestTry=zeros(LengthOPt_StandardPoints,nvar);
parfor id=1:LengthOPt_StandardPoints
    Point=lb+(ub-lb).*OPt_StandardPoints(id,:);
    [XBestTry(id,:),fBestTry(id,1)]= patternsearch(ObjS,Point,[],[],[],[],lb,ub,[],options);
end

[~,minidx]=min(fBestTry);
Parin=XBestTry(minidx,:) ;
[Objfval,Mu,Sigma,Theta,invR,invRRes,condR,FT,R]=ObjS(Parin);
PDF=mvnpdf(Resph,FT'*Mu,Sigma*R);
end
%%
function [fval,Mu,Sigma,Theta,invR,invRRes,condR,FT,R]=Fun_NegLogLikelihood(D,Resph,Theta,n,nugget)

FT=ones(1,n) ;
R=ComputeRmatrix2(D,Theta,nugget);
[invR,logdetR,condR]=invandlogdet(R);
FTinvR=FT*invR;
Mu=(FTinvR*Resph)/sum(invR,'all');

Res=Resph-Mu ;
invRRes=invR*Res;
Sigma=Res'*invRRes/ (n ) ;
fval=n*log(Sigma)+logdetR;

if any(isnan(invR),'all') || any(isinf(invR),'all') || Sigma==0
    fval=Inf;Mu=[];Sigma=[];Theta=[];invR=[];invRRes=[];condR=[];
    return
end
end
%%
function [invR,invRRes,condR,FT,R]=Fun_NegLogLikelihood_Same(Dh,Resph,Theta,n,nugget,Mu,Sigma)

FT=ones(1,n) ;
R=ComputeRmatrix2(Dh,Theta,nugget);
[invR,~,condR]=invandlogdet(R);

Res=Resph-Mu ;
invRRes=invR*Res;

if any(isnan(invR),'all') || any(isinf(invR),'all') || Sigma==0
    invR=[];invRRes=[];condR=[];
    return
end
end
%%
function [ZhPreds,ZhCovs,rT]=Fun_GPPrediction(TeD,Dh,Resph,Theta,Mu,Sigma,invR,invRRes,nugget)
rT=ComputeRmatrix(TeD,Dh,Theta);
rT_invR=rT*invR;
ZhPreds=Mu + rT*invRRes;
ZhCovs =Sigma * ( 1+nugget - sum(rT_invR.*rT,2) ) ;
ZhCovs=max(ZhCovs,0);
end

function [ZhPreds,ZhCovs]=Fun_GPPredictionGrid(TeD,Dh,Resph,Theta,Mu,Sigma,invR,invRRes,nugget,rT)
% rT=ComputeRmatrix(TeD,Dh,Theta);
rT_invR=rT*invR;
ZhPreds=Mu + rT*invRRes;
ZhCovs =Sigma * ( 1+nugget - sum(rT_invR.*rT,2) ) ;
ZhCovs=max(ZhCovs,0);
end
%%
function [RespQuantile]=Fun_GPQuantilePrediction(TeD,Dh,Resph,Theta,Mu,Sigma,invR,invRRes,nugget)
[RespPreds,RespCovs]=Fun_GPPrediction(TeD,Dh,Resph,Theta,Mu,Sigma,invR,invRRes,nugget);
RespQuantile=RespPreds+norminv(0.9)*RespCovs.^0.5;

end
%%
function Fval_MinusEI=Fun_MinusEI(TeD,Dh,Resph,Theta,Mu,Sigma,invR,invRRes,minResph,nugget)
[ZhPreds,ZhCovs]=Fun_GPPrediction(TeD,Dh,Resph,Theta,Mu,Sigma,invR,invRRes,nugget);
SD=ZhCovs.^0.5;
Bias=minResph-ZhPreds;
BOSD=Bias./SD;
EIfval=Bias.*normcdf(BOSD)+SD.*normpdf(BOSD);
EIfval(   logical( [ZhCovs==0] |  min(pdist2(TeD,Dh),[],2) < 10^(-3)  )    )=0;
Fval_MinusEI=-EIfval;
end

function Fval_MinusEI=Fun_MinusEIGrid(TeD,Dh,Resph,Theta,Mu,Sigma,invR,invRRes,minResph,nugget,ZhPreds,ZhCovs)
% [ZhPreds,ZhCovs]=Fun_GPPrediction(TeD,Dh,Resph,Theta,Mu,Sigma,invR,invRRes,nugget);
SD=ZhCovs.^0.5;
Bias=minResph-ZhPreds;
BOSD=Bias./SD;
EIfval=Bias.*normcdf(BOSD)+SD.*normpdf(BOSD);
EIfval(   logical( [ZhCovs==0] |  min(pdist2(TeD,Dh),[],2) < 10^(-3)  )    )=0;
Fval_MinusEI=-EIfval;
end