function [RecordTable,RecordData]=CalibrationSRGP(DataInput)
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
    nlevel=2501;
elseif Dim==3
    nlevel=201;
end
AFPoints=(fullfact(nlevel*ones(1,Dim))-1)/(nlevel-1);
XhatPoints=AFPoints;


lb=0*ones(1,Dim);ub=1*ones(1,Dim);
options=optimoptions('patternsearch','disp','off');

HistoryXhats=[];
AFs(n,:)=0;
%Bayesian optimization
ZFit=1;
while (1)
    if ZFit==1
        [Mu,Sigma,Theta,invR,invRRes,condR,Objfval,PDF]=GPFit(Dh,Resph,nugget);
        ZFit=0;
        Sigmas(n,:)=Sigma;    Thetas(n,:)=Theta;    Mus(n,:)=Mu;
        PDFs(n,:)=PDF;
        
    else
        [invR,invRRes,condR,FT,R]=Fun_NegLogLikelihood_Same(Dh,Resph,Theta,n,nugget,Mu,Sigma);
        PDF=mvnpdf(Resph,FT'*Mu,Sigma*R);
    end
    NewXhatPoints=[XhatPoints;Dh;HistoryXhats];
    ObjResphQuantilePreds=@(x) Fun_GPQuantilePrediction(x,Dh,Resph,Theta,Mu,Sigma,invR,invRRes,nugget);
    
    [~,Sortidx]=sort(ObjResphQuantilePreds(NewXhatPoints));
    parfor id=1:90
        [XBestTry(id,:),fBestTry(id,1)]= patternsearch(ObjResphQuantilePreds,NewXhatPoints(Sortidx(id),:),[],[],[],[],lb,ub,[],options);
    end
    [Resphminhats(n,:),minidx]=min(fBestTry);
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
    MinResphRefEI=Fun_GPPrediction(Xhats(n,:),Dh,Resph,Theta,Mu,Sigma,invR,invRRes,nugget);
    %%%%%%%Maximizes the EI AF and adds a follow-up design point
    MinusEIObj=@(TeD) Fun_MinusEI(TeD,Dh,Resph,Theta,Mu,Sigma,invR,invRRes,MinResphRefEI,nugget);
    fvals=MinusEIObj(NewAFPoints);
    [~,Sortidx]=sort(fvals);
    
    parfor id=1:90
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
RecordData.ZNBC=[];
RecordData.Yh_Xhats=Yh_Xhats;

%Stores GP model parameters and the MLE of the calibration parameter vector at all iterations with a table
RecordTable=table(D,Level,Resplh,Sigmas,Thetas,Mus,Xhats,Resphminhats,SSETrue_Xhats,condRs,PDFs,AFs);
end
%%
function [Mu,Sigma,Theta,invR,invRRes,condR,Objfval,PDF]=GPFit(Dh,Resph,nugget)
[n,Dim]=size(Dh);

lb=[ (0.25)*ones(1,Dim) ];
ub=[ (15)*ones(1,Dim) ];
ObjS=@(Par) Fun_NegLogLikelihood(Dh,Resph,Par,n,nugget);
nvar=numel(ub);
options=optimoptions('patternsearch','disp','off');

Sobolset=sobolset(nvar,'Skip',1e3,'Leap',1e2);
StardardPoints=[net(Sobolset,7000*nvar)];

SIZE=size(StardardPoints,1);
fvals=zeros(SIZE,1);
parfor id=1:SIZE
    Point=lb + (ub-lb).*StardardPoints(id,:);
    fvals(id,1)=ObjS(Point);
end
[sortfvals,Sortfvalsidx]=sort(fvals);

OPt_StardardPoints=[ StardardPoints(Sortfvalsidx(1:90),:)];
Remain_StardardPoints=StardardPoints(Sortfvalsidx((90+1):end),:);
Count=0;
for kd=1:size(Remain_StardardPoints,1)
    if min(pdist2(Remain_StardardPoints(kd,:),OPt_StardardPoints),[],2)>sqrt(nvar*0.1^2)
        OPt_StardardPoints=[OPt_StardardPoints; Remain_StardardPoints(kd,:)];
        Count=Count+1;
        if Count==90
            break
        end
    end
end

parfor id=1:size(OPt_StardardPoints,1)
    Point=lb+(ub-lb).*OPt_StardardPoints(id,:);
    [XBestTry(id,:),fBestTry(id,1)]= patternsearch(ObjS,Point,[],[],[],[],lb,ub,[],options);
end

[fBest,minidx]=min(fBestTry);
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
[invR,logdetR,condR]=invandlogdet(R);

Res=Resph-Mu ;
invRRes=invR*Res;

if any(isnan(invR),'all') || any(isinf(invR),'all') || Sigma==0
    fval=Inf;Mu=[];Sigma=[];Theta=[];invR=[];invRRes=[];condR=[];
    return
end
end
%%
function [RespPreds,RespCovs]=Fun_GPPrediction(TeD,Dh,Resph,Theta,Mu,Sigma,invR,invRRes,nugget)
rT=ComputeRmatrix(TeD,Dh,Theta);
rT_invR=rT*invR;
RespPreds=Mu + rT*invRRes;
RespCovs =Sigma * ( 1+nugget - sum(rT_invR.*rT,2) ) ;
RespCovs=max(RespCovs,0);
end
%%
function [RespQuantile]=Fun_GPQuantilePrediction(TeD,Dh,Resph,Theta,Mu,Sigma,invR,invRRes,nugget)
[RespPreds,RespCovs]=Fun_GPPrediction(TeD,Dh,Resph,Theta,Mu,Sigma,invR,invRRes,nugget);
RespQuantile=RespPreds+norminv(0.9)*RespCovs.^0.5;

end
%%
function Fval_MinusEI=Fun_MinusEI(TeD,Dh,Resph,Theta,Mu,Sigma,invR,invRRes,minResph,nugget)

[RespPreds,RespCovs]=Fun_GPPrediction(TeD,Dh,Resph,Theta,Mu,Sigma,invR,invRRes,nugget);
SD=RespCovs.^0.5;
Bias=minResph-RespPreds;
EIfval=Bias.*normcdf(Bias./SD)+SD.*normpdf(Bias./SD);
EIfval(   logical( [RespCovs==0] |  min(pdist2(TeD,Dh),[],2) < 10^(-3)  )    )=0;
Fval_MinusEI=-EIfval;
end