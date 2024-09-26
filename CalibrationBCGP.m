function [RecordTable,RecordData]=CalibrationBCGP(DataInput,Val)
% Single fidelity calibration method : BC-GP
%           Resph: is a vector of HF SSE,
%           Xhats: is MLE of the calibration parameter vector at all iterations
%           Resphminhats: is posterior means/medians of the Resph at Xhats at all iterations
%           SSETrue_Xhats: is the true values of SSE evaluated at Xhats at all iterations
nugget=1e-6;
%Computes the untransformed HF response
Dh=DataInput.Dh;
Yh=DataInput.Yh;
XTrue=DataInput.XTrue;
PhysData=DataInput.PhysData;
RatioCost=DataInput.RatioCost;
Budget=DataInput.Budget;
Case=DataInput.Case;
[n,Dim]=size(Dh);
Level=2*ones(n,1) ;
ZNBC=1;

Resph=sum((Yh-PhysData).^2,2);
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
ZFit=1;

HistoryXhats=[];
AFs(n,:)=0;
if(Val==1)
NoS=90;
else
NoS=100;    
end
while (1)
    
    %Fits the GP model and finds the minimum posterior mean
    if ZFit==1
        [Mu,Sigma,Theta,invR,invRRes,condR,~,phi,PDF]=GPFit(Dh,Resph,ZNBC,nugget,Val);
        ZFit=0;
        Sigmas(n,:)=Sigma;    Thetas(n,:)=Theta;    Mus(n,:)=Mu;
        PDFs(n,:)=PDF;  phis(n,:)=phi;
        
        [ZhPreds1,ZhCovs1,rT]=Fun_PredictionZ(AFPoints,Dh,Resph,Theta,Mu,Sigma,invR,invRRes,nugget);

    else
        [invR,invRRes,condR,ZData,FT,R]=Fun_NegLogLikelihood_Same(Dh,Resph,Theta,phi,n,ZNBC,nugget,Mu,Sigma);
        [~,~,absdZ]=TransformData(Resph,phi,ZNBC);
        Prod_dZ=prod(absdZ);
        PDF=mvnpdf(ZData,FT'*Mu,Sigma*R)*Prod_dZ;

        rT=[rT, ComputeRmatrix(AFPoints,Dh(n,:),Theta)];
        [ZhPreds1,ZhCovs1]=Fun_PredictionZGrid(AFPoints,Dh,Resph,Theta,Mu,Sigma,invR,invRRes,nugget,rT);
    end

    NewXhatPoints=[AFPoints;Dh;HistoryXhats];
    
    [ZhPreds2,ZhCovs2]=Fun_PredictionZ(Dh,Dh,Resph,Theta,Mu,Sigma,invR,invRRes,nugget);
    [ZhPreds3,ZhCovs3]=Fun_PredictionZ(HistoryXhats,Dh,Resph,Theta,Mu,Sigma,invR,invRRes,nugget);
    ZhPreds=[ZhPreds1;ZhPreds2;ZhPreds3]; ZhCovs=[ZhCovs1;ZhCovs2;ZhCovs3];
    Fval_ZhQuantile=ZhPreds+norminv(0.9)*ZhCovs.^0.5;
    [~,Sortidx]=sort(Fval_ZhQuantile);

    ObjZhQuantile=@(x) Fun_ZhQuantile(x,Dh,Resph,Theta,Mu,Sigma,invR,invRRes,nugget);
    XBestTry=zeros(NoS,Dim);
    fBestTry=zeros(NoS,1);    
    parfor id=1:NoS
        [XBestTry(id,:),fBestTry(id,1)]= patternsearch(ObjZhQuantile,NewXhatPoints(Sortidx(id),:),[],[],[],[],lb,ub,[],options);
    end
    [minObjZhPreds,minidx]=min(fBestTry);
    Xhat_new=XBestTry(minidx,:) ;
    
    condRs(n,:)=condR;
    Xhats(n,:)=Xhat_new;
    Resphminhats(n,:)= TransformData_inv(minObjZhPreds,phi,ZNBC);
    
    %Evaluate the true SSE at the MLE of the calibration parameter vector
    Yh_Xhats(n,:)=Simulator(Xhat_new,2,Case);
    SSETrue_Xhats(n,:)=sum( (Yh_Xhats(n,:)-PhysData).^2);
    
    disp(['Xhat_new and SSETrue_Xhats at  ' num2str(n) ' -iter '  num2str(Xhats(n,:))  '    ' num2str(SSETrue_Xhats(n,:))  ])
    
    if Budget<RatioCost
        break
    end
    HistoryXhats=[HistoryXhats; Xhat_new];
    NewAFPoints=[AFPoints;HistoryXhats];

    [minZDataRef,ZhCovs4]=Fun_PredictionZ(Xhat_new,Dh,Resph,Theta,Mu,Sigma,invR,invRRes,nugget);
    ZhPredsAF=[ZhPreds1;ZhPreds3;minZDataRef]; ZhCovsAF=[ZhCovs1;ZhCovs3;ZhCovs4];

    fvals=Fun_MinusEIGrid(NewAFPoints,Dh,Resph,Theta,Mu,Sigma,invR,invRRes,minZDataRef,nugget,ZhPredsAF,ZhCovsAF);
    [~,Sortidx]=sort(fvals);
    
    %Adds a follow-up design point by Maximizing the AF
    MinusEIObj=@(TeD) Fun_MinusEI(TeD,Dh,Resph,Theta,Mu,Sigma,invR,invRRes,minZDataRef,nugget);
    XBestTry=zeros(NoS,Dim);
    fBestTry=zeros(NoS,1);    
    parfor id=1:NoS
        [XBestTry(id,:),fBestTry(id,1)]= patternsearch(MinusEIObj,NewAFPoints(Sortidx(id),:),[],[],[],[],lb,ub,[],options);
    end
    [fBest,minidx]=min(fBestTry);
    NextPoint=XBestTry(minidx,:);
    AF=-fBest;AFs(n+1,:)=AF;
    disp(['Current budget=' num2str(Budget) '. The ' num2str(n+1) '-th run will be at point ' num2str(NextPoint,' %1.3f ')     ])
    n=n+1;
    Dh(n,:)=NextPoint;
    Yh(n,:)=Simulator(NextPoint,2,Case);
    Level(n,:)=2;
    
    Resph(n,:)=sum((Yh(n,:)-PhysData).^2);
    Budget=Budget-RatioCost ;
end

%Stores design points and their corresponding simulator output
Sigmas(n,:)=Sigma;    Thetas(n,:)=Theta;    Mus(n,:)=Mu;    condRs(n,:)=condR;
PDFs(n,:)=PDF;  phis(n,:)=phi;

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
RecordTable=table(D,Level,Resplh,Sigmas,phis,Thetas,Mus,Xhats,Resphminhats,SSETrue_Xhats,condRs,PDFs,AFs);
end
%%
function [Mu,Sigma,Theta,invR,invRRes,condR,Objfval,phi,PDF]=GPFit(Dh,Resph,ZNBC,nugget,Val)
[n,Dim]=size(Dh);

lb=[ (0.25)*ones(1,Dim) -2];
ub=[ (15)*ones(1,Dim) 2];
ObjS=@(Par) Fun_NegLogLikelihood(Dh,Resph,Par(1:Dim),Par(Dim+1),n,ZNBC,nugget);
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
    StandardPoint=StandardPoints(id,:);
    Point=lb + (ub-lb).*StandardPoint;
    Point(end)=-0.5+1.5*StandardPoint(end);
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
    OPt_StandardPoint=OPt_StandardPoints(id,:);
    Point=lb+(ub-lb).*OPt_StandardPoint;
    Point(end)=-0.5+1.5*OPt_StandardPoint(end);
    [XBestTry(id,:),fBestTry(id,1)]= patternsearch(ObjS,Point,[],[],[],[],lb,ub,[],options);
end

[~,minidx]=min(fBestTry);
Parin=XBestTry(minidx,:) ;
[Objfval,Mu,Sigma,Theta,invR,invRRes,condR,phi,ZData,FT,R]=ObjS(Parin);
[~,~,absdZ]=TransformData(Resph,phi,ZNBC);
Prod_dZ=prod(absdZ);
PDF=mvnpdf(ZData,FT'*Mu,Sigma*R)*Prod_dZ;
end
%%
function [invR,invRRes,condR,ZData,FT,R]=Fun_NegLogLikelihood_Same(Dh,Resph,Theta,phi,n,ZNBC,nugget,Mu,Sigma)
[ZData,~]=TransformData(Resph,phi,ZNBC);
FT=ones(1,n) ;
R=ComputeRmatrix2(Dh,Theta,nugget);
[invR,~,condR]=invandlogdet(R);

Res=ZData-Mu ;
invRRes=invR*Res;
if any(isnan(invR),'all') || any(isinf(invR),'all') || Sigma==0
    invR=[];invRRes=[];condR=[];phi=[];
    return
end
end
%%
function [fval,Mu,Sigma,Theta,invR,invRRes,condR,phi,ZData,FT,R]=Fun_NegLogLikelihood(D,Resph,Theta,phi,n,ZNBC,nugget)

[ZData,SumLogdZ]=TransformData(Resph,phi,ZNBC);
FT=ones(1,n) ;
R=ComputeRmatrix2(D,Theta,nugget);
[invR,logdetR,condR]=invandlogdet(R);
FTinvR=FT*invR;
Mu=(FTinvR*ZData)/sum(invR,'all');

Res=ZData-Mu ;
invRRes=invR*Res;
Sigma=Res'*invRRes/ (n ) ;
fval=n*log(Sigma)+logdetR-2*SumLogdZ;

if any(isnan(invR),'all') || any(isinf(invR),'all') || Sigma==0
    fval=Inf;Mu=[];Sigma=[];Theta=[];invR=[];invRRes=[];condR=[];phi=[];
    return
end
end
%%
function [Fval_ZQuantile]=Fun_ZhQuantile(TeD,Dh,Resph,Theta,Mu,Sigma,invR,invRRes,nugget)
[ZPreds,ZCovs]=Fun_PredictionZ(TeD,Dh,Resph,Theta,Mu,Sigma,invR,invRRes,nugget);
Fval_ZQuantile=ZPreds+norminv(0.9)*ZCovs.^0.5;
end
%%
function [ZPreds,ZCovs,rT]=Fun_PredictionZ(TeD,Dh,Resph,Theta,Mu,Sigma,invR,invRRes,nugget)
rT=ComputeRmatrix(TeD,Dh,Theta);
rT_invR=rT*invR;
ZPreds=Mu + rT*invRRes;
ZCovs =Sigma * ( 1+nugget - sum(rT_invR.*rT,2) ) ;
ZCovs=max(ZCovs,0);
end
%%
function [ZPreds,ZCovs]=Fun_PredictionZGrid(TeD,Dh,Resph,Theta,Mu,Sigma,invR,invRRes,nugget,rT)
rT_invR=rT*invR;
ZPreds=Mu + rT*invRRes;
ZCovs =Sigma * ( 1+nugget - sum(rT_invR.*rT,2) ) ;
ZCovs=max(ZCovs,0);
end

%%
function Fval_MinusEI=Fun_MinusEI(TeD,Dh,Resph,Theta,Mu,Sigma,invR,invRRes,minZData,nugget)
[ZPreds,ZCovs]=Fun_PredictionZ(TeD,Dh,Resph,Theta,Mu,Sigma,invR,invRRes,nugget);
SD=ZCovs.^0.5;
Bias=minZData-ZPreds;
BOSD=Bias./SD;
EIfval=Bias.*normcdf(BOSD)+SD.*normpdf(BOSD);
EIfval(   logical( [ZCovs==0] |  min(pdist2(TeD,Dh),[],2) < 10^(-3)  )    )=0;
Fval_MinusEI=-EIfval;
end

function Fval_MinusEI=Fun_MinusEIGrid(TeD,Dh,Resph,Theta,Mu,Sigma,invR,invRRes,minZDataRef,nugget,ZPreds,ZCovs)
% [ZPreds,ZCovs]=Fun_PredictionZ(TeD,Dh,Resph,Theta,Mu,Sigma,invR,invRRes,nugget);
SD=ZCovs.^0.5;
Bias=minZDataRef-ZPreds;
BOSD=Bias./SD;
EIfval=Bias.*normcdf(BOSD)+SD.*normpdf(BOSD);
EIfval(   logical( [ZCovs==0] |  min(pdist2(TeD,Dh),[],2) < 10^(-3)  )    )=0;
Fval_MinusEI=-EIfval;
end