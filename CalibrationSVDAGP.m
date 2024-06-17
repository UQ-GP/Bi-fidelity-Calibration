function [RecordTable,RecordData]=CalibrationSVDAGP(DataInput,percentage)
% percentage=0.99;
nugget=1e-6;
Dl=DataInput.Dl;
Yl=(DataInput.Yl)';%column vec
Dh=DataInput.Dh;
Yh=(DataInput.Yh)';%column vec
nh=size(Dh,1);[nl,Dim]=size(Dl);
XTrue=DataInput.XTrue;
PhysData=(DataInput.PhysData)';%column vec
RatioCost=DataInput.RatioCost;
Budget=DataInput.Budget;
Case=DataInput.Case;

Yl_Mean=mean(Yl,2);
Yl_normalize=(Yl-Yl_Mean);
[Ul,Svl,~] = svd(Yl_normalize);
Svl=diag(Svl).^2;
temp1=cumsum(Svl)/sum(Svl);
p_etaL=find(temp1>=percentage,1)
Kl1=Ul(:,1:p_etaL);

Deltah=Yh-Yl(:,1:nh);
Deltah_Mean=mean(Deltah,2);
Deltah_normalize=(Deltah-Deltah_Mean);
[Uh,Svh,~] = svd(Deltah_normalize);
Svh=diag(Svh).^2;
temp1=cumsum(Svh)/sum(Svh);
p_etah=find(temp1>=percentage,1)
Kh1=Uh(:,1:p_etah);

Wl=((Kl1')*Yl_normalize)';
WDeltah=((Kh1')*Deltah_normalize)';

WlError=Yl_normalize-Kl1*(Wl');
VarWlError=var(WlError(:),1);

WDeltahError=Deltah_normalize-Kh1*(WDeltah');
VarWDeltahError=var(WDeltahError(:),1);


n=nl+nh;
D=[Dl; Dh];
Level=[ ones(nl,1) ; 2*ones(nh,1) ];

SSEhs=sum((Yh-PhysData).^2,1);
Budget=Budget-(nl*1+nh*RatioCost);
MinBudget=RatioCost+1;

DlnotDh=Dl(nh+1:end,:);
WlnotDh=Wl(nh+1:end,:);
YlnotDh=Yl(:,nh+1:end);

if Dim==2
    nlevel=2501;
elseif Dim==3
    nlevel=201;
end
AFPoints=(fullfact(nlevel*ones(1,Dim))-1)/(nlevel-1);
XhatPoints=AFPoints;

lb=0*ones(1,Dim);ub=1*ones(1,Dim);
options=optimoptions('patternsearch','disp','off');
ZFit=1;
%Bayesian optimization
HistoryXhats=[];
while (1)
    if ZFit==1
        [Thetal,Mul,Sigmal,invRl,~,condRl,invRlRes]=GPFitWh(Dl,Wl,nugget);
        [Thetah,Muh,Sigmah,invRh,~,condRh,invRhRes]=GPFitWh(Dh,WDeltah,nugget);
        Thetals(n,:)=(Thetal(:))';
        Sigmals(n,:)=(Sigmal(:))';
        Thetahs(n,:)=(Thetah(:))';
        Sigmahs(n,:)=(Sigmah(:))';
        
        p_etaLs(n,:)=p_etaL;
        p_etaHs(n,:)=p_etah;
        ZFit=0;
    else
        [invRl,condRl,invRlRes]=Fun_NegLogLikelihoodWh_Same(Thetal,Dl,Wl,nugget);
        [invRh,condRh,invRhRes]=Fun_NegLogLikelihoodWh_Same(Thetah,Dh,WDeltah,nugget);
    end
    
    NewXhatPoints=[XhatPoints;Dh;Dl;HistoryXhats];
    
    ESh_Obj = @(x) ESh_Fun(x,Dh,WDeltah,Thetah,Muh,Sigmah,invRh,invRhRes,Dl,Wl,Thetal,Mul,Sigmal,invRl,invRlRes,nugget,Kl1,Kh1,Yl_Mean,Deltah_Mean,PhysData,VarWlError,VarWDeltahError);
    
    [~,Sortidx]=sort(ESh_Obj(NewXhatPoints));
    XBestTry=zeros(90,Dim);
    fBestTry=zeros(90,1);
    parfor id=1:90
        [XBestTry(id,:),fBestTry(id,1)]= patternsearch(ESh_Obj,NewXhatPoints(Sortidx(id),:),[],[],[],[],lb,ub,[],options);
    end
    [~,minidx]=min(fBestTry);
    Xhat_new=XBestTry(minidx,:) ;
    condRlRhs(n-1:n,:)=[condRl,condRh;condRl,condRh];
    %Stores GP model parameters and the MLE of the calibration parameter vector
    Xhats(n-1:n,:)=[ Xhat_new; Xhat_new];
    Yh_Xhat_new=Simulator(Xhat_new,2,Case);
    SSEXhat_new=sum( ((Yh_Xhat_new')-PhysData).^2);
    Yh_Xhats(n-1:n,:)=[ Yh_Xhat_new ; Yh_Xhat_new];
    SSETrue_Xhats(n-1:n,:)=[SSEXhat_new; SSEXhat_new  ];
    
    HistoryXhats=[HistoryXhats;Xhat_new];
    disp(['Xhat_new and SSETrue_Xhats at  ' num2str(n) ' -iter '  num2str(Xhats(n,:))  '    ' num2str(SSETrue_Xhats(n,:))  ])
    
    if Budget<(MinBudget)
        break
    end
    %Adds a follow-up design point by Maximizing the AF
    NewAFPoints=[AFPoints;HistoryXhats];
    
    f=(Yh-PhysData).^2;
    [fmin,~]=min(f,[],2);
    srfmin=fmin.^0.5;
    
    Obj_MinusEI_Wh= @(x) Fun_MinusEI_Wh(x,Dh,WDeltah,Thetah,Muh,Sigmah,invRh,invRhRes,Dl,Wl,Thetal,Mul,Sigmal,invRl,invRlRes,nugget,Kl1,Kh1,Yl_Mean,Deltah_Mean,fmin,srfmin,PhysData,VarWlError,VarWDeltahError);
    fval_MinusEIObj= Obj_MinusEI_Wh(NewAFPoints);
    [~,Sortidx]=sort(fval_MinusEIObj);
    
    clear XBestTry fBestTry
    parfor id=1:90
        StartPoint=NewAFPoints(Sortidx(id),:);
        [XBestTry(id,:),fBestTry(id,1)]=patternsearch(Obj_MinusEI_Wh,StartPoint,[],[],[],[],lb,ub,[],options);
    end
    [~,minidx]=min(fBestTry);
    NextPoint=XBestTry(minidx,:);
    
    disp(['Current budget=' num2str(Budget) '. The ' num2str(n+1) '-th run will be at point ' num2str(NextPoint,' %1.3f ')     ])
    Dh(nh+1,:)=NextPoint;
    Yh(:,nh+1)=(Simulator(NextPoint,2,Case))';
    
    Dl=[Dh;DlnotDh];
    Yl=[Yl(:,1:nh), (Simulator(NextPoint,1,Case))',YlnotDh];
    Yl_normalize_New=(Yl(:,nh+1)-Yl_Mean);%column vec
    Wl_New=(Yl_normalize_New')*Kl1;
    Wl=[Wl(1:nh,:); Wl_New; WlnotDh];
    
    Deltah(:,nh+1)=Yh(:,nh+1)-Yl(:,nh+1);
    Deltah_normalize(:,nh+1)=Deltah(:,nh+1)-Deltah_Mean;
    WDeltah(nh+1,:)=(Deltah_normalize(:,nh+1)')*Kh1;
    
    D=[D ; NextPoint; NextPoint ];
    Level=[ Level ; 2 ;1];
    
    SSEhs=[SSEhs sum((Yh(:,end)-PhysData).^2)];
    
    Budget=Budget-RatioCost-1;
    nh=nh+1;
    nl=nl+1;
    n=n+2;
end
Thetals(n,:)=(Thetal(:))';
Sigmals(n,:)=(Sigmal(:))';
Thetahs(n,:)=(Thetah(:))';
Sigmahs(n,:)=(Sigmah(:))';

p_etaLs(n,:)=p_etaL;
p_etaHs(n,:)=p_etah;

%Stores design points and their corresponding simulator output
SSEhs=SSEhs';
RecordData.Dl=Dl;
RecordData.Yl=Yl;
RecordData.Respl=[];
RecordData.Dh=Dh;
RecordData.Yh=Yh;
RecordData.XTrue=XTrue;
RecordData.PhysData=PhysData;
RecordData.RatioCost=RatioCost;
RecordData.Budget=Budget;
RecordData.ZNBC=[];
RecordData.Yh_Xhats=Yh_Xhats;

RecordTable=table(D,Level,Thetals,Sigmals,Thetahs,Sigmahs,Xhats,SSETrue_Xhats,p_etaLs,p_etaHs,condRlRhs);
end
%%
function fval =         ESh_Fun(TeD,Dh,WDeltah,Thetah,Muh,Sigmah,invRh,invRhRes,Dl,Wl,Thetal,Mul,Sigmal,invRl,invRlRes,nugget,Kl1,Kh1,Yl_Mean,Deltah_Mean,PhysData,VarWlError,VarWDeltahError)
[WlPreds,WlCovs]=GPPredictionW0(TeD,Dl,Wl,Thetal,Mul,Sigmal,invRl,nugget,invRlRes);
[WDeltaPreds,WDeltaCovs]=GPPredictionW0(TeD,Dh,WDeltah,Thetah,Muh,Sigmah,invRh,nugget,invRhRes);
fval=zeros(size(TeD,1),1);
parfor id=1:size(TeD,1)
    YlPreds_id=Kl1*(WlPreds(id,:)')+Yl_Mean;
    DeltaPreds_id=Kh1*(WDeltaPreds(id,:)')+Deltah_Mean;
    YhPreds_id=YlPreds_id+DeltaPreds_id;
    YhCovs_id= (Kl1.^2)*(WlCovs(id,:)')+VarWlError+(Kh1.^2)*(WDeltaCovs(id,:)')+VarWDeltahError;
    
    A=YhPreds_id-PhysData;
    Part1=sum(A.^2,1);%1-n
    Part2=sum(YhCovs_id,1);%1-n
    fval(id)=Part1+Part2;
end
end
%%
function [YhPreds,YhCovs,WlPreds,WlCovs,WDeltaPreds,WDeltaCovs]=GPPredictionYh(TeD,Dh,WDeltah,Thetah,Muh,Sigmah,invRh,invRhRes,Dl,Wl,Thetal,Mul,Sigmal,invRl,invRlRes,nugget,Kl1,Kh1,Yl_Mean,Deltah_Mean,VarWlError,VarWDeltahError)
[WlPreds,WlCovs]=GPPredictionW0(TeD,Dl,Wl,Thetal,Mul,Sigmal,invRl,nugget,invRlRes);
[WDeltaPreds,WDeltaCovs]=GPPredictionW0(TeD,Dh,WDeltah,Thetah,Muh,Sigmah,invRh,nugget,invRhRes);
YlPreds=Kl1*(WlPreds')+Yl_Mean;
DeltaPreds=Kh1*(WDeltaPreds')+Deltah_Mean;
YhPreds=YlPreds+DeltaPreds;
YhCovs= (Kl1.^2)*(WlCovs')+VarWlError+(Kh1.^2)*(WDeltaCovs')+VarWDeltahError;
end
%%
function [WPreds,WCov]=GPPredictionW0(TeD,Dh,Wh,Thetah,Muh,Sigmah,invRh,nugget,invRRes)
[~,N]=size(Wh);
SIZE=size(TeD,1);
WPreds=zeros(SIZE,N);
WCov=zeros(SIZE,N);
for td=1:N
    Thetah0=Thetah(td,:);
    invRhRes0=invRRes{td};
    invRh0=invRh{td};
    rhT=ComputeRmatrix(TeD,Dh,Thetah0);
    WPreds(:,td)=  rhT*invRhRes0;
    WCov(:,td)=max(0,Sigmah(td)*(nugget+1-sum((rhT*invRh0).*rhT,2)));
end
end

%%
function MinusEI_fval= Fun_MinusEI_Wh(TeD,Dh,WDeltah,Thetah,Muh,Sigmah,invRh,invRhRes,Dl,Wl,Thetal,Mul,Sigmal,invRl,invRlRes,nugget,Kl1,Kh1,Yl_Mean,Deltah_Mean,fmin,srfmin,PhysData,VarWlError,VarWDeltahError)
[WlPreds,WlCovs]=GPPredictionW0(TeD,Dl,Wl,Thetal,Mul,Sigmal,invRl,nugget,invRlRes);%nTest-p_etaL
[WDeltaPreds,WDeltaCovs]=GPPredictionW0(TeD,Dh,WDeltah,Thetah,Muh,Sigmah,invRh,nugget,invRhRes);%%nTest-p_etaH

EI_fval=zeros(size(TeD,1),1);
parfor id=1:size(TeD,1)
    
    YlPreds_id=Kl1*(WlPreds(id,:)')+Yl_Mean;
    DeltaPreds_id=Kh1*(WDeltaPreds(id,:)')+Deltah_Mean;
    YhPreds_id=YlPreds_id+DeltaPreds_id;
    YhCovs_id= (Kl1.^2)*(WlCovs(id,:)')+VarWlError+(Kh1.^2)*(WDeltaCovs(id,:)')+VarWDeltahError;
    
    YhCovs2=YhCovs_id;
    srYhCovs=sqrt(YhCovs2);
    Res=PhysData-YhPreds_id;
    Qplus=(Res+srfmin)./srYhCovs;
    Qminus=(Res-srfmin)./srYhCovs;
    EIPart1=(fmin-Res.^2-YhCovs2).*(normcdf(Qplus)-normcdf(Qminus));
    EIPart2=(srfmin-Res).*normpdf(Qplus)+(srfmin+Res).*normpdf(Qminus);
    EIs=EIPart1+srYhCovs.*EIPart2;
    EI_fval(id,:)=sum(EIs);
end
EI_fval(min(pdist2(TeD,Dh),[],2) < 10^(-3))=0;
MinusEI_fval=-EI_fval;
end
%%
function [Thetah,Muh,Sigmah,invRh,Objectiveh0,condRh0,invRhRes]=GPFitWh(Dh,Wh,nugget)

[~,Dim]=size(Dh);
[~,N]=size(Wh);

lb=[ (0.25)*ones(1,Dim) ];
ub=[ (15)*ones(1,Dim) ];
nvar=numel(lb) ;
options=optimoptions('patternsearch','disp','off');
Sobolset=sobolset(nvar,'Skip',1e3,'Leap',1e2);
StandardPoints=[net(Sobolset,7000*nvar)];
SIZE=size(StandardPoints,1);
fvals=zeros(SIZE,1);
invRh=cell(1,N);
invRhRes=cell(1,N);
Muh=zeros(1,N);
Thetah=zeros(N,Dim);
Sigmah=zeros(1,N);

for td=1:N
    Wh0=Wh(:,td);
    ObjS=@(Thetah0) Fun_NegLogLikelihoodWh0(Thetah0,Dh,Wh0,nugget);
    
    parfor id=1:SIZE
        Point=lb + (ub-lb).*StandardPoints(id,:);
        fvals(id,1)=ObjS(Point);
    end
    [~,Sortfvalsidx]=sort(fvals);
    
    OPt_StandardPoints=[ StandardPoints(Sortfvalsidx(1:90),:)];
    Remain_StandardPoints=StandardPoints(Sortfvalsidx((90+1):end),:);
    Count=0;
    for kd=1:size(Remain_StandardPoints,1)
        if min(pdist2(Remain_StandardPoints(kd,:),OPt_StandardPoints),[],2)>sqrt(nvar*0.1^2)
            OPt_StandardPoints=[OPt_StandardPoints; Remain_StandardPoints(kd,:)];
            Count=Count+1;
            if Count==90
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
    Thetah0=XBestTry(minidx,:) ;
    [Objectiveh0,Muh0,Sigmah0,invRh0,condRh0,invRhRes0]=ObjS(Thetah0);%Fun_NegLogLikelihoodWh0(Thetah,Dh,Wh0,Fh,nugget);
    
    Muh(td)=Muh0;
    Thetah(td,:)=Thetah0;
    Sigmah(td)=Sigmah0;
    invRh{td}=invRh0;
    invRhRes{td}=invRhRes0;
end

end

function [Objectiveh0,Muh0,Sigmah0,invRh0,condRh0,invRhRes0]=Fun_NegLogLikelihoodWh0(Thetah0,Dh,Wh0,nugget)
[nh,~]=size(Wh0);
Rh0=ComputeRmatrix(Dh,Dh,Thetah0)+nugget*eye(nh);
[invRh0, logdetRh0,condRh0]=invandlogdet(Rh0);
Muh0=0;
Res0=Wh0;
invRhRes0=invRh0*Res0;
Sigmah0=(Res0')*invRhRes0/nh;
Objectiveh0=(nh)*log(Sigmah0) + logdetRh0;
end

function [invRh,condRh0,invRhRes]=Fun_NegLogLikelihoodWh_Same(Thetah,Dh,Wh,nugget)
[nh,N]=size(Wh);
invRhRes=cell(1,N);
invRh=cell(1,N);
for td=1:N
    Wh0=Wh(:,td);
    Thetah0=Thetah(td,:);
    Rh0=ComputeRmatrix(Dh,Dh,Thetah0)+nugget*eye(nh);
    [invRh0, ~,condRh0]=invandlogdet(Rh0);
    Res0=Wh0;
    invRhRes0=invRh0*Res0;
    invRhRes{td}=invRhRes0;
    invRh{td}=invRh0;
end
end