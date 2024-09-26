function [RecordTable,RecordData]=CalibrationSVD(DataInput,Val,percentage)
%percentage=0.99;
nugget=1e-6;
Dh=DataInput.Dh;
Yh=(DataInput.Yh)';%column vec

XTrue=DataInput.XTrue;
PhysData=(DataInput.PhysData)';%column vec
RatioCost=DataInput.RatioCost;
Budget=DataInput.Budget;
Case=DataInput.Case;
[n,Dim]=size(Dh);
Level=2*ones(n,1) ;

Yh_Mean=mean(Yh,2);
PhysData2=PhysData-Yh_Mean;
Yh_normalize=(Yh-Yh_Mean);
[U1,S1,~] = svd(Yh_normalize);
S1=diag(S1).^2;
temp1=cumsum(S1)/sum(S1);
p_eta=find(temp1>=percentage,1)
K1=U1(:,1:p_eta);
K2=U1(:,(p_eta+1):end);

Wh=((K1')*Yh_normalize)';%Yh_normalize'*K1;

WhError=Yh_normalize-K1*(Wh');
VarWhError=var(WhError(:),1);

c_col=(K1')*(Yh_Mean-PhysData);
d_col=(K2')*(Yh_Mean-PhysData);
dTd=(d_col')*d_col;

SSEs=sum((Yh-PhysData).^2,1);
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
%Bayesian optimization
HistoryXhats=[];
if(Val==1)
NoS=90;
else
NoS=100;    
end
while (1)
    %Fits the GP model and finds the minimum posterior mean
    if ZFit==1
        [Thetah,Muh,Sigmah,invRh,~,CondRh,invRhRes]=GPFitWh(Dh,Wh,nugget,Val);
        Sigmahs(n,:)=(Sigmah(:))';
        Thetahs(n,:)=(Thetah(:))';
        p_etas(n,:)=p_eta;
        ZFit=0;
        
        [fval1,WhPreds1,WhCovs1,rhT] = ESh_Fun(AFPoints,Dh,Wh,Thetah,Muh,Sigmah,invRh,nugget,invRhRes,c_col,dTd,VarWhError,p_eta);

    else
        [invRh,CondRh,invRhRes]=Fun_NegLogLikelihoodWh_Same(Thetah,Dh,Wh,nugget);
        
        for td=1:p_eta
            rhT{td}=[rhT{td}, ComputeRmatrix(AFPoints,Dh(n,:),Thetah(td,:))];
        end
        [fval1,WhPreds1,WhCovs1] = ESh_FunGrid(AFPoints,Dh,Wh,Thetah,Muh,Sigmah,invRh,nugget,invRhRes,c_col,dTd,VarWhError,p_eta,rhT);
    end
    
    NewXhatPoints=[AFPoints;Dh;HistoryXhats];
        
    ESh_Obj = @(x) ESh_Fun(x,Dh,Wh,Thetah,Muh,Sigmah,invRh,nugget,invRhRes,c_col,dTd,VarWhError,p_eta);
    [fval2,~,~] = ESh_Obj(Dh);
    [fval3,~,~] = ESh_Obj(HistoryXhats);
    
    [~,Sortidx]=sort([fval1(:);fval2(:);fval3(:)]);
    
    XBestTry=zeros(NoS,Dim);
    fBestTry=zeros(NoS,1);
    parfor id=1:NoS
        [XBestTry(id,:),fBestTry(id,1)]= patternsearch(ESh_Obj,NewXhatPoints(Sortidx(id),:),[],[],[],[],lb,ub,[],options);
    end
    [~,minidx]=min(fBestTry);
    Xhat_new=XBestTry(minidx,:) ;
    
    %Stores GP model parameters and the MLE of the calibration parameter vector
    Xhats(n,:)=Xhat_new;
    CondRhs(n,:)=CondRh;
    Yh_Xhats(n,:)=Simulator(Xhat_new,2,Case);
    SSETrue_Xhats(n,:)=sum( (  (Yh_Xhats(n,:)')-PhysData).^2);
    HistoryXhats=[HistoryXhats;Xhat_new];
    disp(['Xhat_new and SSETrue_Xhats at  ' num2str(n) ' -iter '  num2str(Xhats(n,:))  '    ' num2str(SSETrue_Xhats(n,:))  ])
    
    if Budget<RatioCost
        break
    end
    
    NewAFPoints=[AFPoints;HistoryXhats];
    f=(Yh-PhysData).^2;
    [fmin,~]=min(f,[],2);
    srfmin=fmin.^0.5;
    
    Obj_MinusEI_Wh= @(x) Fun_MinusEI_Wh(x,Dh,Wh,Thetah,Muh,Sigmah,invRh,nugget,invRhRes,PhysData2,K1,fmin,srfmin,VarWhError);
        
    MinusEI_fval1 = Fun_MinusEI_WhGrid(AFPoints,Dh,Wh,Thetah,Muh,Sigmah,invRh,nugget,invRhRes,PhysData2,K1,fmin,srfmin,VarWhError,WhPreds1,WhCovs1);
    MinusEI_fval4 = Obj_MinusEI_Wh(HistoryXhats);
    [~,Sortidx]=sort([MinusEI_fval1(:); MinusEI_fval4(:) ]);

    XBestTry=zeros(NoS,Dim);
    fBestTry=zeros(NoS,1);
    parfor id=1:NoS
        StartPoint=NewAFPoints(Sortidx(id),:);
        [XBestTry(id,:),fBestTry(id,1)]=patternsearch(Obj_MinusEI_Wh,StartPoint,[],[],[],[],lb,ub,[],options);
    end
    [~,minidx]=min(fBestTry);
    NextPoint=XBestTry(minidx,:);
    
    disp(['Current budget=' num2str(Budget) '. The ' num2str(n+1) '-th run will be at point ' num2str(NextPoint,' %1.3f ')     ])
    n=n+1;
    Dh(n,:)=NextPoint;
    Yh(:,n)=Simulator(NextPoint,2,Case);
    Yh_normalize(:,n)=(Yh(:,n)-Yh_Mean);%column vec
    
    Wh(n,:)=(Yh_normalize(:,n)')*K1;
    
    SSEs(n)=sum((Yh(:,n)-PhysData).^2);
    Level(n,:)=2;
    
    Budget=Budget-RatioCost ;
end

Sigmahs(n,:)=(Sigmah(:))';
Thetahs(n,:)=(Thetah(:))';
p_etas(n,:)=p_eta;
%Stores design points and their corresponding simulator output
SSEs=SSEs';
D=Dh;
RecordData.Dl=[];
RecordData.Yl=[];
RecordData.Respl=[];
RecordData.Dh=Dh;
RecordData.Yh=Yh;
RecordData.XTrue=XTrue;
RecordData.PhysData=PhysData;
RecordData.RatioCost=RatioCost;
RecordData.Budget=Budget;
RecordData.Yh_Xhats=Yh_Xhats;

%Stores GP model parameters and the MLE of the calibration parameter vector at all iterations with a table
RecordTable=table(D,Level,Sigmahs,Thetahs,SSEs,Xhats,SSETrue_Xhats,p_etas,CondRhs);
end
%%
function [fval,WhPreds,WhCovs,rhT] = ESh_Fun(TeD,Dh,Wh,Thetah,Muh,Sigmah,invRh,nugget,invRhRes,c_col,dTd,VarWhError,N)
[WhPreds,WhCovs,rhT]=GPPredictionWh(TeD,Dh,Wh,Thetah,Muh,Sigmah,invRh,nugget,invRhRes);
Part1=sum(( (WhPreds')+c_col).^2,1);%1-n
Part2=(sum(WhCovs,2))'+N*VarWhError;%1-n
fval=Part1+Part2+dTd;
end

function [fval,WhPreds,WhCovs] = ESh_FunGrid(TeD,Dh,Wh,Thetah,Muh,Sigmah,invRh,nugget,invRhRes,c_col,dTd,VarWhError,N,rhT)
[WhPreds,WhCovs]=GPPredictionWhGrid(TeD,Dh,Wh,Thetah,Muh,Sigmah,invRh,nugget,invRhRes,rhT);
Part1=sum(( (WhPreds')+c_col).^2,1);%1-n
Part2=(sum(WhCovs,2))'+N*VarWhError;%1-n
fval=Part1+Part2+dTd;
end
%%
function [Thetah,Muh,Sigmah,invRh,Objectiveh0,condRh0max,invRhRes]=GPFitWh(Dh,Wh,nugget,Val)

[~,Dim]=size(Dh);
[~,p_eta]=size(Wh);
lb=[ (0.25)*ones(1,Dim) ];
ub=[ (15)*ones(1,Dim) ];
nvar=numel(lb) ;
options=optimoptions('patternsearch','disp','off');
if(Val==1)
HNoS=90;
else
HNoS=100;
end
Sobolset=sobolset(nvar,'Skip',1e3,'Leap',1e2);
if(Val==1)
StandardPoints=[net(Sobolset,7000*nvar)];%[0 1]
else
StandardPoints=[net(Sobolset,8000*nvar)];%[0 1]    
end

SIZE=size(StandardPoints,1);
fvals=zeros(SIZE,1);
Muh=zeros(1,p_eta);
Thetah=zeros(p_eta,Dim);
Sigmah=zeros(1,p_eta);
invRh=cell(1,p_eta);
invRhRes=cell(1,p_eta);
condRh0max=-Inf;
for td=1:p_eta
    Wh0=Wh(:,td);
    ObjS=@(Thetah) Fun_NegLogLikelihoodWh0(Thetah,Dh,Wh0,nugget);
    
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
    Thetah0=XBestTry(minidx,:) ;
    [Objectiveh0,~,Sigmah0,invRh0,condRh0,invRhRes0]=ObjS(Thetah0);
    
    %Muh(td)=Muh0;=0;
    Thetah(td,:)=Thetah0;
    Sigmah(td)=Sigmah0;
    invRh{td}=invRh0;
    invRhRes{td}=invRhRes0;
    if(condRh0>condRh0max)
        condRh0max=condRh0;
    end
end

end

function [Objectiveh0,Muh0,Sigmah0,invRh0,condRh0,invRhRes0]=Fun_NegLogLikelihoodWh0(Thetah0,Dh,Wh0,nugget)
[nh,~]=size(Wh0);%p_eta=1
Rh0=ComputeRmatrix2(Dh,Thetah0,nugget);
[invRh0, logdetRh0,condRh0]=invandlogdet(Rh0);
Muh0=0;
Res0=Wh0;
invRhRes0=invRh0*Res0;
Sigmah0=(Res0')*invRhRes0/(nh );
Objectiveh0=(nh)*log(Sigmah0)+logdetRh0;
end

function [invRh,condRh0max,invRhRes]=Fun_NegLogLikelihoodWh_Same(Thetah,Dh,Wh,nugget)
[nh,p_eta]=size(Wh);
invRhRes=cell(1,p_eta);
invRh=cell(1,p_eta);
condRh0max=-Inf;
for td=1:p_eta
    Wh0=Wh(:,td);
    Thetah0=Thetah(td,:);
    Rh0=ComputeRmatrix2(Dh,Thetah0,nugget);
    [invRh0, ~,condRh0]=invandlogdet(Rh0);
    Res0=Wh0;
    invRhRes0=invRh0*Res0;
    invRhRes{td}=invRhRes0;
    invRh{td}=invRh0;
    if(condRh0>condRh0max)
        condRh0max=condRh0;
    end
end

end
%%
function [WhPreds,WhCovs,rhT]=GPPredictionWh(TeD,Dh,Wh,Thetah,Muh,Sigmah,invRh,nugget,invRhRes)
[~,p_eta]=size(Wh);
SIZE=size(TeD,1);
WhPreds=zeros(SIZE,p_eta);
WhCovs=zeros(SIZE,p_eta);
for td=1:p_eta
    Thetah0=Thetah(td,:);
    invRhRes0=invRhRes{td};
    invRh0=invRh{td};
    rhT{td}=ComputeRmatrix(TeD,Dh,Thetah0);
    WhPreds(:,td)=0 +rhT{td}*invRhRes0;
    Cov0=Sigmah(td)*(nugget+1-sum((rhT{td}*invRh0).*rhT{td},2));
    WhCovs(:,td)=max(Cov0,0);
end
end

function [WhPreds,WhCovs]=GPPredictionWhGrid(TeD,Dh,Wh,Thetah,Muh,Sigmah,invRh,nugget,invRhRes,rhT)
[~,p_eta]=size(Wh);
SIZE=size(TeD,1);
WhPreds=zeros(SIZE,p_eta);
WhCovs=zeros(SIZE,p_eta);

for td=1:p_eta
    invRhRes0=invRhRes{td};
    invRh0=invRh{td};
    % rhT{td}=ComputeRmatrix(TeD,Dh,Thetah0);
    WhPreds(:,td)=0 +rhT{td}*invRhRes0;
    Cov0=Sigmah(td)*(nugget+1-sum((rhT{td}*invRh0).*rhT{td},2));
    WhCovs(:,td)=max(Cov0,0);
end
end

%%
function [MinusEI_fval]= Fun_MinusEI_Wh(TeD,Dh,Wh,Thetah,Muh,Sigmah,invRh,nugget,invRhRes,PhysData2,K1,fmin,srfmin,VarWhError)
[WhPreds,WhCovs]=GPPredictionWh(TeD,Dh,Wh,Thetah,Muh,Sigmah,invRh,nugget,invRhRes);
EI_fval=zeros(size(TeD,1),1);
K12=K1.^2;

parfor id=1:size(TeD,1)
    YhPreds_id=K1*(WhPreds(id,:)');
    YhCovs_id=K12*(WhCovs(id,:)');
    
    YhCovs2=YhCovs_id+VarWhError;
    srYhCovs=sqrt(YhCovs2);
    Res=PhysData2-YhPreds_id;
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

function [MinusEI_fval] = Fun_MinusEI_WhGrid(TeD,Dh,Wh,Thetah,Muh,Sigmah,invRh,nugget,invRhRes,PhysData2,K1,fmin,srfmin,VarWhError,WhPreds,WhCovs)
% [WhPreds,WhCovs]=GPPredictionWh(TeD,Dh,Wh,Thetah,Muh,Sigmah,invRh,nugget,invRhRes);
EI_fval=zeros(size(TeD,1),1);
K12=K1.^2;

parfor id=1:size(TeD,1)
    YhPreds_id=K1*(WhPreds(id,:)');
    YhCovs_id=K12*(WhCovs(id,:)');
    
    YhCovs2=YhCovs_id+VarWhError;
    srYhCovs=sqrt(YhCovs2);
    Res=PhysData2-YhPreds_id;
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