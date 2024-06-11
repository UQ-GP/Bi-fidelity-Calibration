function [RecordTable,RecordData]=CalibrationAGP(DataInput,ZNBC,ZMLFSSEorZLFSSE)
% Bi-fidelity calibration methods that include MBC-AGP, BC-AGP, MID-AGP and SR-AGP methods
%           Resph: is the HF SSE
%           Respl: is the LF SSE or modified LF SSE
%           Xhats: is MLE of the calibration parameter vector at all iterations
%           Resphminhats: is posterior means/medians of the Resph at Xhats at all iterations
%           SSETrue_Xhats: is the true values of SSE evaluated at Xhats at all iterations
nugget=1e-6;
%Computes the untransformed HF and the LF responses
Dl=DataInput.Dl;
Yl=DataInput.Yl;
Dh=DataInput.Dh;
Yh=DataInput.Yh;
XTrue=DataInput.XTrue;
PhysData=DataInput.PhysData;
RatioCost=DataInput.RatioCost;
Budget=DataInput.Budget;
Case=DataInput.Case;

[nl,Dim]=size(Dl);
[nh,~]=size(Dh);
n=nl+nh;
D=[Dl; Dh];
Level=[ ones(nl,1) ; 2*ones(nh,1) ];
if ZMLFSSEorZLFSSE==1
    NTimePoints=numel(PhysData);
    [~,idxinDl,idxinDh]=intersect(Dl,Dh,'rows','stable');
    SameYl=Yl(idxinDl,:);
    SameYh=Yh(idxinDh,:);
    Ones2=ones(numel(idxinDh),1);
    for kd=1:NTimePoints
        Ylkd=SameYl(:,kd);
        Yhkd=SameYh(:,kd);
        MatrixX=[Ones2,Ylkd];
        if Case==3
            if all(Ylkd<10^(-12)) 
                ai_bi(:,kd)=[0,1];
                Sum_ErrorYlYh0=sum(abs(Yhkd-Ylkd)) ;
                if Sum_ErrorYlYh0>0
                    return
                end
            else
                ai_bi(:,kd)=regress(Yhkd,MatrixX);
            end
        else
            ai_bi(:,kd)=regress(Yhkd,MatrixX);
        end
    end
    YlModified=ai_bi(1,:)+Yl.*ai_bi(2,:);
    Respl=sum((YlModified-PhysData).^2,2);
    Resph=sum((Yh-PhysData).^2,2);
elseif ZMLFSSEorZLFSSE==0
    Resph=sum((Yh-PhysData).^2,2);
    Respl=sum((Yl-PhysData).^2,2);
end

Resplh= [Respl;Resph ];
Budget=Budget-(nl*1+nh*RatioCost);

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
AFs(n,1)=0;
while(1)
    %Fits the GP model and finds the minimum posterior mean
    if ZFit==1
        [Objfval,Sigmal,Thetal,Rho,RatioSigmah,Sigmah,Thetah,phi,Mu,V,invV,invVRes,condV,PDF]= AGPFit(Dl,Respl,Dh,Resph,ZNBC,nugget);%@@@
        ZFit=0;
        Objfvals(n,:)=Objfval;    Sigmals(n,:)=Sigmal;    Thetals(n,:)=Thetal;    Rhos(n,:)=Rho;
        RatioSigmahs(n,:)=RatioSigmah;    Sigmahs(n,:)=Sigmah;    Thetahs(n,:)=Thetah;
        phis(n,:)=phi;      Mus(n,:)=Mu;
        PDFs(n,:)=PDF;
    else
        [V,invV,invVRes,condV,invVprime,invVprimeRes,ZData,F]=Fun_NegLogLikelihoodSame(Dl,Respl,Dh,Resph,ZNBC,Thetal,Thetah,Rho,RatioSigmah,phi,nl,nh,nugget,Sigmal,Sigmah,Mu);
        [~,~,absdZ]=TransformData([Respl;Resph],phi,ZNBC);
        Prod_dZ=prod(absdZ);
        PDF=mvnpdf(ZData,F*Mu,V)*Prod_dZ;
    end
    NewXhatPoints=[XhatPoints;Dl;Dh;HistoryXhats];
    Obj_Fun_ZhQuantile=@(x) Fun_ZhQuantile(x,2,Dl,Respl,Dh,Resph,Sigmal,Thetal,Rho,Sigmah,Thetah,phi,ZNBC,Mu,invV,invVRes,nugget);
    [~,Sortidx]=sort(Obj_Fun_ZhQuantile(NewXhatPoints));
    parfor id=1:90
        [XBestTry(id,:),fBestTry(id,1)]= patternsearch(Obj_Fun_ZhQuantile,NewXhatPoints(Sortidx(id),:),[],[],[],[],lb,ub,[],options);
    end
    [MinZhQuantile,minidx]=min(fBestTry);
    Xhat_new=XBestTry(minidx,:) ;
    
    %Stores GP model parameters and the MLE of the calibration parameter vector
    Xhats(n,:)=Xhat_new;
    Resphminhats(n,:)=TransformData_inv(MinZhQuantile,phi,ZNBC);
    condVs(n,:)=condV;
    %Evaluate the true SSE at the MLE of the calibration parameter vector
    Yh_Xhats(n,:)=Simulator(Xhat_new,2,Case);
    SSETrue_Xhats(n,:)=sum( (Yh_Xhats(n,:)-PhysData).^2);
    
    disp(['Xhat_new and SSETrue_Xhats at  ' num2str(n) ' -iter '  num2str(Xhats(n,:))  '    ' num2str(SSETrue_Xhats(n,:))  ])
    
    if Budget<1
        break
    end
    HistoryXhats=[HistoryXhats; Xhat_new];
    %Adds a follow-up design point by Maximizing the AF
    [NextPoint,NextLevel,AF]=SequentialRun_AEI(Budget,RatioCost,Dl,Respl,Dh,Resph,Sigmal,Thetal,Rho,Sigmah,Thetah,phi,ZNBC,Mu,invV,invVRes,[AFPoints;HistoryXhats],Xhat_new,nugget);
    AF;AFs(n+1,:)=AF;
    disp(['Current budget=' num2str(Budget) '. The ' num2str(n+1) '-th run will be at point ' num2str(NextPoint,' %1.3f ') ' at Level ' num2str(NextLevel)  ])
    n=n+1;
    if NextLevel==1%Adds a LF design point
        nl=nl+1;
        Dl(nl,:)=NextPoint;
        Yl(nl,:)=Simulator(NextPoint,1,Case);
        if ZMLFSSEorZLFSSE==1
            YlModified=ai_bi(1,:)+Yl.*ai_bi(2,:);
            Respl=sum((YlModified-PhysData).^2,2);
        elseif ZMLFSSEorZLFSSE==0
            Respl=sum((Yl-PhysData).^2,2);
        end
        Budget=Budget-1;
        Resplh(n,1)=Respl(nl,:);
        
    else%Adds a HF design point
        nh=nh+1;
        Dh(nh,:)=NextPoint;
        Yh(nh,:)=Simulator(NextPoint,2,Case);
        Resph=sum((Yh-PhysData).^2,2);
        Budget=Budget-RatioCost;
        Resplh(n,1)=Resph(nh,:);
    end
    
    D(n,:)=NextPoint;
    Level(n,:)=NextLevel;
end
Objfvals(n,:)=Objfval;    Sigmals(n,:)=Sigmal;    Thetals(n,:)=Thetal;    Rhos(n,:)=Rho;
RatioSigmahs(n,:)=RatioSigmah;    Sigmahs(n,:)=Sigmah;    Thetahs(n,:)=Thetah;
phis(n,:)=phi;      Mus(n,:)=Mu;    condVs(n,:)=condV;
PDFs(n,:)=PDF;

%Stores design points and their corresponding simulator output
RecordData.Dl=Dl;
RecordData.Yl=Yl;
RecordData.Respl=Respl;
RecordData.Dh=Dh;
RecordData.Yh=Yh;
RecordData.Resph=Resph;
RecordData.XTrue=XTrue;
RecordData.PhysData=PhysData;
RecordData.RatioCost=RatioCost;
RecordData.Budget=Budget;
RecordData.Yh_Xhats=Yh_Xhats;

%Stores GP model parameters and the MLE of the calibration parameter vector at all iterations with a table
RecordTable=table(D,Level,Resplh,Objfvals,Sigmals,Thetals,Rhos,RatioSigmahs,Sigmahs,Thetahs,phis,Mus,Xhats,Resphminhats,SSETrue_Xhats,condVs,PDFs,AFs);
end
%%
function [Objfval,Sigmal,Thetal,Rho,RatioSigmah,Sigmah,Thetah,phi,Mu,V,invV,invVRes,condV,PDF]=AGPFit(Dl,Respl,Dh,Resph,ZNBC,nugget)%@@@
[nl,Dim]=size(Dl);
[nh,~]=size(Dh);

if ZNBC==1
    lb=[    (0.25)*ones(1,Dim)       (0.25)*ones(1,Dim)     0         -6  -2 ];
    ub=[      (15)*ones(1,Dim)    (15)*ones(1,Dim)         2         2    2 ];
    Obj_NegLogL=@(Par) Fun_NegLogLikelihood(Dl,Respl,Dh,Resph,ZNBC,Par(1:Dim),Par((Dim+1):(2*Dim)),Par(2*Dim+1),10^Par(2*Dim+2),Par(2*Dim+3),nl,nh,nugget);
    nvar=numel(lb);
    Sobolset=sobolset(nvar,'Skip',1e3,'Leap',1e2);
    
    StandardPoints=[net(Sobolset,7000*nvar)];%[0 1]
    
    idxcolumn2=[2*Dim+1 2*Dim+2 2*Dim+3];
    idxcolumn1=setdiff(1:nvar,idxcolumn2);
    Points=StandardPoints;
    Points(:,idxcolumn1)=lb(idxcolumn1)+(ub(idxcolumn1)-lb(idxcolumn1)).*StandardPoints(:,idxcolumn1);   %[lb ub]
    Points(:,idxcolumn2) =[0.5 -6 -0.5] +([1.5 0 1]-[0.5 -6 -0.5]).*StandardPoints(:,idxcolumn2);   %[lb ub]      rho log10 gamma phi   [0.5 1.5] and [-6 0 ] [-.5 1]
    
elseif ZNBC==0 || ZNBC==2
    lb=[    (0.25)*ones(1,Dim)       (0.25)*ones(1,Dim)          0         -6       ];
    ub=[      (15)*ones(1,Dim)    (15)*ones(1,Dim)             2          2     ];
    Obj_NegLogL=@(Par) Fun_NegLogLikelihood(Dl,Respl,Dh,Resph,ZNBC,Par(1:Dim),Par((Dim+1):(2*Dim)),Par(2*Dim+1),10^Par(2*Dim+2),0,nl,nh,nugget);
    nvar=numel(lb);
    Sobolset=sobolset(nvar,'Skip',1e3,'Leap',1e2);
    
    
    StandardPoints=[net(Sobolset,7000*nvar)];
    
    idxcolumn2=[2*Dim+1 2*Dim+2];
    idxcolumn1=setdiff(1:nvar,idxcolumn2);
    Points=StandardPoints;
    Points(:,idxcolumn1)=lb(idxcolumn1)+(ub(idxcolumn1)-lb(idxcolumn1)).*StandardPoints(:,idxcolumn1);   %[lb ub]
    Points(:,idxcolumn2) =[0.5 -6] +([1.5 0]-[0.5 -6]).*StandardPoints(:,idxcolumn2);   %[lb ub]      rho log10 gamma   [0.5 1.5] and [-6 0 ]
    
end
options=optimoptions('patternsearch','disp','off','MaxFunctionEvaluations',2000*nvar,'MaxIterations',100*nvar); %Default optimoptions

SIZE=size(StandardPoints,1);
fvals=zeros(SIZE,1);
parfor id=1:SIZE
    fvals(id,1)=Obj_NegLogL(Points(id,:));
end
[sortfvals,Sortfvalsidx]=sort(fvals);
clear Point

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

Points=OPt_StandardPoints;
Points(:,idxcolumn1)=lb(idxcolumn1)+(ub(idxcolumn1)-lb(idxcolumn1)).*OPt_StandardPoints(:,idxcolumn1);   %[lb ub]
if ZNBC==1
    Points(:,idxcolumn2) =[0.5 -6 -0.5] +([1.5 0 1]-[0.5 -6 -0.5]).*OPt_StandardPoints(:,idxcolumn2);   %[lb ub]      rho log10 gamma phi   [0.5 1.5] and [-6 0 ] [-.5 1]
else
    Points(:,idxcolumn2) =[0.5 -6] +([1.5 0]-[0.5 -6]).*OPt_StandardPoints(:,idxcolumn2);   %[lb ub]      rho log10 gamma   [0.5 1.5] and [-6 0 ]
end

parfor id=1:size(OPt_StandardPoints,1)
    [XBestTry(id,:),fBestTry(id,1)]= patternsearch(Obj_NegLogL,Points(id,:),[],[],[],[],lb,ub,[],options); %#ok<PFBNS>
end
[~,minidx]=min(fBestTry);
ParIn=XBestTry(minidx,:) ;

[Objfval,Sigmal,Thetal,Rho,RatioSigmah,Sigmah,Thetah,phi,Mu,V,invV,invVRes,condV,invVprime,invVprimeRes,ZData,F]=Obj_NegLogL(ParIn);
[~,~,absdZ]=TransformData([Respl;Resph],phi,ZNBC);
Prod_dZ=prod(absdZ);
PDF=mvnpdf(ZData,F*Mu,V)*Prod_dZ;
end
%%
function [Objfval,Sigmal,Thetal,Rho,RatioSigmah,Sigmah,Thetah,phi,Mu,V,invV,invVRes,condV,invVprime,invVprimeRes,ZData,F]=Fun_NegLogLikelihood(Dl,Respl,Dh,Resph,ZNBC,Thetal,Thetah,Rho,RatioSigmah,phi,nl,nh,nugget)%@@@

[ZData,SumLogdZ]=TransformData([Respl;Resph],phi,ZNBC);
VprimeDl_Dh=Rho*ComputeRmatrix(Dl,Dh,Thetal);
Vprime=[ ComputeRmatrix2(Dl,Thetal,nugget) ,           VprimeDl_Dh;
    VprimeDl_Dh' , Rho^2*ComputeRmatrix(Dh,Dh,Thetal)+RatioSigmah*ComputeRmatrix2(Dh,Thetah,nugget) ];

F=[                ones(nl,1) , zeros(nl,1);
    Rho*ones(nh,1) , ones(nh,1)];
[invVprime,logdetVprime,condV]=invandlogdet(Vprime);
if any(isnan(invVprime),'all')  ||  any(isinf(invVprime),'all')
    Objfval=Inf;Sigmal=[];Thetal=[];Rho=[];RatioSigmah=[];Sigmah=[];Thetah=[];phi=[];Mu=[];V=[];invV=[];invVRes=[];condV=[];invVprime=[];invVprimeRes=[];ZData=[];
    return
end
FT_invVprime=F'*invVprime;
[inv1,~]=invandlogdet(FT_invVprime*F);

Mu=inv1*FT_invVprime*ZData;
Res=ZData-F*Mu;

invVprimeRes=invVprime*Res;
Sigmal=(Res'*invVprimeRes)/(nl+nh);

Objfval=logdetVprime+(nl+nh)*log(Sigmal)     - 2*SumLogdZ;

Sigmah=RatioSigmah*Sigmal;
V=Vprime*Sigmal;
invV=invVprime/Sigmal;
invVRes=invVprimeRes/Sigmal;
end
%%
function [V,invV,invVRes,condV,invVprime,invVprimeRes,ZData,F]=Fun_NegLogLikelihoodSame(Dl,Respl,Dh,Resph,ZNBC,Thetal,Thetah,Rho,RatioSigmah,phi,nl,nh,nugget,Sigmal,Sigmah,Mu)

[ZData]=TransformData([Respl;Resph],phi,ZNBC);
VprimeDl_Dh=Rho*ComputeRmatrix(Dl,Dh,Thetal);
Vprime=[ ComputeRmatrix2(Dl,Thetal,nugget) ,           VprimeDl_Dh;
    VprimeDl_Dh' , Rho^2*ComputeRmatrix(Dh,Dh,Thetal)+RatioSigmah*ComputeRmatrix2(Dh,Thetah,nugget) ];
V=Vprime*Sigmal;

F=[    ones(nl,1) , zeros(nl,1);
    Rho*ones(nh,1) , ones(nh,1)];
[invVprime,logdetVprime,condV]=invandlogdet(Vprime);
if any(isnan(invVprime),'all')  ||  any(isinf(invVprime),'all')
    Objfval=Inf;Sigmal=[];Thetal=[];Rho=[];RatioSigmah=[];Sigmah=[];Thetah=[];phi=[];Mu=[];V=[];invV=[];invVRes=[];condV=[];invVprime=[];invVprimeRes=[];ZData=[];
    return
end

Res=ZData-F*Mu;

invVprimeRes=invVprime*Res;
invVRes=invVprimeRes/Sigmal;
invV=invVprime/Sigmal;


end

%%
function Fval_ZhQuantile=Fun_ZhQuantile(TeD,Level,Dl,Respl,Dh,Resph,Sigmal,Thetal,Rho,Sigmah,Thetah,phi,ZNBC,Mu,invV,invVRes,nugget)
[ZhPreds,ZhCovs]=Fun_PredictionZ(TeD,Level,Dl,Respl,Dh,Resph,Sigmal,Thetal,Rho,Sigmah,Thetah,phi,ZNBC,Mu,invV,invVRes,nugget);
Fval_ZhQuantile=ZhPreds+norminv(0.9)*ZhCovs.^0.5;
end

%%
function [ZhPreds,ZhVars,ZlPreds,ZlVars,Corrs,CrossCovs]=Fun_PredictionZ(TeD,Level,Dl,Respl,Dh,Resph,Sigmal,Thetal,Rho,Sigmah,Thetah,phi,ZNBC,Mu,invV,invVRes,nugget,ForAEI)
ZhPreds=[];
ZhVars=[];
ZlPreds=[];
ZlVars=[];
Corrs=[];
CrossCovs=[];

if Level==1
    fl=[1 ,0 ] ;
    rlT= [ Sigmal*ComputeRmatrix(TeD,Dl,Thetal) , Rho*Sigmal*ComputeRmatrix(TeD,Dh,Thetal)];
    rlT_invV=rlT*invV;
    ZlPreds=fl*Mu + rlT*invVRes;
    ZlVars=Sigmal*(1+nugget) - sum(rlT_invV.*rlT,2);%%column vector
    ZlVars=max(ZlVars,0);
    
    if(nargin>17)
        %Covariance and correlation coefficient between Zh(TeD)|Z and Zl(TeD)|Z
        fh=[Rho ,1];
        rhT=[ Rho*Sigmal*ComputeRmatrix(TeD,Dl,Thetal) , Rho^2*Sigmal*ComputeRmatrix(TeD,Dh,Thetal)+Sigmah*ComputeRmatrix(TeD,Dh,Thetah) ] ;
        rhT_invV=rhT*invV;
        ZhPreds=fh*Mu + rhT*invVRes;
        ZhVars= Rho^2*Sigmal + Sigmah*(1+nugget) - sum(rhT_invV.*rhT,2);%column vector
        ZhVars=max(ZhVars,0);
        
        CrossCovs=Rho*Sigmal*(1) - sum(rlT_invV.*rhT,2);%column vector
        Corrs=CrossCovs./ ((ZlVars).^0.5) ./ ((ZhVars).^0.5);%column vector
        Corrs(((ZhVars==0) | (ZlVars==0)  ))=0;%column vector
        Corrs=max(min(Corrs,1),-1);
    end
    
end
if Level==2
    fh=[Rho ,1];
    rhT=[ Rho*Sigmal*ComputeRmatrix(TeD,Dl,Thetal) , Rho^2*Sigmal*ComputeRmatrix(TeD,Dh,Thetal)+Sigmah*ComputeRmatrix(TeD,Dh,Thetah) ] ;
    ZhPreds=fh*Mu + rhT*invVRes;
    rhT_invV=rhT*invV;
    ZhVars= Rho^2*Sigmal + Sigmah*(1+nugget) - sum(rhT_invV.*rhT,2);%column vector
    ZhVars=max(ZhVars,0);
end

end
%%
function [NextPoint,NextLevel,AF]=SequentialRun_AEI(Budget,RatioCost,Dl,Respl,Dh,Resph,Sigmal,Thetal,Rho,Sigmah,Thetah,phi,ZNBC,Mu,invV,invVRes,NewAFPoints,Xhat_new,nugget)
Dim=size(Dl,2);
lb=0*ones(1,Dim);ub=1*ones(1,Dim);
options=optimoptions('patternsearch','disp','off');
minZhRef=Fun_PredictionZ(Xhat_new,2,Dl,Respl,Dh,Resph,Sigmal,Thetal,Rho,Sigmah,Thetah,phi,ZNBC,Mu,invV,invVRes,nugget);
if Budget<RatioCost
    
    NextLevel=1;
    Obj_Fun_MinusAEI=@(x) Fun_MinusAEI(x,minZhRef,RatioCost,NextLevel,Dl,Respl,Dh,Resph,Sigmal,Thetal,Rho,Sigmah,Thetah,phi,ZNBC,Mu,invV,invVRes,nugget);
    [fBest,Sortidx]=sort(Obj_Fun_MinusAEI(NewAFPoints));
    parfor id=1:90
        [XBestTry(id,:),fBestTry(id,1)]= patternsearch(Obj_Fun_MinusAEI,NewAFPoints(Sortidx(id),:),[],[],[],[],lb,ub,[],options);
    end
    [fBest,minidx]=min(fBestTry);
    NextPoint=XBestTry(minidx,:) ;
    AF=-fBest;
    
else
    Level1=1;
    Obj_Fun_MinusAEI1=@(x) Fun_MinusAEI(x,minZhRef,RatioCost,Level1,Dl,Respl,Dh,Resph,Sigmal,Thetal,Rho,Sigmah,Thetah,phi,ZNBC,Mu,invV,invVRes,nugget);
    [fval_MinusAEI,fval_MinusEI]=Obj_Fun_MinusAEI1(NewAFPoints);
    [~,Sortidx]=sort(fval_MinusAEI);
    parfor id=1:90
        [XBestTry(id,:),fBestTry(id,1)]= patternsearch(Obj_Fun_MinusAEI1,NewAFPoints(Sortidx(id),:),[],[],[],[],lb,ub,[],options);
    end
    [fval_MinusAEI1,minidx]=min(fBestTry);
    NextPoint1=XBestTry(minidx,:) ;
    
    Level2=2;
    Obj_Fun_MinusAEI2=@(x) Fun_MinusAEI(x,minZhRef,RatioCost,Level2,Dl,Respl,Dh,Resph,Sigmal,Thetal,Rho,Sigmah,Thetah,phi,ZNBC,Mu,invV,invVRes,nugget);
    [~,Sortidx]=sort(fval_MinusEI);
    parfor id=1:90
        [XBestTry(id,:),fBestTry(id,1)]= patternsearch(Obj_Fun_MinusAEI2,NewAFPoints(Sortidx(id),:),[],[],[],[],lb,ub,[],options);
    end
    [fval_MinusAEI2,minidx]=min(fBestTry);
    NextPoint2=XBestTry(minidx,:) ;
    
    if fval_MinusAEI1<= fval_MinusAEI2
        NextPoint=NextPoint1;        NextLevel=1;AF=-fval_MinusAEI1;
    else
        NextPoint=NextPoint2;        NextLevel=2;AF=-fval_MinusAEI2;
    end
end

end
%%
function [fval_MinusAEI,fval_MinusEIs,Corrs,ZhCovs]=Fun_MinusAEI(TeD,MinZhRef,RatioCost,AEILevel,Dl,Respl,Dh,Resph,Sigmal,Thetal,Rho,Sigmah,Thetah,phi,ZNBC,Mu,invV,invVRes,nugget)

if AEILevel==2
    Corrs=1;RatioCost1=1;
    TeDLevel2=2;
    [ZhPreds,ZhCovs]=Fun_PredictionZ(TeD,TeDLevel2,Dl,Respl,Dh,Resph,Sigmal,Thetal,Rho,Sigmah,Thetah,phi,ZNBC,Mu,invV,invVRes,nugget);
    Bias=MinZhRef-ZhPreds;
    SD=ZhCovs.^0.5;
    EIs=Bias.*normcdf(Bias./SD) +SD.*normpdf(Bias./SD);
    fvalAEIs=EIs.*Corrs*RatioCost1;
    fvalAEIs(  [ZhCovs==0] |  (min(pdist2(TeD,Dh),[],2) < 10^(-3))         ) =0;
elseif AEILevel==1
    ForAEI=1;
    TeDLevel1=1;
    [ZhPreds,ZhCovs,~,~,Corrs]=Fun_PredictionZ(TeD,TeDLevel1,Dl,Respl,Dh,Resph,Sigmal,Thetal,Rho,Sigmah,Thetah,phi,ZNBC,Mu,invV,invVRes,nugget,ForAEI);
    Bias=MinZhRef-ZhPreds;
    SD=ZhCovs.^0.5;
    EIs=Bias.*normcdf(Bias./SD) +SD.*normpdf(Bias./SD);
    fvalAEIs=EIs.*Corrs*RatioCost;
    fvalAEIs(  [ min(pdist2(TeD,Dl),[],2) < 10^(-3)] |   [min(pdist2(TeD,Dh),[],2) < 10^(-3)] |  [Corrs < 0]      )=0; %1 and 3 and 4 and 5
end
fval_MinusAEI=-fvalAEIs;
fval_MinusEIs=-EIs;
end