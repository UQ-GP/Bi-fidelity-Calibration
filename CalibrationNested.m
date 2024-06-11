function [RecordTable,RecordData]=CalibrationNested(DataInput)
% Bi-fidelity calibration method : Nested proposed by Pang et al. (2017)
%           Resph: is the square root of HF SSE
%           Respl: is the square root of LF SSE
%           Xhats: is MLE of the calibration parameter vector at all iterations
%           Resphminhats: is posterior mean of the Resph at Xhats at all iterations
%           SSETrue_Xhats: is the true values of SSE evaluated at Xhats at all iterations
% Reference:
% Le Gratiet L, Garnier J (2014) Recursive co-kriging model for design of computer experiments with multiple levels of fidelity. Int J Uncertain Quantif, 4(5):365-386
% Pang, G., Perdikaris, P., Cai, W., & Karniadakis, G. E. (2017). Discovering variable fractional orders of advection-dispersion equations from field data using multi-fidelity Bayesian optimization. Journal of Computational Physics, 348, 694-714.
nugget=1e-6;
%Computes the HF response and the LF response
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
Resph=(sum((Yh-PhysData).^2,2)).^0.5;
Respl=(sum((Yl-PhysData).^2,2)).^0.5;

Resplh=[Respl;Resph];
Budget=Budget-(nl*1+nh*RatioCost);
DlnotDh=Dl(nh+1:nl,:);
ResplnotDh=Respl(nh+1:nl,:);
YlnotDh=Yl(nh+1:nl,:);
if Dim==2
    nlevel=2501;
elseif Dim==3
    nlevel=201;
end
AFPoints=(fullfact(nlevel*ones(1,Dim))-1)/(nlevel-1);

lb=0*ones(1,Dim);ub=1*ones(1,Dim);
options=optimoptions('patternsearch','disp','off');
%Bayesian optimization
ZFit=1;
AFs(n,:)=0;
while (1)
    %Fits the GP model and finds the minimum posterior mean using the recursive formulations by Le Gratiet L and Garnier J (2014)
    if ZFit==1
        [Thetal,Mul,Sigmal,invRl,Objfvall,condRl,invRlRes,PDFl]=GPFitLF(Dl,Respl,nugget);
        [Thetah,Rho,Muh,Sigmah,invRh,Objfvalh,condRh,invRhRes,PDFh]=GPFitHF(Dh,Resph,Dl,Respl,nugget);
        ZFit=0;
        Sigmals(n-1:n,:)=[Sigmal;Sigmal];
        Thetals(n-1:n,:)=[Thetal;Thetal];
        Rhos(n-1:n,:)=[Rho;Rho];
        Sigmahs(n-1:n,:)=[ Sigmah;Sigmah];
        Thetahs(n-1:n,:)=[Thetah;Thetah];
        Mus(n-1:n,:)=[Mul Muh; Mul Muh];
        
        PDFhs(n-1:n,:)=PDFh;
        PDFls(n-1:n,:)=PDFl;
        PDFs(n-1:n,:)=PDFh*PDFl;
        
    else
        [invRl,condRl,invRlRes,FlT,Rl]=Fun_NegLogLikehoodLF_Same(Dl,Respl,Thetal,nugget,Mul,Sigmal);
        PDFl=mvnpdf(Respl,FlT'*Mul,Rl*Sigmal);
        
        [invRh,condRh,invRhRes,Fh,Betah,Rh]=Fun_NegLogLikehoodHF_Same(Dh,Resph,Dl,Respl,Thetah,nugget,Rho,Muh,Sigmah);
        PDFh = mvnpdf(Resph,Fh*Betah,Rh*Sigmah);
    end
    [minResph,minidx]=min(Resph);
    Xhat_new=Dh(minidx,:) ;
    Resphminhat_new=minResph;
    
    %Stores GP model parameters and the MLE of the calibration parameter vector
    
    Xhats(n-1:n,:)=[ Xhat_new; Xhat_new];
    Resphminhats(n-1:n,:)=[ Resphminhat_new; Resphminhat_new];
    
    %Evaluate the true SSE at the MLE of the calibration parameter vector
    Yh_Xhat_new=Yh(minidx,:);
    SSEXhat_new=sum( (Yh_Xhat_new-PhysData).^2);
    Yh_Xhats(n-1:n,:)=[ Yh_Xhat_new ; Yh_Xhat_new];
    SSETrue_Xhats(n-1:n,:)=[SSEXhat_new; SSEXhat_new  ];
    condRlRhs(n-1:n,:)=[condRl,condRh;condRl,condRh];
    disp(['Xhat_new and SSETrue_Xhats at  ' num2str(n) ' -iter '  num2str(Xhats(n,:))  '    ' num2str(SSETrue_Xhats(n,:))  ])
    
    if Budget<1
        break
    end
    
    %Adds the follow-up design points by maximizing the EI for the HF response
    Obj_MinusEI_HF= @(x) Fun_MinusEI_HF(x,Dh,Resph,Dl,Respl,Thetal,Mul,Sigmal,invRl,invRlRes,Thetah,Rho,Muh,Sigmah,invRh,invRhRes,minResph,nugget);
    fvalEIObj=Obj_MinusEI_HF(AFPoints);
    [~,Sortidx]=sort(fvalEIObj);
    parfor id=1:90
        StartPoint=AFPoints(Sortidx(id),:);
        [XBestTry(id,:),fBestTry(id,1)]=patternsearch(Obj_MinusEI_HF,StartPoint,[],[],[],[],lb,ub,[],options);
    end
    [fBest,minidx]=min(fBestTry);
    NextPoint=XBestTry(minidx,:);
    AF=-fBest;
    AFs([n+1 n+2],:)=AF;
    disp(['Current budget=' num2str(Budget) '. The ' num2str(n+1:n+2) '-th run will be at point ' num2str(NextPoint,' %1.3f ')    ])
    
    Dh(nh+1,:)=NextPoint;
    Yh(nh+1,:)=Simulator(NextPoint,2,Case);
    Resph(nh+1,:) =sum((Yh(nh+1,:)-PhysData).^2)^0.5;
    
    Dl=[Dh;DlnotDh];
    Yl=[Yl(1:nh,:);Yh(nh+1,:);YlnotDh];
    Yl(nl+1,:)=Simulator(NextPoint,1,Case);
    Respl_new=sum((Yl(nl+1,:)-PhysData).^2)^0.5;
    Respl=[Respl(1:nh,:); Respl_new; ResplnotDh];
    
    D=[D ; NextPoint; NextPoint ];
    Level=[ Level ; 2 ;1];
    Resplh=[Resplh; Resph(nh+1,:) ; Respl_new];
    
    Budget=Budget-RatioCost-1;
    nh=nh+1;
    nl=nl+1;
    n=n+2;
end
Sigmals(n-1:n,:)=[Sigmal;Sigmal];
Thetals(n-1:n,:)=[Thetal;Thetal];
Rhos(n-1:n,:)=[Rho;Rho];
Sigmahs(n-1:n,:)=[ Sigmah;Sigmah];
Thetahs(n-1:n,:)=[Thetah;Thetah];
Mus(n-1:n,:)=[Mul Muh; Mul Muh];
condRlRhs(n-1:n,:)=[condRl,condRh;condRl,condRh];
PDFhs(n-1:n,:)=PDFh;
PDFls(n-1:n,:)=PDFl;
PDFs(n-1:n,:)=PDFh*PDFl;

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
RecordData.ZNBC=[];
RecordData.Yh_Xhats=Yh_Xhats;

%Stores GP model parameters and the MLE of the calibration parameter vector at all iterations with a table
RecordTable=table(D,Level,Resplh,Sigmals,Thetals,Rhos,Sigmahs,Thetahs,Mus,Xhats,Resphminhats,SSETrue_Xhats,condRlRhs,PDFhs,PDFls,PDFs,AFs);
end
%%
function [Thetal,Mul,Sigmal,invRl,Objfvall,condRl,invRlRes,PDFl]=GPFitLF(Dl,Respl,nugget)

[nl,Dim]=size(Dl);
lb=[(0.25)*ones(1,Dim) ];
ub=[(15)*ones(1,Dim ) ];
nvar=numel(lb) ;
Objl=@(Parl) Fun_NegLogLikehoodLF(Dl,Respl,Parl,nugget);
options=optimoptions('patternsearch','disp','off','MaxFunctionEvaluations',2000*nvar,'MaxIterations',100*nvar);

Sobolset=sobolset(nvar,'Skip',1e3,'Leap',1e2);
StardardPoints=[net(Sobolset,7000*nvar);];

SIZE=size(StardardPoints,1);
fvals=zeros(SIZE,1);
parfor id=1:SIZE
    Point=lb + (ub-lb).*StardardPoints(id,:);
    fvals(id,1)=Objl(Point);
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
    [XBestTry(id,:),fBestTry(id,1)]= patternsearch(Objl,Point,[],[],[],[],lb,ub,[],options); %#ok<PFBNS>
end

[fBest,minidx]=min(fBestTry);
Parl=XBestTry(minidx,:) ;
[Objfvall,Mul,Sigmal,invRl,Thetal,condRl,invRlRes,FlT,Rl]=Objl(Parl);
PDFl=mvnpdf(Respl,FlT'*Mul,Rl*Sigmal);
end

function [Objfvall,Mul,Sigmal,invRl,Thetal,condRl,invRlRes,FlT,Rl]=Fun_NegLogLikehoodLF(Dl,Respl,Thetal,nugget)
[nl,~]=size(Respl);
FlT=ones(1,nl);
Rl=ComputeRmatrix2(Dl,Thetal,nugget);
[invRl, logdetRl,condRl]=invandlogdet(Rl);
FlT_invRl=FlT*invRl;

Mul=(FlT_invRl*Respl)/sum(invRl,'all');
Res=Respl-Mul;
invRlRes=invRl*Res;
Sigmal=Res'*invRlRes/ (nl );
Objfvall=(nl)*log(Sigmal)+ logdetRl;
end
function [invRl,condRl,invRlRes,FlT,Rl]=Fun_NegLogLikehoodLF_Same(Dl,Respl,Thetal,nugget,Mul,Sigmal)
[nl,~]=size(Respl);
FlT=ones(1,nl);
Rl=ComputeRmatrix2(Dl,Thetal,nugget);
[invRl, logdetRl,condRl]=invandlogdet(Rl);
FlT_invRl=FlT*invRl;

Res=Respl-Mul;
invRlRes=invRl*Res;
end

%%
function [Thetah,Rho,Muh,Sigmah,invRh,Objfvalh,condRh,invRhRes,PDFh]=GPFitHF(Dh,Resph,Dl,Respl,nugget)

[nh,Dim]=size(Dh);
lb=[(0.25)*ones(1,Dim) ];
ub=[(15)*ones(1,Dim ) ];
nvar=numel(lb) ;
options=optimoptions('patternsearch','disp','off','MaxFunctionEvaluations',2000*nvar,'MaxIterations',100*nvar);
Objh= @(Parh) Fun_NegLogLikehoodHF(Dh,Resph,Dl,Respl,Parh,nugget);

Sobolset=sobolset(nvar,'Skip',1e3,'Leap',1e2);
StardardPoints=[net(Sobolset,7000*nvar)];

SIZE=size(StardardPoints,1);
fvals=zeros(SIZE,1);
parfor id=1:SIZE
    Point=lb + (ub-lb).*StardardPoints(id,:);
    fvals(id,1)=Objh(Point);
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
    [XBestTry(id,:),fBestTry(id,1)]= patternsearch(Objh,Point,[],[],[],[],lb,ub,[],options);
end

[~,minidx]=min(fBestTry);
Parh=XBestTry(minidx,:) ;
[Objfvalh,Rho,Muh,Sigmah,invRh,Thetah,condRh,invRhRes,Fh,Betah,Rh]=Objh(Parh);
PDFh = mvnpdf(Resph,Fh*Betah,Rh*Sigmah);
end

function [Objfvalh,Rho,Muh,Sigmah,invRh,Thetah,condRh,invRhRes,Fh,Betah,Rh]=Fun_NegLogLikehoodHF(Dh,Resph,Dl,Respl,Thetah,nugget)
[nh,~]=size(Resph);
Fh=[Respl(1:nh,:) ones(nh,1)];
Rh=ComputeRmatrix2(Dh,Thetah,nugget);
[invRh, logdetRh,condRh]=invandlogdet(Rh);
FhT_invRh=Fh'*invRh;

Betah=(FhT_invRh * Fh)\( FhT_invRh *Resph);
Rho=Betah(1);
Muh=Betah(2);
Res=Resph-Fh*Betah;
invRhRes=invRh*Res;
Sigmah=Res'*invRh*Res/ (nh);
Objfvalh=nh*log(Sigmah)+ logdetRh;
end

function [invRh,condRh,invRhRes,Fh,Betah,Rh]=Fun_NegLogLikehoodHF_Same(Dh,Resph,Dl,Respl,Thetah,nugget,Rho,Muh,Sigmah)
[nh,~]=size(Resph);
Fh=[Respl(1:nh,:) ones(nh,1)];
Rh=ComputeRmatrix2(Dh,Thetah,nugget);
[invRh, logdetRh,condRh]=invandlogdet(Rh);

Betah=[Rho,Muh]';
Res=Resph-Fh*Betah;
invRhRes=invRh*Res;
end
%%
function [Predh,Covh,Predl,Covl]=Fun_Prediction(TeDh,Dh,Resph,Dl,Respl,Thetal,Mul,Sigmal,invRl,invRlRes,Thetah,Rho,Muh,Sigmah,invRh,invRhRes,nugget)
[nTest,~]=size(TeDh);
rlT=ComputeRmatrix(TeDh,Dl,Thetal);
Predl=Mul+ rlT*invRlRes;

rhT=ComputeRmatrix(TeDh,Dh,Thetah);
fh=[ Predl ones(nTest,1) ] ;
Betah=[Rho; Muh];
Predh=fh*Betah+rhT*invRhRes;

Covl=Sigmal*( 1+nugget- sum(rlT*invRl.*rlT,2));
Covl2=max(Covl-nugget*Sigmal ,0);
Covl=max(Covl ,0);
Covh=Rho^2*(Covl2)+Sigmah*( 1+nugget- sum(rhT*invRh.*rhT,2));
Covh=max(Covh,0);

end
%%
function Fval_MinusEI_HF=Fun_MinusEI_HF(AFPoints,Dh,Resph,Dl,Respl,Thetal,Mul,Sigmal,invRl,invRlRes,Thetah,Rho,Muh,Sigmah,invRh,invRhRes,minResph,nugget)
[Predh,Covh]=Fun_Prediction(AFPoints,Dh,Resph,Dl,Respl,Thetal,Mul,Sigmal,invRl,invRlRes,Thetah,Rho,Muh,Sigmah,invRh,invRhRes,nugget);
SD=Covh.^0.5;
Bias=minResph-Predh;
fvalEI=Bias.*normcdf(Bias./SD)+SD.*normpdf(Bias./SD);
fvalEI( (Covh==0) |  (min(pdist2(AFPoints,Dh),[],2) < 10^(-3)))=0;
Fval_MinusEI_HF=-fvalEI;
end