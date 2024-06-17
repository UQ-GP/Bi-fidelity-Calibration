%Section 1: Sets parameters for all calibration methods 
clc,clear,format compact
Dim=2;
Case=3;
nl=12;
nh=4;
nh0=8;
RatioCost=3;
InitialBudget=nl*1+nh*RatioCost;
InitialBudget0=nh0*RatioCost;
Budget=InitialBudget+12;
load Example3Size2.mat MultiDataInput SingleDataInput XTrue yh_XTrue PhysData
% XTrue=[0.1 0.4];
% [yh_XTrue]= Simulator(XTrue,2,Case);
% std_error=(var(yh_XTrue)*0.0001)^0.5;
% PhysData=yh_XTrue+normrnd(0,std_error,size(yh_XTrue));
SSE_XTrue=sum([Simulator(XTrue,2,Case)-PhysData].^2);

[X1,X2]=meshgrid(linspace(0,1,51)');
TestPoints= [X1(:) X2(:)];
for id=1:size(TestPoints,1)
    TrueSh(id,1)=sum((Simulator(TestPoints(id,:),2,Case)-PhysData).^2); 
end
[~,sortidx]=sort(TrueSh);

lb=0*ones(1,Dim);
ub=1*ones(1,Dim);
options=optimoptions('patternsearch','MaxIterations',10^6,'MeshTolerance',10^-6,'TolFun',10^-8,'TolX',10^-8,'MaxFunEvals',10^8);

SSHFun=@(x) sum([Simulator(x,2,Case)-PhysData].^2);
for id=1:50
    StartPoint= TestPoints(sortidx(id),:);
    [XMLETry(id,:),fval(id,:)]=patternsearch(SSHFun,StartPoint,[],[],[],[],lb,ub,[],options)  ;
end
[~,minidx]=min(fval);
XMLE=XMLETry(minidx,:);
SSE_XMLE=min(fval);
%{
parfor id=1:100
    id
    [Dl,Dh]=GenerateNestedLHD(nl,nh,Dim,1e5);     
    [Dh0]=GenerateNestedLHD(nh0,nh0,Dim,1e5);     
    
    Dls(:,:,id)=Dl;
    Dhs(:,:,id)=Dh;
    Dh0s(:,:,id)=Dh0;    
end


for id=1:100
    id

    Dl=Dls(:,:,id);
    Dh=Dhs(:,:,id);
    Dh0=Dh0s(:,:,id);
    
    clear Yl Yh
    for jd=1:nl
        Yl(jd,:)=Simulator(Dl(jd,:),1,Case);
    end
    for jd=1:nh
        Yh(jd,:)=Simulator(Dh(jd,:),2,Case);
    end
    clear  Yh0
    for jd=1:nh0
        Yh0(jd,:)=Simulator(Dh0(jd,:),2,Case);
    end
    
    MultiDataInput(id).Dl=Dl;       MultiDataInput(id).Yl=Yl;
    MultiDataInput(id).Dh= Dh;    MultiDataInput(id).Yh=Yh;
    MultiDataInput(id).XTrue=XTrue;
    MultiDataInput(id).PhysData=PhysData;    MultiDataInput(id).RatioCost=RatioCost;
    MultiDataInput(id).Budget=Budget;           MultiDataInput(id).Case=Case;
    
    SingleDataInput(id).Dl =[] ;       SingleDataInput(id).Yl=[];
    SingleDataInput(id).Dh= Dh0;    SingleDataInput(id).Yh=Yh0;
    SingleDataInput(id).XTrue=XTrue;
    SingleDataInput(id).PhysData=PhysData;    SingleDataInput(id).RatioCost=RatioCost;
    SingleDataInput(id).Budget=Budget;          SingleDataInput(id).Case=Case;
end
%}



















































































%%
%Section 2: Bayesian optimization

ZNBC_BC=1;   ZNBC_ID=0;   ZNBC_SR=2;
ZMLFSSE=1;   ZLFSSE=0;

for id=1:100
    disp('---')
    
    [T_MBC_AGP{id,1}] =CalibrationAGP(MultiDataInput(id),ZNBC_BC,ZMLFSSE); 'MBC-AGP'
    [T_BC_AGP{id,1}] =CalibrationAGP(MultiDataInput(id),ZNBC_BC,ZLFSSE); 'BC-AGP'
    [T_MID_AGP{id,1}] =CalibrationAGP(MultiDataInput(id),ZNBC_ID,ZMLFSSE); 'MID-AGP'
    [T_SR_AGP{id,1}] =CalibrationAGP(MultiDataInput(id),ZNBC_SR,ZLFSSE); 'SR-AGP'
    [T_Nested{id,1}] =CalibrationNested(MultiDataInput(id)); 'Nested'
    [T_SVDAGP{id,1}] =CalibrationSVDAGP(MultiDataInput(id),0.99);'SVD-AGP'
    [T_BC_GP{id,1}] =CalibrationBCGP(SingleDataInput(id)); 'BC-GP'
    [T_SR_GP{id,1}] =CalibrationSRGP(SingleDataInput(id)); 'SR-GP'
    [T_SVD{id,1}] =CalibrationSVD(SingleDataInput(id),0.99);    'SVD'
    save Example3Size2.mat
end
%%
%Section 3: Show BO results
clc,clear
load('Example3Size2.mat')
idx=(1:100);
BORecordTable=[T_MBC_AGP(idx)  T_BC_AGP(idx)   T_MID_AGP(idx)  T_SR_AGP(idx)   T_Nested(idx) T_SVDAGP(idx)   T_BC_GP(idx)  T_SR_GP(idx)  T_SVD(idx)    ];
Labels={'MBC-AGP','BC-AGP','MID-AGP','SR-AGP', 'Nested','SVD-AGP', 'BC-GP','SR-GP' ,'SVD'}' ;

for Trainidx=1:size(BORecordTable,1)
    for Methodidx=1:9
        Table=BORecordTable{Trainidx,Methodidx} ;
        
        DiffSSETrue_XhatsEnd(Trainidx,Methodidx)=Table.SSETrue_Xhats(end,:)-SSE_XMLE;
        SSETrue_XhatsEnd(Trainidx,Methodidx)=Table.SSETrue_Xhats(end,:);
        XhatsEnd=Table.Xhats(end,:);
        L2End(Trainidx,Methodidx)=norm(XhatsEnd-XMLE);
        if Methodidx<=2 || Methodidx==7
            phiEnd(Trainidx,Methodidx)=Table.phis(end,:);
        end
        
        costs=[1 RatioCost]';
        SSETrue_Xhats_iter=Table.SSETrue_Xhats;
        Xhats_iter=Table.Xhats;
        L2s_iter=sum((Xhats_iter-XMLE).^2,2).^0.5;
        Level_iter=Table.Level;
        Budget_iter=cumsum(costs(Level_iter));
        
        TrueSSE_Xhats_Budget(1:Budget,Methodidx,Trainidx) = interp1(Budget_iter,SSETrue_Xhats_iter,1:Budget);
        
        L2_Budget(1:Budget,Methodidx,Trainidx)=interp1(Budget_iter,L2s_iter,1:Budget);
        
        if Methodidx==5 || Methodidx==6
            deleteLFidx=(nl+nh+1):2:size(Table,1);
            Budget_iter(deleteLFidx,:)=[];
            SSETrue_Xhats_iter(deleteLFidx,:)=[];
            L2s_iter(deleteLFidx,:)=[];
            
            TrueSSE_Xhats_Budget(:,Methodidx,Trainidx) = interp1(Budget_iter,SSETrue_Xhats_iter,1:(Budget));
            L2_Budget(:,Methodidx,Trainidx)=interp1(Budget_iter,L2s_iter,1:(Budget));
        end
        
    end
end

meanTrueSSE_Xhats_Budget_SSEXMLE=mean(TrueSSE_Xhats_Budget,3)-SSE_XMLE;
meanL2_Budget=mean(L2_Budget,3);

idx1=1;
for idx2=1:9
    [ ~, ttest_p_Sh(idx2,1)]=ttest(SSETrue_XhatsEnd(:,idx1),SSETrue_XhatsEnd(:,idx2));
    [ ~, ttest_p_L2(idx2,1)]=ttest(L2End(:,idx1),L2End(:,idx2));
end

ExampleEnvTable =table(Labels,mean(SSETrue_XhatsEnd)'-SSE_XMLE,ttest_p_Sh,mean(L2End)',ttest_p_L2)










































FontSize=24;
figure,clf
tiledlayout(1,2,'Padding','none','TileSpacing','none');
nexttile
htmlGray = [128 128 128]/255;
htmlGreen = [0.4660 0.6740 0.1880];

MarkerSize=15;
linewidth=4;
plot(1:Budget,meanTrueSSE_Xhats_Budget_SSEXMLE(1:Budget,1),'ko-','linewidth',linewidth,'MarkerSize',MarkerSize,'MarkerIndices',[InitialBudget:3:Budget Budget]),hold on
plot(1:Budget,meanTrueSSE_Xhats_Budget_SSEXMLE(1:Budget,2),'b:o','linewidth',linewidth,'MarkerSize',MarkerSize,'MarkerFaceColor','b','MarkerIndices',[InitialBudget (InitialBudget+2):2:Budget Budget]),
plot(1:Budget,meanTrueSSE_Xhats_Budget_SSEXMLE(1:Budget,3),'k^-','linewidth',linewidth,'MarkerSize',MarkerSize,'MarkerIndices',[InitialBudget (InitialBudget):3:Budget Budget])
plot(1:Budget,meanTrueSSE_Xhats_Budget_SSEXMLE(1:Budget,4),'--v','linewidth',linewidth,'Color', htmlGray,'MarkerSize',MarkerSize,'MarkerIndices',[(InitialBudget+1):2:Budget Budget]),hold on
plot(1:Budget,meanTrueSSE_Xhats_Budget_SSEXMLE(1:Budget,5),':s','linewidth',linewidth,'color',htmlGreen,'MarkerFaceColor',htmlGreen,'MarkerSize',MarkerSize,'MarkerIndices',[InitialBudget:4:Budget ])
plot(1:Budget,meanTrueSSE_Xhats_Budget_SSEXMLE(1:Budget,6),'b-x','linewidth',linewidth,'MarkerSize',MarkerSize+10,'MarkerIndices',[InitialBudget (InitialBudget+1):3:Budget Budget])
plot(1:Budget,meanTrueSSE_Xhats_Budget_SSEXMLE(1:Budget,7),':s','linewidth',linewidth,'Color', 'r','MarkerSize',MarkerSize,'MarkerIndices',[InitialBudget (InitialBudget+2):3:Budget Budget ]),hold on
plot(1:Budget,meanTrueSSE_Xhats_Budget_SSEXMLE(1:Budget,8),'--h','linewidth',linewidth,'MarkerFaceColor','none','MarkerSize',MarkerSize,'MarkerIndices',[InitialBudget (InitialBudget+1):3:Budget Budget ]),hold on
plot(1:Budget,meanTrueSSE_Xhats_Budget_SSEXMLE(1:Budget,9),':d','linewidth',linewidth,'MarkerSize',MarkerSize,'MarkerIndices',[InitialBudget (InitialBudget+3):3:Budget Budget ]),hold on

ylim([0.5 200])
set(gca,'FontWeight','bold','FontSize',FontSize,'YScale','log')
xlim([InitialBudget,Budget])
xlabel('Computational cost','FontWeight','normal')
ylabel('Average  $S_h(\hat{\textbf{x}}^*_{\mathbf{ML}})-0.093939$','Interpreter','latex','FontSize',32);
leg = legend(Labels,'NumColumns',3,'Location','northeast');
leg.ItemTokenSize = [74,50];
title('(a)','FontWeight','bold')
yticks([0.3 10.^[0:2] ])
yticklabels({'0.3 ', '10^0 ','10^1 ','10^2 '})
 set(gca,'TickLabelInterpreter', 'tex');
 
nexttile
plot(1:Budget,meanL2_Budget(1:Budget,1),'ko-','linewidth',linewidth,'MarkerSize',MarkerSize,'MarkerIndices',[InitialBudget:3:Budget Budget]),hold on
plot(1:Budget,meanL2_Budget(1:Budget,2),'b:o','linewidth',linewidth,'MarkerSize',MarkerSize,'MarkerFaceColor','b','MarkerIndices',[InitialBudget (InitialBudget):2:Budget Budget]),
plot(1:Budget,meanL2_Budget(1:Budget,3),'k^-','linewidth',linewidth,'MarkerSize',MarkerSize,'MarkerIndices',[InitialBudget (InitialBudget):3:Budget Budget])
plot(1:Budget,meanL2_Budget(1:Budget,4),'--v','linewidth',linewidth,'Color', htmlGray,'MarkerSize',MarkerSize,'MarkerIndices',[(InitialBudget+1):2:Budget Budget]),hold on
plot(1:Budget,meanL2_Budget(1:Budget,5),':s','linewidth',linewidth,'color',htmlGreen,'MarkerFaceColor',htmlGreen,'MarkerSize',MarkerSize,'MarkerIndices',[InitialBudget:4:Budget Budget])
plot(1:Budget,meanL2_Budget(1:Budget,6),'b-x','linewidth',linewidth,'MarkerSize',MarkerSize+10,'MarkerIndices',[InitialBudget (InitialBudget+1):3:Budget Budget])
plot(1:Budget,meanL2_Budget(1:Budget,7),':s','linewidth',linewidth,'Color', 'r','MarkerSize',MarkerSize,'MarkerIndices',[InitialBudget (InitialBudget+2):3:Budget Budget ]),hold on
plot(1:Budget,meanL2_Budget(1:Budget,8),'--h','linewidth',linewidth,'MarkerSize',MarkerSize,'MarkerIndices',[InitialBudget (InitialBudget+1):3:Budget Budget ]),hold on
plot(1:Budget,meanL2_Budget(1:Budget,9),':d','linewidth',linewidth,'MarkerSize',MarkerSize,'MarkerIndices',[InitialBudget (InitialBudget+3):3:Budget Budget ]),hold on
xlim([InitialBudget,Budget])
ylim([0.019    0.56])

xlabel('Computational cost','FontWeight','normal')
title('(b)','FontWeight','bold')
set(gca,'FontWeight','bold','FontSize',FontSize)
ylabel('Average  $L_2(\hat{\textbf{x}}^*_{\mathbf{ML}})$','Interpreter','latex','FontSize',32);
leg = legend(Labels,'NumColumns',3,'Location','northeast');
leg.ItemTokenSize = [74,50];
set(findobj(gcf,'type','axes'),'FontWeight','Bold', 'LineWidth', 3);
set(gca,'YScale','log')
yticks(0.01*2.^[1:5])
yticklabels({'0.02 ','0.04 ','0.08 ','0.16 ','0.32 '})
set(gcf,'Position',[          0         100        1920         615])





















































