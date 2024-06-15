%%%%%%%% %Section 1: Sets parameters for all calibration methods 
clc,clear,format compact 
Dim=3;
Case=2;
nl=20;
nh=10;
nh0=15;
RatioCost=4;
InitialBudget=nl*1+nh*RatioCost
InitialBudget0=nh0*RatioCost
Budget=InitialBudget0+20
Budget1=Budget
load Example2Size2.mat  Dls Dhs Dh0s XTrue yh_XTrue PhysData SSE_XTrue XMLE SSE_XMLE  MultiDataInput SingleDataInput
%{
XTrue=rand(1,Dim);
yh_XTrue= Simulator(XTrue,2,Case)
std_error=(var(yh_XTrue)*0.0001)^0.5;
PhysData=yh_XTrue+normrnd(0,std_error,size(yh_XTrue));
SSE_XTrue=sum([yh_XTrue-PhysData].^2)

lb=0*ones(1,Dim);ub=1*ones(1,Dim);
options=optimoptions('patternsearch','MaxIterations',10^6,'MeshTolerance',10^-3,'MaxFunEvals',10^8);
SSHFun=@(x) sum((Simulator(x,2,Case)-PhysData).^2);
[XMLE,SSE_XMLE]=patternsearch(SSHFun,XTrue,[],[],[],[],lb,ub,[],options)

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
    parfor jd=1:nl
        tic
        Yl(jd,:)=Simulator(Dl(jd,:),1,Case);
        timel(id,jd)=toc;
    end
    parfor jd=1:nh
        tic
        Yh(jd,:)=Simulator(Dh(jd,:),2,Case);
        timeh(id,jd)=toc;
    end
    clear  Yh0
    parfor jd=1:nh0
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
    SingleDataInput(id).Budget=Budget1;          SingleDataInput(id).Case=Case;
    
    
end
%}
%%
%Section 2: Bayesian optimization
ZNBC_BC=1;   ZNBC_ID=0;   ZNBC_SR=2;
ZMLFSSE=1;   ZLFSSE=0;
for id=1:100
    id
    T_MBC_AGP{id,1} =CalibrationAGP(MultiDataInput(id),ZNBC_BC,ZMLFSSE); 'MBC-AGP'
    T_BC_AGP{id,1} =CalibrationAGP(MultiDataInput(id),ZNBC_BC,ZLFSSE); 'BC-AGP'
    T_MID_AGP{id,1} =CalibrationAGP(MultiDataInput(id),ZNBC_ID,ZMLFSSE); 'MID-AGP'
    T_SR_AGP{id,1} =CalibrationAGP(MultiDataInput(id),ZNBC_SR,ZLFSSE); 'SR-AGP'
    T_Nested{id,1} =CalibrationNested(MultiDataInput(id)); 'Nested'
    T_SVDAGP{id,1} =CalibrationSVDAGP(MultiDataInput(id),0.99);'SVD-AGP'
    T_BC_GP{id,1} =CalibrationBCGP(SingleDataInput(id)); 'BC-GP'
    T_SR_GP{id,1} =CalibrationSRGP(SingleDataInput(id)); 'SR-GP'
    T_SVD{id,1} =CalibrationSVD(SingleDataInput(id),0.99);'SVD'
    save Example2Size2.mat
end
%%
%%%%%% %Section 3: Show BO results
clear,format compact,clc;
load Example2Size2.mat
idx=[1:100];
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
        
        if Methodidx==5 || Methodidx==6%Nested
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
Labels1={'(i) vs (i) ','(i) vs BC-AGP ','(i) vs ID-AGP ','(i) vs SR-AGP ' ' (i) vs Nested' ,'(i) vs SVD-AGP',' (i) vs BC-GP',' (i) vs SR-GP', '(i) vs SVD'}';

Table3 =table(Labels,mean(DiffSSETrue_XhatsEnd)',ttest_p_Sh,mean(L2End)',ttest_p_L2)

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
ylim([  1000     3000000])
set(gca,'FontWeight','bold','FontSize',FontSize,'YScale','log')
xlabel('Computational cost','FontWeight','normal')
ylabel('Average $S_h(\hat{\textbf{x}}^*_{\mathbf{ML}})-2553$','Interpreter','latex','FontSize',32);
leg = legend(Labels,'NumColumns',3,'Location','northeast');
leg.ItemTokenSize = [74,50];
title('(a)','FontWeight','bold')
yticks([10.^[0:8] ])
 set(gca,'TickLabelInterpreter', 'tex');
 xticks(InitialBudget:4:Budget)
 xlim([InitialBudget-0.1,Budget+0.1])

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
xlim([InitialBudget-0.1,Budget+0.1])
ylim( [0.2400    0.66000 ])
xlabel('Computational cost','FontWeight','normal')
title('(b)','FontWeight','bold')
set(gca,'FontWeight','bold','FontSize',FontSize)
ylabel('Average  $L_2(\hat{\textbf{x}}^*_{\mathbf{ML}})$','Interpreter','latex','FontSize',32);
leg = legend(Labels,'NumColumns',3,'Location','northeast');
leg.ItemTokenSize = [74,50];
set(findobj(gcf,'type','axes'),'FontWeight','Bold', 'LineWidth', 3);
 xticks(InitialBudget:4:Budget)
set(gcf,'Position',[          0         100        1920         615])










































































































































