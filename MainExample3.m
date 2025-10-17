%Example 3
%Section 1: Sets input data for all model calibration methods.
%Dl is the initial LF design (an nl x d matrix) and Dh is the initial HF design (an nh x d matrix) for each bi-fidelity method in a trial;
%Yl is the initial LF output data (an nl x N matrix) and Yh is the initial HF output data (an nh x N matrix) for each bi-fidelity method in a trial; 
%Dh0 is the initial HF design (an nh0 x d matrix) for each single-fidelity method in a trial;
%Yh0 is the initial HF output data (an nh0 x N matrix) for each single-fidelity method in a trial;
%CostRatio is c_h/c_l; 
%Budget is the total budget for each bi-fidelity method in a trial; 
%Budget is also the total budget for each single-fidelity method in a trial;
%w is the vector of field/physical data (a 1 x N vector).

%%Uncomment the commands below in this section if you want to generate a different set of 100 initial designs for each method and a different vector 
%of field data, rather than use those stored in the .mat file loaded by this script, to perform another set of 100 trials.
%{ 
clear all;clc,format compact 
d=2;
Case=3;
nl=18;
nh=6;
nh0=12;
CostRatio=3;
InitialBudget=nl*1+nh*CostRatio;
InitialBudget0=nh0*CostRatio;
if(InitialBudget~=InitialBudget0)
    return
end
Budget=InitialBudget0+12;
xstar=[0.1 0.4];
yhxstar=Simulator(xstar,2,Case);
stddev=(var(yhxstar)*0.0001)^0.5;
w=yhxstar+normrnd(0,stddev,size(yhxstar));
Shxstar=sum([yhxstar-w].^2);

parfor id=1:100
    id
    [Dl,Dh]=GenerateNestedLHD(nl,nh,d,1e5);     
    [Dh0]=GenerateNestedLHD(nh0,nh0,d,1e5);     
    
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
    clear Yh0
    for jd=1:nh0
        Yh0(jd,:)=Simulator(Dh0(jd,:),2,Case);
    end
    
    MultiDataInput(id).Dl=Dl;               MultiDataInput(id).Yl=Yl;
    MultiDataInput(id).Dh=Dh;               MultiDataInput(id).Yh=Yh;
    MultiDataInput(id).xstar=xstar;
    MultiDataInput(id).w=w;                 MultiDataInput(id).CostRatio=CostRatio;
    MultiDataInput(id).Budget=Budget;       MultiDataInput(id).Case=Case;
    
    SingleDataInput(id).Dl=[];              SingleDataInput(id).Yl=[];
    SingleDataInput(id).Dh=Dh0;             SingleDataInput(id).Yh=Yh0;
    SingleDataInput(id).xstar=xstar;
    SingleDataInput(id).w=w;                SingleDataInput(id).CostRatio=CostRatio;
    SingleDataInput(id).Budget=Budget;      SingleDataInput(id).Case=Case;

end
%}
%save Example3InputData.mat
%clear all
%load Example3InputData.mat xstar yhxstar w Shxstar MultiDataInput SingleDataInput   
%save Example3InputData.mat  
clear all;clc,format compact 
d=2;
Case=3;
nl=18;
nh=6;
nh0=12;
CostRatio=3;
InitialBudget=nl*1+nh*CostRatio;
InitialBudget0=nh0*CostRatio;
if(InitialBudget~=InitialBudget0)
    return
end
Budget=InitialBudget0+12;
load Example3InputData.mat xstar yhxstar w Shxstar MultiDataInput SingleDataInput

[X1,X2]=meshgrid(linspace(0,1,51));
GridPoints=[X1(:) X2(:)];
ShFun=@(x) sum((Simulator(x,2,Case)-w).^2);
NGP=size(GridPoints,1);
ShVals=zeros(NGP,1);
for id=1:NGP
    ShVals(id,1)=ShFun(GridPoints(id,:)); 
end
[~,sortidx]=sort(ShVals);
lb=0*ones(1,d); ub=1*ones(1,d);
options=optimoptions('patternsearch','Display','off');
xstarMLTry=zeros(50,d); ShFunBestValTry=zeros(50,1); exitflagTry=zeros(50,1);
for id=1:50
    StartPoint=GridPoints(sortidx(id),:);
    [xstarMLTry(id,:),ShFunBestValTry(id,:),exitflagTry(id)]=patternsearch(ShFun,StartPoint,[],[],[],[],lb,ub,[],options);
end
[ShxstarML,minidx]=min(ShFunBestValTry);
xstarML=xstarMLTry(minidx,:);

idxTrial=45;
YlDh=MultiDataInput(idxTrial).Yl(1:nh,:);
Yh=MultiDataInput(idxTrial).Yh;

N=numel(w);
OnesVec=ones(nh,1);
ahati_bhati=zeros(2,N);
for kd=1:N
    ModelMatrix=[OnesVec,YlDh(:,kd)];
    lastwarn('');
    ahati_bhati(:,kd)=regress(Yh(:,kd),ModelMatrix);
    [warnMsg,~]=lastwarn;
    if contains(warnMsg,'X is rank deficient to within machine precision.')
        return
    end
end

[X1,X2]=meshgrid(linspace(0,1,501));
ModifiedLFSSEVals=zeros(501,501); LFSSEVals=zeros(501,501); HFSSEVals=zeros(501,501);
for id=1:501
    for jd=1:501
        yl0=Simulator([X1(id,jd),X2(id,jd)],1,Case);
        yh0=Simulator([X1(id,jd),X2(id,jd)],2,Case);
        YlModifiedGrid=ahati_bhati(1,:)+yl0.*ahati_bhati(2,:);
        
        ModifiedLFSSEVals(id,jd)=sum((YlModifiedGrid-w).^2);         
        LFSSEVals(id,jd)=sum((yl0-w).^2); 
        HFSSEVals(id,jd)=sum((yh0-w).^2); 
    end
end

Levels=[3 10 25 50 100 250 500 1000 1500 2.5e3 6e3 12e3 24e3 40e3];
FontSize0=32;
FontSizeLevel=30;
figure,clf
tiledlayout(1,3,'Padding','none','TileSpacing','none');
nexttile
[C,h]=contour(X1,X2,LFSSEVals,Levels,'TextStep',2,'linewidth',4);
clabel(C,h,'LabelSpacing',155,'FontWeight','bold','FontSize',FontSizeLevel,'Color','k','linewidth',2)
text(0.47,-0.17,'x_1','FontSize',FontSize0,'FontWeight','Bold')
ylabel('x_2','FontSize',FontSize0,'Rotation',0,'HorizontalAlignment','right')
title('(a)','FontSize',FontSize0,'FontWeight','Bold')
xticks([0:0.2:1])
yticks([0:0.2:1])
set(gca,'FontWeight','bold','FontSize',FontSize0)
grid on
 
nexttile
[C,h]=contour(X1,X2,HFSSEVals,Levels,'linewidth',4);
clabel(C,h,'LabelSpacing',200,'FontWeight','bold','FontSize',FontSizeLevel,'Color','k','linewidth',2)
text(0.47,-0.17,'x_1','FontSize',FontSize0,'FontWeight','Bold')
ylabel('x_2','FontSize',FontSize0,'Rotation',0,'HorizontalAlignment','right')
title('(b)','FontSize',FontSize0,'FontWeight','Bold')
xticks([0:0.2:1])
yticks([0:0.2:1])
set(gca,'FontWeight','bold','FontSize',FontSize0)
grid on

nexttile
[C,h]=contour(X1,X2,ModifiedLFSSEVals,Levels,'linewidth',4);
clabel(C,h,'LabelSpacing',190,'FontWeight','bold','FontSize',FontSizeLevel,'Color','k','linewidth',2)
xlabel(' ','FontSize',FontSize0)
text(0.47,-0.17,'x_1','FontSize',FontSize0,'FontWeight','Bold')
ylabel('x_2','FontSize',FontSize0,'Rotation',0,'HorizontalAlignment','right')
title('(c)','FontSize',FontSize0,'FontWeight','Bold')
xticks([0:0.2:1])
yticks([0:0.2:1])
set(gca,'FontWeight','bold','FontSize',FontSize0)
grid on
set(findobj(gca,'type','line'),'LineWidth',4)
set(gcf,'position',[0 150 1886 631])
set(findobj(gcf,'type','axes'),'FontWeight','Bold','LineWidth',3); 
%%
save Example3InputData.mat
clear all
load Example3InputData.mat xstar yhxstar w Shxstar xstarML ShxstarML MultiDataInput SingleDataInput    
save Example3InputData.mat
clear all
d=2;
Case=3;
nl=18;
nh=6;
nh0=12;
CostRatio=3;
InitialBudget=nl*1+nh*CostRatio;
InitialBudget0=nh0*CostRatio;
if(InitialBudget~=InitialBudget0)
    return
end
Budget=InitialBudget0+12;
load Example3InputData.mat xstar yhxstar w Shxstar xstarML ShxstarML MultiDataInput SingleDataInput    
%%
%Section 2: Runs all model calibration methods.
Z_BC=1;    Z_ID=0;   Z_SR=2;
ZMLFSSE=1; ZLFSSE=0; AccuracyLevel=1; t=0.99; 
for id=1:100
    id
    [T_MBC_AGP{id,1},~,RunTime_MBC_AGP(id)]=CalibrationAGP(MultiDataInput(id),Z_BC,ZMLFSSE,AccuracyLevel); 'MBC-AGP'
    [T_BC_AGP{id,1},~,RunTime_BC_AGP(id)]=CalibrationAGP(MultiDataInput(id),Z_BC,ZLFSSE,AccuracyLevel); 'BC-AGP'
    [T_MID_AGP{id,1},~,RunTime_MID_AGP(id)]=CalibrationAGP(MultiDataInput(id),Z_ID,ZMLFSSE,AccuracyLevel); 'MID-AGP'
    [T_SR_AGP{id,1},~,RunTime_SR_AGP(id)]=CalibrationAGP(MultiDataInput(id),Z_SR,ZLFSSE,AccuracyLevel); 'SR-AGP'
    [T_Nested{id,1},~,RunTime_Nested(id)]=CalibrationNested(MultiDataInput(id),AccuracyLevel); 'Nested'
    [T_SVD_AGP{id,1},~,RunTime_SVD_AGP(id)]=CalibrationSVDAGP(MultiDataInput(id),AccuracyLevel,t); 'SVD-AGP'
    [T_BC_GP{id,1},~,RunTime_BC_GP(id)]=CalibrationBCGP(SingleDataInput(id),AccuracyLevel); 'BC-GP'
    [T_SR_GP{id,1},~,RunTime_SR_GP(id)]=CalibrationSRGP(SingleDataInput(id),AccuracyLevel); 'SR-GP'
    [T_SVD{id,1},~,RunTime_SVD(id)]=CalibrationSVD(SingleDataInput(id),AccuracyLevel,t); 'SVD'
    save Example3.mat
end
%%
%Section 3: Constructs figures and a table to illustrate results.
clear all;clc,format compact 
load Example3.mat
idx=(1:100);
RecordTable=[T_MBC_AGP(idx) T_BC_AGP(idx) T_MID_AGP(idx) T_SR_AGP(idx) T_Nested(idx) T_SVD_AGP(idx) T_BC_GP(idx) T_SR_GP(idx) T_SVD(idx)];
Labels={'MBC-AGP','BC-AGP','MID-AGP','SR-AGP','Nested','SVD-AGP','BC-GP','SR-GP','SVD'}';

for idxMethod=1:9
    
    for idxTrial=1:numel(idx)
        Table=RecordTable{idxTrial,idxMethod};
        
        if idxMethod<=2 || idxMethod==7
            phiEnd(idxTrial,idxMethod)=Table.phis(end,:);
        end
        
        costs=[1 CostRatio]';
        ShxhatstarMLs=Table.ShxhatstarMLs;
        xhatstarMLs=Table.xhatstarMLs;
        
        L2xhatstarMLs=sum((xhatstarMLs-xstarML).^2,2).^0.5;
        Levels=Table.Level;
        Costs=cumsum(costs(Levels));

        ShxhatstarMLsEnd(idxTrial,idxMethod)=ShxhatstarMLs(end,:);
        L2xhatstarMLsEnd(idxTrial,idxMethod)=L2xhatstarMLs(end);
        ShxhatstarMLsEndminusShxstarML(idxTrial,idxMethod)=ShxhatstarMLs(end,:)-ShxstarML;
        
        if idxMethod~=5 && idxMethod~=6
        InterpolatedShxhatstarMLs(1:Budget,idxMethod,idxTrial)=interp1(Costs,ShxhatstarMLs,1:Budget);
        
        InterpolatedL2xhatstarMLs(1:Budget,idxMethod,idxTrial)=interp1(Costs,L2xhatstarMLs,1:Budget);
        
        elseif idxMethod==5 || idxMethod==6
            deleteidx=(nl+nh+1):2:size(Table,1);
            Costs(deleteidx,:)=[];
            ShxhatstarMLs(deleteidx,:)=[];
            L2xhatstarMLs(deleteidx,:)=[];
            
            InterpolatedShxhatstarMLs(1:Budget,idxMethod,idxTrial)=interp1(Costs,ShxhatstarMLs,1:Budget);
            InterpolatedL2xhatstarMLs(1:Budget,idxMethod,idxTrial)=interp1(Costs,L2xhatstarMLs,1:Budget);
        end
        
    end
end
meanInterpolatedShxhatstarMLs=mean(InterpolatedShxhatstarMLs,3);
meanInterpolatedShxhatstarMLsminusShxstarML=meanInterpolatedShxhatstarMLs-ShxstarML;
meanInterpolatedL2xhatstarMLs=mean(InterpolatedL2xhatstarMLs,3);

idx1=1;
for idx2=1:9
    [~,ttest_pval_Sh(idx2,1)]=ttest(ShxhatstarMLsEnd(:,idx1),ShxhatstarMLsEnd(:,idx2));
    [~,ttest_pval_L2(idx2,1)]=ttest(L2xhatstarMLsEnd(:,idx1),L2xhatstarMLsEnd(:,idx2));
end
AverageSh=mean(ShxhatstarMLsEnd)'; AverageL2=mean(L2xhatstarMLsEnd)';
Table3=table(Labels,AverageSh,ttest_pval_Sh,AverageL2,ttest_pval_L2)

Labels2={'MBC-AGP','\newline BC-AGP','MID-AGP','\newline SR-AGP','Nested','\newline SVD-AGP','BC-GP','\newline SR-GP','SVD'}';
figure,clf
subplot(121)
boxplot(ShxhatstarMLsEndminusShxstarML,'Labels',Labels2)
bp=gca; bp.FontSize=20;
bpXAxisFontSize=21;
bp.XAxis.FontWeight='bold'; bp.XAxis.FontSize=bpXAxisFontSize;
bp.YAxis.FontWeight='bold'; bp.YAxis.FontSize=23;
ylim([0.00005,60])
set(gca,'YScale','log')
ylabel('$S_h(\hat{\textbf{x}}^*_{\mathbf{ML}})-$0.093939','Interpreter','latex','FontSize',28);
set(findobj(gca,'type','line'),'LineWidth',2)
title('(a)','FontSize',25,'FontWeight','bold')
set(gca,'Position',[0.07 0.155 0.42 0.765])
set(gca,'yGrid','on','GridLineStyle','--')
yticks([10.^[-4:1] 50])
yticklabels({'10^{-4}','10^{-3}','10^{-2}','10^{-1}','10^{0}','10^{1}','50'})
set(gca,'TickLabelInterpreter','tex');

subplot(122)
boxplot(L2xhatstarMLsEnd,'Labels',Labels2)
bp=gca; bp.FontSize=20;
bp.XAxis.FontWeight='bold'; bp.XAxis.FontSize=bpXAxisFontSize;
bp.YAxis.FontWeight='bold'; bp.YAxis.FontSize=23;
ylim([0.00008,0.6])
set(gca,'YScale','log')
ylabel('$L_2(\hat{\textbf{x}}^*_{\mathbf{ML}})$','Interpreter','latex','FontSize',28);
set(findobj(gca,'type','line'),'LineWidth',2)
title('(b)','FontSize',25,'FontWeight','bold')
set(gca,'Position',[0.575 0.155 0.42 0.765])
set(gca,'yGrid','on','GridLineStyle','--')
yticks([10.^[-4:-1] 0.5])
yticklabels({'10^{-4}','10^{-3}','10^{-2}','10^{-1}','0.5'})
set(gca,'TickLabelInterpreter','tex');
set(findobj(gcf,'type','axes'),'FontWeight','Bold','LineWidth',2);
set(gcf,'position',[0 386 1920 510]) 

htmlGray=[128 128 128]/255;
htmlGreen=[0.4660 0.6740 0.1880];

Jump=4;
for i=1:9
    JJ(i)= mean(log(meanInterpolatedShxhatstarMLsminusShxstarML(InitialBudget:Budget,i)));
end
[~,indicesJJ]=sort(JJ);
for i=1:9
    Shift=mod(find(indicesJJ==i),Jump);
    II{i}=unique([InitialBudget (InitialBudget+Shift):Jump:Budget Budget]);        
end

figure,clf
tiledlayout(1,2,'Padding','none','TileSpacing','none');
nexttile
FontSize1=24;
linewidth=4;
MarkerSize1=15;
plot(1:Budget,meanInterpolatedShxhatstarMLsminusShxstarML(1:Budget,1),'ko-','linewidth',linewidth,'MarkerSize',MarkerSize1,'MarkerIndices',II{1}),hold on
plot(1:Budget,meanInterpolatedShxhatstarMLsminusShxstarML(1:Budget,2),'b:o','linewidth',linewidth,'MarkerSize',MarkerSize1,'MarkerFaceColor','b','MarkerIndices',II{2})
plot(1:Budget,meanInterpolatedShxhatstarMLsminusShxstarML(1:Budget,3),'k^-','linewidth',linewidth,'MarkerSize',MarkerSize1,'MarkerIndices',II{3})
plot(1:Budget,meanInterpolatedShxhatstarMLsminusShxstarML(1:Budget,4),'--v','linewidth',linewidth,'color',htmlGray,'MarkerSize',MarkerSize1,'MarkerIndices',II{4})
plot(1:Budget,meanInterpolatedShxhatstarMLsminusShxstarML(1:Budget,5),':s','linewidth',linewidth,'color',htmlGreen,'MarkerFaceColor',htmlGreen,'MarkerSize',MarkerSize1,'MarkerIndices',II{5})
plot(1:Budget,meanInterpolatedShxhatstarMLsminusShxstarML(1:Budget,6),'b-x','linewidth',linewidth,'MarkerSize',MarkerSize1+10,'MarkerIndices',II{6})
plot(1:Budget,meanInterpolatedShxhatstarMLsminusShxstarML(1:Budget,7),':rs','linewidth',linewidth,'MarkerSize',MarkerSize1,'MarkerIndices',II{7})
plot(1:Budget,meanInterpolatedShxhatstarMLsminusShxstarML(1:Budget,8),'--h','linewidth',linewidth,'color',[0.00,0.45,0.74],'MarkerFaceColor',[0.00,0.45,0.74],'MarkerSize',MarkerSize1,'MarkerIndices',II{8})
plot(1:Budget,meanInterpolatedShxhatstarMLsminusShxstarML(1:Budget,9),':d','linewidth',linewidth,'MarkerSize',MarkerSize1,'MarkerIndices',II{9})
xlabel('Computational cost');
set(gca,'YScale','log','FontSize',FontSize1,'FontWeight','bold','LineWidth',3);
ylabel('Average $S_h(\hat{\textbf{x}}^*_{\mathbf{ML}})-$0.093939','Interpreter','latex','FontSize',32);
leg=legend(Labels,'NumColumns',3,'Location','northeast');
leg.ItemTokenSize=[74,50];
ax = gca; ax.YMinorTick = 'off';
yticks([0.3 10.^[0:2]])
ylim([0.25 170])
yticklabels({'0.3 ','10^0 ','10^1 ','10^2 '})
xticks(InitialBudget:2:Budget)
xlim([InitialBudget,Budget])
set(gca,'TickLabelInterpreter','tex');
title('(a)','FontWeight','bold')

for i=1:9
    JJ(i)= mean(log(meanInterpolatedL2xhatstarMLs(InitialBudget:Budget,i)));
end
[~,indicesJJ]=sort(JJ);
for i=1:9
    Shift=mod(find(indicesJJ==i),Jump);
    II{i}=unique([InitialBudget (InitialBudget+Shift):Jump:Budget Budget]);
end

nexttile
plot(1:Budget,meanInterpolatedL2xhatstarMLs(1:Budget,1),'ko-','linewidth',linewidth,'MarkerSize',MarkerSize1,'MarkerIndices',II{1}),hold on
plot(1:Budget,meanInterpolatedL2xhatstarMLs(1:Budget,2),'b:o','linewidth',linewidth,'MarkerSize',MarkerSize1,'MarkerFaceColor','b','MarkerIndices',II{2})
plot(1:Budget,meanInterpolatedL2xhatstarMLs(1:Budget,3),'k^-','linewidth',linewidth,'MarkerSize',MarkerSize1,'MarkerIndices',II{3})
plot(1:Budget,meanInterpolatedL2xhatstarMLs(1:Budget,4),'--v','linewidth',linewidth,'color',htmlGray,'MarkerSize',MarkerSize1,'MarkerIndices',II{4})
plot(1:Budget,meanInterpolatedL2xhatstarMLs(1:Budget,5),':s','linewidth',linewidth,'color',htmlGreen,'MarkerFaceColor',htmlGreen,'MarkerSize',MarkerSize1,'MarkerIndices',II{5})
plot(1:Budget,meanInterpolatedL2xhatstarMLs(1:Budget,6),'b-x','linewidth',linewidth,'MarkerSize',MarkerSize1+10,'MarkerIndices',II{6})
plot(1:Budget,meanInterpolatedL2xhatstarMLs(1:Budget,7),':rs','linewidth',linewidth,'MarkerSize',MarkerSize1,'MarkerIndices',II{7})
plot(1:Budget,meanInterpolatedL2xhatstarMLs(1:Budget,8),'--h','linewidth',linewidth,'color',[0.00,0.45,0.74],'MarkerFaceColor',[0.00,0.45,0.74],'MarkerSize',MarkerSize1,'MarkerIndices',II{8})
plot(1:Budget,meanInterpolatedL2xhatstarMLs(1:Budget,9),':d','linewidth',linewidth,'MarkerSize',MarkerSize1,'MarkerIndices',II{9})
xlabel('Computational cost');
set(gca,'FontWeight','bold','FontSize',FontSize1);
ylabel('Average $L_2(\hat{\textbf{x}}^*_{\mathbf{ML}})$','Interpreter','latex','FontSize',32);
leg=legend(Labels,'NumColumns',3,'Location','northeast');
leg.ItemTokenSize=[74,50];
set(findobj(gcf,'type','axes'),'FontWeight','Bold','LineWidth',3);
set(gca,'YScale','log')
ax = gca; ax.YMinorTick = 'off';
yticks(0.01*2.^[1:5])
ylim([0.019 0.52])
yticklabels({'0.02 ','0.04 ','0.08 ','0.16 ','0.32 '})
xticks(InitialBudget:2:Budget)
xlim([InitialBudget,Budget])
title('(b)','FontWeight','bold')
set(gcf,'Position',[0 100 1920 615])

figure,clf
Labels2Method={'MBC-AGP','BC-AGP','BC-GP'};
boxplot(phiEnd(:,[1 2 7]),'Labels',Labels2Method,'OutlierSize',10,'Widths',0.8*[1 1 1])
set(findobj(gca,'type','line'),'LineWidth',2)
set(findobj(gcf,'type','axes'),'FontSize',27,'FontWeight','Bold','LineWidth',3);
ylabel('$ \hat \varphi$','Interpreter','latex','FontSize',50,'Rotation',0,'HorizontalAlignment','right','VerticalAlignment','baseline')
set(gca,'Position',[0.15 0.15 0.83 0.83])
yticks([-0.1:0.1:0.6])
set(gca,'yGrid','on','GridLineStyle','--')
ylim([-0.17 0.64])
set(gcf,'Position',[109 159 900 372])
medians=median(phiEnd(:,[1 2 7]));
FontSize2=20;
text(1,1.11*medians(1),['Median=' num2str(medians(1),2)],'HorizontalAlignment','center','FontSize',FontSize2,'FontWeight','Bold')
text(2,1.1*medians(2),['Median=' num2str(medians(2),2)],'HorizontalAlignment','center','FontSize',FontSize2,'FontWeight','Bold')
text(3,1.14*medians(3),['Median=' num2str(medians(3),2)],'HorizontalAlignment','center','FontSize',FontSize2,'FontWeight','Bold')
xlim([0.45 3.55])

figure;clf
idxTrial=88
tiledlayout(2,10,'Padding','none','TileSpacing','none');
pd1=1;
pd2=2;
FontSize3=18;
linewidth=1;
Markersize1=15;
for idxMethod=[1:6]
    if idxMethod==6
        nexttile([1 1])
        axis off

        nexttile([1 2])
    else
        
        nexttile([1 2])
    end
    
    Table=RecordTable{idxTrial,idxMethod};
    Level=Table.Level;
    AllDesignPoints=Table.D;
    FinalDh=AllDesignPoints(Level==2,:);
    FinalDl=AllDesignPoints(Level==1,:);

    xhatstarMLsEnd=Table.xhatstarMLs(end,:);
    
    InitialDh=FinalDh(1:nh,:);
    InitialDl=FinalDl(1:nl,:);
    FollowDh=FinalDh(nh+1:end,:);
    FollowDl=FinalDl(nl+1:end,:);
    
    plot(InitialDh(:,pd1),InitialDh(:,pd2),'bs','linewidth',linewidth,'markersize',Markersize1)
    hold on
    plot(InitialDl(:,pd1),InitialDl(:,pd2),'bx','linewidth',linewidth,'markersize',Markersize1)
    
    plot(FollowDh(:,pd1),FollowDh(:,pd2),'ko','linewidth',linewidth,'markersize',Markersize1)

    plot(FollowDl(:,pd1),FollowDl(:,pd2),'k+','linewidth',linewidth,'markersize',Markersize1)

    plot(xstarML(:,pd1),xstarML(:,pd2),'kp','MarkerSize',25)    
    plot(xhatstarMLsEnd(:,pd1),xhatstarMLsEnd(:,pd2),'k^','MarkerSize',25)
    
    xlabel('x_1','FontSize',FontSize3)
    ylabel('x_2','FontSize',FontSize3,'Rotation',0,'HorizontalAlignment','right')
    
    xticks([0:0.2:1])
    yticks([0:0.2:1])
    
    lim0=0.03;
    xlim([-lim0 1+lim0])
    ylim([-lim0 1+lim0])
    
    title(Labels(idxMethod),'FontWeight','Bold')
end

for idxMethod=7:9
    Table=RecordTable{idxTrial,idxMethod};
    FinalDh=Table.D;
    xhatstarMLsEnd=Table.xhatstarMLs(end,:);
    
    InitialDh=FinalDh(1:nh0,:);
    FollowDh=FinalDh(nh0+1:end,:);
    
    nexttile([1 2])
    
    plot(InitialDh(:,pd1),InitialDh(:,pd2),'bs','linewidth',linewidth,'markersize',Markersize1)
    hold on
    plot(FollowDh(:,pd1),FollowDh(:,pd2),'ko','linewidth',linewidth,'markersize',Markersize1)
 
    plot(xstarML(:,pd1),xstarML(:,pd2),'kp','MarkerSize',25)
    plot(xhatstarMLsEnd(:,pd1),xhatstarMLsEnd(:,pd2),'k^','MarkerSize',25)
    
    xlabel('x_1','FontSize',FontSize3)
    ylabel('x_2','FontSize',FontSize3,'Rotation',0,'HorizontalAlignment','right')
    
    xticks([0:0.2:1])
    yticks([0:0.2:1])
    
    lim0=0.03;
    xlim([-lim0 1+lim0])
    ylim([-lim0 1+lim0])
    
    title(Labels(idxMethod),'FontWeight','Bold')
end
set(findobj(gcf,'type','axes'),'FontSize',FontSize3,'FontWeight','Bold','LineWidth',1);
set(gcf,'Position',[0 0 1600 700])