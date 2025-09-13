%Example 2
%Section 1: Sets input data for all model calibration methods.
%Dl is the initial LF design (an nl x d matrix) and Dh is the initial HF design (an nh x d matrix) for each bi-fidelity method in a trial;
%Yl is the initial LF output data (an nl x N matrix) and Yh is the initial HF output data (an nh x N matrix) for each bi-fidelity method in a trial; 
%Dh0 is the initial HF design (an nh0 x d matrix) for each single-fidelity method in a trial;
%Yh0 is the initial HF output data (an nh0 x N matrix) for each single-fidelity method in a trial;
%CostRatio is c_h/c_l; 
%Budget is the total budget for each bi-fidelity method in a trial; 
%Budget1 is the total budget for each single-fidelity method in a trial;
%w is the vector of field/physical data (a 1 x N vector).

%%Uncomment the commands below in this section if you want to generate a different set of 100 initial designs for each method and a different vector 
%of field data, rather than use those stored in the .mat file loaded by this script, to perform another set of 100 trials.
%{
clear all;clc,format compact 
d=3;
Case=2;
nl=16;
nh=8;
nh0=12;
CostRatio=4;
InitialBudget=nl*1+nh*CostRatio;
InitialBudget0=nh0*CostRatio;
if(InitialBudget~=InitialBudget0)
    return
end
Budget=InitialBudget0+20;
Budget1=Budget;
xstar=[0.8535,0.1333,0.4258];
yhxstar=Simulator(xstar,2,Case);
stddev=(var(yhxstar)*0.0001)^0.5;
w=yhxstar+normrnd(0,stddev,size(yhxstar));
Shxstar=sum([yhxstar-w].^2);
[xstarML,ShxstarML]=Example2FindTrueMLE(w,Shxstar,xstar);

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
    clear Yh0
    parfor jd=1:nh0
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
    SingleDataInput(id).Budget=Budget1;     SingleDataInput(id).Case=Case;
       
end
%}
%save Example2InputData.mat
%clear all
%load Example2InputData.mat xstar yhxstar w Shxstar xstarML ShxstarML MultiDataInput SingleDataInput timeh timel
%save Example2InputData.mat
clear all;clc,format compact 
d=3;
Case=2;
nl=16;
nh=8;
nh0=12;
CostRatio=4;
InitialBudget=nl*1+nh*CostRatio;
InitialBudget0=nh0*CostRatio;
if(InitialBudget~=InitialBudget0)
    return
end
Budget=InitialBudget0+20;
Budget1=Budget;
load Example2InputData.mat xstar yhxstar w Shxstar xstarML ShxstarML MultiDataInput SingleDataInput timeh timel
%%
%Section 2: Runs all model calibration methods.
Z_BC=1;    Z_ID=0;   Z_SR=2;
ZMLFSSE=1; ZLFSSE=0; AccuracyLevel=1; t=0.99; t2=0.998;
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

    [T_SVD_AGP_2{id,1},~,RunTime_SVD_AGP_2(id)]=CalibrationSVDAGP(MultiDataInput(id),AccuracyLevel,t2); 'SVD-AGP, t=0.998'
    [T_SVD_2{id,1},~,RunTime_SVD_2(id)]=CalibrationSVD(SingleDataInput(id),AccuracyLevel,t2); 'SVD, t=0.998'
    save Example2.mat
end
%%
%Section 3: Constructs figures and a table to illustrate results.
clear all;clc,format compact 
load Example2.mat
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
disp('Table I.1')
Table2=table(Labels,AverageSh,ttest_pval_Sh,AverageL2,ttest_pval_L2)

htmlGray=[128 128 128]/255;
htmlGreen=[0.4660 0.6740 0.1880];

Jump=5;
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
ylabel('Average $S_h(\hat{\textbf{x}}^*_{\mathbf{ML}})-$2551.9','Interpreter','latex','FontSize',32);
leg=legend(Labels,'NumColumns',3,'Location','northeast');
leg.ItemTokenSize=[74,50];
yticks([10.^[3:6]])
ylim([1000 3000000])
xticks(InitialBudget:2:Budget)
xlim([InitialBudget-0.1,Budget+0.1])
set(gca,'TickLabelInterpreter','tex');
title('(a)','FontWeight','Bold')

for i=1:9
    JJ(i)= mean(meanInterpolatedL2xhatstarMLs(InitialBudget:Budget,i));
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
yticks([0.3:0.05:0.6])
ylim([0.2700 0.64000])
xticks(InitialBudget:2:Budget)
xlim([InitialBudget-0.1,Budget+0.1])
title('(b)','FontWeight','Bold')
set(gcf,'Position',[0 100 1920 615])

figure,clf
Labels2Method={'MBC-AGP','BC-AGP','BC-GP'};
boxplot(phiEnd(:,[1 2 7]),'Labels',Labels2Method,'OutlierSize',10,'Widths',0.8*[1 1 1])
set(findobj(gca,'type','line'),'LineWidth',2)
set(findobj(gcf,'type','axes'),'FontSize',27,'FontWeight','Bold','LineWidth',3);
ylabel('$ \hat \varphi$','Interpreter','latex','FontSize',50,'Rotation',0,'HorizontalAlignment','right','VerticalAlignment','baseline')
set(gca,'Position',[0.15 0.15 0.83 0.81])
yticks([-0.1:0.1:0.9])
set(gca,'yGrid','on','GridLineStyle','--')
ylim([0.24 0.7])
set(gcf,'Position',[109 159 900 372])
medians=median(phiEnd(:,[1 2 7]));
FontSize2=19;
text(1,1.045*medians(1),['Median=' num2str(medians(1),2)],'HorizontalAlignment','center','FontSize',FontSize2,'FontWeight','Bold')
text(2,0.96*medians(2),['Median=' num2str(medians(2),2)],'HorizontalAlignment','center','FontSize',FontSize2,'FontWeight','Bold')
text(3,1.05*medians(3),['Median=' num2str(medians(3),2)],'HorizontalAlignment','center','FontSize',FontSize2,'FontWeight','Bold')
xlim([0.45 3.55])

idxTrial=82;
pairindices=[1 2; 1 3; 2 3];
linewidth=1.5;
Markersize1=15;
for idxMethod=[1 3 6]
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
    
    figure,clf
    tiledlayout(1,3,'Padding','none','TileSpacing','none');
    
    for kd=1:3
        pd1=pairindices(kd,1);
        pd2=pairindices(kd,2);
        nexttile

        plot(InitialDh(:,pd1),InitialDh(:,pd2),'bs','linewidth',linewidth,'markersize',Markersize1)
        hold on
        plot(InitialDl(:,pd1),InitialDl(:,pd2),'bx','linewidth',linewidth,'markersize',Markersize1)
        
        plot(FollowDh(:,pd1),FollowDh(:,pd2),'ko','linewidth',linewidth,'markersize',Markersize1)
 
        plot(FollowDl(:,pd1),FollowDl(:,pd2),'k+','linewidth',linewidth,'markersize',Markersize1)
        
        plot(xstarML(:,pd1),xstarML(:,pd2),'kp','MarkerSize',25,'linewidth',linewidth)
        plot(xhatstarMLsEnd(:,pd1),xhatstarMLsEnd(:,pd2),'k^','MarkerSize',25,'linewidth',linewidth)
        
        xlabel(['x_' num2str(pd1)])
        ylabel(['x_' num2str(pd2)],'Rotation',0,'HorizontalAlignment','right')

        xticks([0:0.2:1])
        yticks([0:0.2:1])

        lim0=0.02;
        xlim([-lim0 1+lim0])
        ylim([-lim0 1+lim0])
                         
    end
    
    set(findobj(gcf,'type','axes'),'FontSize',17,'FontWeight','Bold','LineWidth',1);
    if idxMethod==1
        sgtitle('(a) MBC-AGP','fontsize',25,'FontWeight','Bold')
        set(gcf,'Position',[0 650 1600 350])        
    elseif idxMethod==3
        sgtitle('(b) MID-AGP','fontsize',25,'FontWeight','Bold')
        set(gcf,'Position',[0 450 1600 350])        
    elseif idxMethod==6
        sgtitle('(c) SVD-AGP','fontsize',25,'FontWeight','Bold')
        set(gcf,'Position',[0 250 1600 350])
    end
    
end

for idxMethod=9
    Table=RecordTable{idxTrial,idxMethod};
    Level=Table.Level;
    AllDesignPoints=Table.D;
    FinalDh=AllDesignPoints(Level==2,:);

    xhatstarMLsEnd=Table.xhatstarMLs(end,:);
    
    InitialDh=FinalDh(1:nh0,:);
    FollowDh=FinalDh(nh0+1:end,:);
    
    figure,clf
    tiledlayout(1,3,'Padding','none','TileSpacing','none');
    
    for kd=1:3
        pd1=pairindices(kd,1);
        pd2=pairindices(kd,2);
        nexttile
        Markersize1=15;
        linewidth=1.5;
        plot(InitialDh(:,pd1),InitialDh(:,pd2),'bs','linewidth',linewidth,'markersize',Markersize1)
        hold on       
        plot(FollowDh(:,pd1),FollowDh(:,pd2),'ko','linewidth',linewidth,'markersize',Markersize1)
        
        plot(xstarML(:,pd1),xstarML(:,pd2),'kp','MarkerSize',25,'linewidth',linewidth)
        plot(xhatstarMLsEnd(:,pd1),xhatstarMLsEnd(:,pd2),'k^','MarkerSize',25,'linewidth',linewidth)
        
        xlabel(['x_' num2str(pd1)])
        ylabel(['x_' num2str(pd2)],'Rotation',0,'HorizontalAlignment','right')

        xticks([0:0.2:1])
        yticks([0:0.2:1])

        lim0=0.02;
        xlim([-lim0 1+lim0])
        ylim([-lim0 1+lim0])
            
    end
    
    set(findobj(gcf,'type','axes'),'FontSize',17,'FontWeight','Bold','LineWidth',1);

    sgtitle('(d) SVD','fontsize',25,'FontWeight','Bold')
    set(gcf,'Position',[0 0 1600 350])        
    
end