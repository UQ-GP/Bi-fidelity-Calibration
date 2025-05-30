%Example 2: Larger initial design size
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
nl=20;
nh=10;
nh0=15;
CostRatio=4;
InitialBudget=nl*1+nh*CostRatio;
InitialBudget0=nh0*CostRatio;
if(InitialBudget~=InitialBudget0)
    return
end
Budget=InitialBudget0+20;
Budget1=Budget;
%xstar=[0.8535,0.1333,0.4258];
%yhxstar=Simulator(xstar,2,Case);
%stddev=(var(yhxstar)*0.0001)^0.5;
%w=yhxstar+normrnd(0,stddev,size(yhxstar));
%Shxstar=sum([yhxstar-w].^2);
%[xstarML,ShxstarML]=Example2FindTrueMLE(w,Shxstar,xstar);
%Uncomment the six rows above and comment away the one row below if you want to generate and use a different vector of field data than that used in the MainExample2.m script.
load Example2.mat xstar yhxstar w Shxstar xstarML ShxstarML

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
%save Example2Size2InputData.mat
%clear all
%load Example2Size2InputData.mat xstar yhxstar w Shxstar xstarML ShxstarML MultiDataInput SingleDataInput timeh timel
%save Example2Size2InputData.mat
clear all;clc,format compact 
d=3;
Case=2;
nl=20;
nh=10;
nh0=15;
CostRatio=4;
InitialBudget=nl*1+nh*CostRatio;
InitialBudget0=nh0*CostRatio;
if(InitialBudget~=InitialBudget0)
    return
end
Budget=InitialBudget0+20;
Budget1=Budget;
load Example2Size2InputData.mat xstar yhxstar w Shxstar xstarML ShxstarML MultiDataInput SingleDataInput timeh timel
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
    save Example2Size2.mat
end
%%
%Section 3: Constructs figures and a table to illustrate results.
clear all;clc,format compact 
load Example2Size2.mat
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
Table2=table(Labels,AverageSh,ttest_pval_Sh,AverageL2,ttest_pval_L2)

htmlGray=[128 128 128]/255;
htmlGreen=[0.4660 0.6740 0.1880];

figure,clf
tiledlayout(1,2,'Padding','none','TileSpacing','none');
nexttile
FontSize1=24;
linewidth=4;
MarkerSize1=15;
plot(1:Budget,meanInterpolatedShxhatstarMLsminusShxstarML(1:Budget,1),'ko-','linewidth',linewidth,'MarkerSize',MarkerSize1,'MarkerIndices',[InitialBudget:3:Budget Budget]),hold on
plot(1:Budget,meanInterpolatedShxhatstarMLsminusShxstarML(1:Budget,2),'b:o','linewidth',linewidth,'MarkerSize',MarkerSize1,'MarkerFaceColor','b','MarkerIndices',[InitialBudget (InitialBudget+2):2:Budget Budget])
plot(1:Budget,meanInterpolatedShxhatstarMLsminusShxstarML(1:Budget,3),'k^-','linewidth',linewidth,'MarkerSize',MarkerSize1,'MarkerIndices',[InitialBudget (InitialBudget):3:Budget Budget])
plot(1:Budget,meanInterpolatedShxhatstarMLsminusShxstarML(1:Budget,4),'--v','linewidth',linewidth,'color',htmlGray,'MarkerSize',MarkerSize1,'MarkerIndices',[(InitialBudget+1):2:Budget Budget])
plot(1:Budget,meanInterpolatedShxhatstarMLsminusShxstarML(1:Budget,5),':s','linewidth',linewidth,'color',htmlGreen,'MarkerFaceColor',htmlGreen,'MarkerSize',MarkerSize1,'MarkerIndices',[InitialBudget:4:Budget])
plot(1:Budget,meanInterpolatedShxhatstarMLsminusShxstarML(1:Budget,6),'b-x','linewidth',linewidth,'MarkerSize',MarkerSize1+10,'MarkerIndices',[InitialBudget (InitialBudget+1):3:Budget Budget])
plot(1:Budget,meanInterpolatedShxhatstarMLsminusShxstarML(1:Budget,7),':rs','linewidth',linewidth,'MarkerSize',MarkerSize1,'MarkerIndices',[InitialBudget (InitialBudget+2):3:Budget Budget])
plot(1:Budget,meanInterpolatedShxhatstarMLsminusShxstarML(1:Budget,8),'--h','linewidth',linewidth,'color',[0.00,0.45,0.74],'MarkerFaceColor',[0.00,0.45,0.74],'MarkerSize',MarkerSize1,'MarkerIndices',[InitialBudget (InitialBudget+1):3:Budget Budget])
plot(1:Budget,meanInterpolatedShxhatstarMLsminusShxstarML(1:Budget,9),':d','linewidth',linewidth,'MarkerSize',MarkerSize1,'MarkerIndices',[InitialBudget (InitialBudget+3):3:Budget Budget])
xlabel('Computational cost');
set(gca,'YScale','log','FontSize',FontSize1,'FontWeight','bold','LineWidth',3);
ylabel('Average $S_h(\hat{\textbf{x}}^*_{\mathbf{ML}})-$2551.9','Interpreter','latex','FontSize',32);
leg=legend(Labels,'NumColumns',3,'Location','northeast');
leg.ItemTokenSize=[74,50];
yticks([10.^[0:8]])
ylim([1000 3000000])
xticks(InitialBudget:2:Budget)
xlim([InitialBudget-0.1,Budget+0.1])
set(gca,'TickLabelInterpreter','tex');
title('(a)','FontWeight','bold')

nexttile
plot(1:Budget,meanInterpolatedL2xhatstarMLs(1:Budget,1),'ko-','linewidth',linewidth,'MarkerSize',MarkerSize1,'MarkerIndices',[InitialBudget:3:Budget Budget]),hold on
plot(1:Budget,meanInterpolatedL2xhatstarMLs(1:Budget,2),'b:o','linewidth',linewidth,'MarkerSize',MarkerSize1,'MarkerFaceColor','b','MarkerIndices',[InitialBudget (InitialBudget):2:Budget Budget])
plot(1:Budget,meanInterpolatedL2xhatstarMLs(1:Budget,3),'k^-','linewidth',linewidth,'MarkerSize',MarkerSize1,'MarkerIndices',[InitialBudget (InitialBudget):3:Budget Budget])
plot(1:Budget,meanInterpolatedL2xhatstarMLs(1:Budget,4),'--v','linewidth',linewidth,'color',htmlGray,'MarkerSize',MarkerSize1,'MarkerIndices',[(InitialBudget+1):2:Budget Budget])
plot(1:Budget,meanInterpolatedL2xhatstarMLs(1:Budget,5),':s','linewidth',linewidth,'color',htmlGreen,'MarkerFaceColor',htmlGreen,'MarkerSize',MarkerSize1,'MarkerIndices',[InitialBudget:4:Budget Budget])
plot(1:Budget,meanInterpolatedL2xhatstarMLs(1:Budget,6),'b-x','linewidth',linewidth,'MarkerSize',MarkerSize1+10,'MarkerIndices',[InitialBudget (InitialBudget+1):3:Budget Budget])
plot(1:Budget,meanInterpolatedL2xhatstarMLs(1:Budget,7),':rs','linewidth',linewidth,'MarkerSize',MarkerSize1,'MarkerIndices',[InitialBudget (InitialBudget+2):3:Budget Budget])
plot(1:Budget,meanInterpolatedL2xhatstarMLs(1:Budget,8),'--h','linewidth',linewidth,'color',[0.00,0.45,0.74],'MarkerFaceColor',[0.00,0.45,0.74],'MarkerSize',MarkerSize1,'MarkerIndices',[InitialBudget (InitialBudget+1):3:Budget Budget])
plot(1:Budget,meanInterpolatedL2xhatstarMLs(1:Budget,9),':d','linewidth',linewidth,'MarkerSize',MarkerSize1,'MarkerIndices',[InitialBudget (InitialBudget+3):3:Budget Budget])
xlabel('Computational cost');
set(gca,'FontWeight','bold','FontSize',FontSize1);
ylabel('Average $L_2(\hat{\textbf{x}}^*_{\mathbf{ML}})$','Interpreter','latex','FontSize',32);
leg=legend(Labels,'NumColumns',3,'Location','northeast');
leg.ItemTokenSize=[74,50];
set(findobj(gcf,'type','axes'),'FontWeight','Bold','LineWidth',3);
yticks([0.25:0.05:0.6])
ylim([0.2400 0.68000])
xticks(InitialBudget:2:Budget)
xlim([InitialBudget-0.1,Budget+0.1])
title('(b)','FontWeight','bold')
set(gcf,'Position',[0 100 1920 615])