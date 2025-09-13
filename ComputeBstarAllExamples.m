%%
clear all;clc,format compact;format long g 

for ii=1:6
    clearvars -except ii
    if(ii==1)
        load Example1GridData26.mat  AllYh AllYl w GridPoints
        nlevels=26;
    elseif(ii==2)
        load Example1GridData26.mat  AllYh AllYl w GridPoints
        nlevels=26;                    
    elseif(ii==3)
        filename = 'Example2GridData.mat';  
        if isfile(filename)==1
        load Example2.mat w          
        load Example2GridData.mat  AllYh AllYl GridPoints 
        nlevels=41; 
        else
        disp('You need to run ApproximateSh_SVD.m to get the Example2GridData.mat data file.');
        break
        end        
    elseif(ii==4)
        load Example2Size2.mat w        
        load Example2GridData.mat  AllYh AllYl GridPoints 
        nlevels=41;                           
    elseif(ii==5)
        load Example3.mat w
        nlevels=51;
        [X1,X2]=meshgrid(linspace(0,1,nlevels));
        GridPoints=[X1(:) X2(:)];
        NGP=size(GridPoints,1);
        for id=1:NGP
            AllYh(id,:)=Simulator(GridPoints(id,:),2,3);
            AllYl(id,:)=Simulator(GridPoints(id,:),1,3);   
        end
    elseif(ii==6)    
        load Example3Size2.mat w
        nlevels=51;
        [X1,X2]=meshgrid(linspace(0,1,nlevels));
        GridPoints=[X1(:) X2(:)];
        NGP=size(GridPoints,1);
        for id=1:NGP
            AllYh(id,:)=Simulator(GridPoints(id,:),2,3);
            AllYl(id,:)=Simulator(GridPoints(id,:),1,3);   
        end
    end

AllSh=sum((AllYh-w).^2,2);
AllSl=sum((AllYl-w).^2,2);

nlevelsm1=nlevels-1;
nlevelsm2=nlevels-2;
cellno=1;
if(ii<=4)
for a=0:(1/nlevelsm1):nlevelsm2*(1/nlevelsm1)
    for b=0:(1/nlevelsm1):nlevelsm2*(1/nlevelsm1)
        for c=0:(1/nlevelsm1):nlevelsm2*(1/nlevelsm1)
            Lowervertex=[a,b,c];
            V=(fullfact([2 2 2])-1)*(1/nlevelsm1);
            Allvertices{cellno}=repmat(Lowervertex,8,1)+V;
            cellno=cellno+1;
        end
    end
end
else
for a=0:(1/nlevelsm1):nlevelsm2*(1/nlevelsm1)
    for b=0:(1/nlevelsm1):nlevelsm2*(1/nlevelsm1)
        Lowervertex=[a,b];
        V=(fullfact([2 2])-1)*(1/nlevelsm1);
        Allvertices{cellno}=repmat(Lowervertex,4,1)+V;
        cellno=cellno+1;
    end
end    
end

if(ii<=4)
for k=1:nlevelsm1^3
    [mindists,minidx]=min(pdist2(Allvertices{k},GridPoints),[],2);
    Cellverticesindices(k,:)=minidx;
    CellPoints=GridPoints(Cellverticesindices(k,:),:);
    vol=prod(range(CellPoints,1));
    if((size(CellPoints,1)~=8) || (abs(vol-(1/nlevelsm1)^3)>10^-12) || (min(pdist(CellPoints))<10^-12) || (max(mindists)>10^-12))
        input('error')
        return
    end
end
else
for k=1:nlevelsm1^2
    [mindists,minidx]=min(pdist2(Allvertices{k},GridPoints),[],2);
    Cellverticesindices(k,:)=minidx;
    CellPoints=GridPoints(Cellverticesindices(k,:),:);
    vol=prod(range(CellPoints,1));
    if((size(CellPoints,1)~=4) || (abs(vol-(1/nlevelsm1)^2)>10^-12) || (min(pdist(CellPoints))<10^-12) || (max(mindists)>10^-12))
        input('error')
        return
    end
end    
end

AllSh_Vertices=AllSh(Cellverticesindices);
AllSl_Vertices=AllSl(Cellverticesindices);

NoTrials=100;

    if(ii==1)
        load Example1.mat MultiDataInput T_MBC_AGP T_SR_AGP
    elseif(ii==2)
        load Example1Size2.mat MultiDataInput T_MBC_AGP T_SR_AGP
    elseif(ii==3)
        load Example2.mat MultiDataInput T_MBC_AGP T_SR_AGP
    elseif(ii==4)
        load Example2Size2.mat MultiDataInput T_MBC_AGP T_SR_AGP
    elseif(ii==5)
        load Example3.mat MultiDataInput T_MBC_AGP T_SR_AGP
    elseif(ii==6)    
        load Example3Size2.mat MultiDataInput T_MBC_AGP T_SR_AGP
    end
    
corr_ZhZlPlus_MBC_AGP=zeros(NoTrials,1); 
phis_MBC_AGP=zeros(NoTrials,1); 
rhos_MBC_AGP=zeros(NoTrials,1); 
Bstar_MBC_AGP=zeros(NoTrials,1);
ID_BC_or_SR=1;

for Trial=1:NoTrials
    Dl=MultiDataInput(Trial).Dl; Dh=MultiDataInput(Trial).Dh;
    Yl=MultiDataInput(Trial).Yl; Yh=MultiDataInput(Trial).Yh;
    [AllYlModified,ahati_bhati]=regress_aibi(Dl,Dh,Yl,Yh,AllYl);
    AllSlPlus=sum((AllYlModified-w).^2,2);
    
    AllSlPlus_Vertices=AllSlPlus(Cellverticesindices);
    
    phis_MBC_AGP(Trial,1)=T_MBC_AGP{Trial,1}.phis(end);
    rhos_MBC_AGP(Trial,1)=T_MBC_AGP{Trial,1}.rhos(end);
    phi=phis_MBC_AGP(Trial,1);
    
    Zh_Vertices=TransformData(AllSh_Vertices,phi,ID_BC_or_SR);
    
    Mean_Zh=mean(Zh_Vertices,'all');
    Zh2_Vertices=Zh_Vertices.^2;
    Mean_Zh2=mean(Zh2_Vertices,'all');
    Var_Zh=Mean_Zh2-(Mean_Zh)^2;

    ZlPlus_Vertices=TransformData(AllSlPlus_Vertices,phi,ID_BC_or_SR);
    
    ZhZlPlus_Vertices=Zh_Vertices.*ZlPlus_Vertices;
    Mean_ZhZlPlus=mean(ZhZlPlus_Vertices,'all');    
    Mean_ZlPlus=mean(ZlPlus_Vertices,'all');    
    ZlPlus2_Vertices=ZlPlus_Vertices.^2;
    Mean_ZlPlus2=mean(ZlPlus2_Vertices,'all');
    Var_ZlPlus=Mean_ZlPlus2-(Mean_ZlPlus)^2;
    
    corr_ZhZlPlus_MBC_AGP(Trial,1)=(Mean_ZhZlPlus-Mean_Zh*Mean_ZlPlus)/(Var_Zh^0.5*Var_ZlPlus^0.5);
    Bstar_MBC_AGP(Trial,1)=corr_ZhZlPlus_MBC_AGP(Trial,1)*sqrt(Var_Zh)/sqrt(Var_ZlPlus);
end

rhos_SR_AGP=zeros(NoTrials,1);
ID_BC_or_SR=2;

Zh_Vertices=TransformData(AllSh_Vertices,0,ID_BC_or_SR);
    
Mean_Zh=mean(Zh_Vertices,'all');
Zh2_Vertices=Zh_Vertices.^2;
Mean_Zh2=mean(Zh2_Vertices,'all');
Var_Zh=Mean_Zh2-(Mean_Zh)^2;

Zl_Vertices=TransformData(AllSl_Vertices,0,ID_BC_or_SR);
    
ZhZl_Vertices=Zh_Vertices.*Zl_Vertices;
Mean_ZhZl=mean(ZhZl_Vertices,'all');    
Mean_Zl=mean(Zl_Vertices,'all');    
Zl2_Vertices=Zl_Vertices.^2;
Mean_Zl2=mean(Zl2_Vertices,'all');
Var_Zl=Mean_Zl2-(Mean_Zl)^2;
    
corr_ZhZl_SR_AGP=(Mean_ZhZl-Mean_Zh*Mean_Zl)/(Var_Zh^0.5*Var_Zl^0.5);
Bstar_SR_AGP=corr_ZhZl_SR_AGP*sqrt(Var_Zh)/sqrt(Var_Zl);

for Trial=1:NoTrials
    rhos_SR_AGP(Trial,1)=T_SR_AGP{Trial,1}.rhos(end);    
end
    if(ii==1)
        disp('Example 1, 100 trials for Figure 2:')
    elseif(ii==2)
        disp('Example 1, 100 trials with larger initial design size:')
    elseif(ii==3)        
        disp('Example 2, 100 trials for Figure 5:')                
    elseif(ii==4)
        disp('Example 2, 100 trials with larger initial design size:')    
    elseif(ii==5)
        disp('Example 3, 100 trials for Figure F.3:')  
    elseif(ii==6)    
        disp('Example 3, 100 trials with smaller initial design size:')
    end
disp('Minimum and maximum of Bstar for the MBC-AGP method over the 100 trials are:')
[min(Bstar_MBC_AGP) max(Bstar_MBC_AGP)]
disp('Minimum and maximum of rho for the MBC-AGP method over the 100 trials are:')
[min(rhos_MBC_AGP) max(rhos_MBC_AGP)]

disp('Bstar for the SR-AGP method is:')
[Bstar_SR_AGP]
disp('Minimum and maximum of rho for the SR-AGP method over the 100 trials are:')
[min(rhos_SR_AGP) max(rhos_SR_AGP)]

end

function [YlTestModified,ahati_bhati]=regress_aibi(Dl,Dh,Yl,Yh,YlTest)
N=size(Yl,2);
[~,idxinDl,idxinDh]=intersect(Dl,Dh,'rows','stable');
SameYl=Yl(idxinDl,:);
SameYh=Yh(idxinDh,:);
OnesVec=ones(numel(idxinDh),1);
ahati_bhati=zeros(2,N);
for kd=1:N
    ModelMatrix=[OnesVec,SameYl(:,kd)];
    lastwarn('');
    ahati_bhati(:,kd)=regress(SameYh(:,kd),ModelMatrix);
    [warnMsg,~]=lastwarn;
    if contains(warnMsg,'X is rank deficient to within machine precision.')
        return
    end
end
YlTestModified=ahati_bhati(1,:)+YlTest.*ahati_bhati(2,:);
end