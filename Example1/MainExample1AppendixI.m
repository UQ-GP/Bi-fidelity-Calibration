clear,clc
load Example1GridData.mat  AllYh AllYl PhysData Points
load Example1.mat MultiDataInput
AllSh=sum((AllYh-PhysData).^2,2);
AllSl=sum((AllYl-PhysData).^2,2);

phi=0.38;
% Points 11^3-by 3
% AllYh 11^3-by 22
% AllYl 11^3-by 22
nlevels=11;
nlevelsm1=nlevels-1;
nlevelsm2=nlevels-2;
cellno=1;
for a=0:(1/nlevelsm1):nlevelsm2*(1/nlevelsm1)
    for b=0:(1/nlevelsm1):nlevelsm2*(1/nlevelsm1)
        for c=0:(1/nlevelsm1):nlevelsm2*(1/nlevelsm1)
            Lowervertex=[a,b,c];
            V=(fullfact([2 2  2])-1)*(1/nlevelsm1);
            Allvertices{cellno}=repmat(Lowervertex,8,1)+V ;
            cellno=cellno+1;
        end
    end
end
% grid

aa=0;
for k=1:nlevelsm1^3
    [aaa,minidx]=min(pdist2(Allvertices{k},Points),[],2);
    Cellverticesindices(k,:)=minidx;
    kkk=Points(Cellverticesindices(k,:),:);
    vol=prod(range(kkk,1));
    if((size(kkk,1)~=8)||(abs(vol-(1/nlevelsm1)^3)>10^-12)||(min(pdist(kkk))<10^-12)||max(aaa)>10^-12)
        k
        size(kkk,1)
        vol
        input('xxx')
    end
end

% L2 (Sh, Sl)
AllSh_Vertice=AllSh(Cellverticesindices);
AllSl_Vertice=AllSl(Cellverticesindices);

% corr (Sh, Sl),  L2(Sh,Sl)
ShSl_Vertice=AllSh_Vertice.*AllSl_Vertice;
Mean_ShSl=mean(ShSl_Vertice,'all');

Sh_Vertice=AllSh_Vertice;
Mean_Sh=mean(Sh_Vertice,'all');

Sl_Vertice=AllSl_Vertice;
Mean_Sl=mean(Sl_Vertice,'all');

Sh2_Vertice=AllSh_Vertice.^2;
Mean_Sh2=mean(Sh2_Vertice,'all');
Var_Sh=Mean_Sh2-(Mean_Sh)^2;

Sl2_Vertice=AllSl_Vertice.^2;
Mean_Sl2=mean(Sl2_Vertice,'all');
Var_Sl=Mean_Sl2-(Mean_Sl)^2;

corr_ShSl=(Mean_ShSl-Mean_Sh*Mean_Sl)/ ( Var_Sh^0.5*Var_Sl^0.5  );
corr_ShSl_codes=corr(AllSh,AllSl);

L2_ShSl_Vertice=(AllSh_Vertice-AllSl_Vertice).^2;
L2_ShSl=mean(L2_ShSl_Vertice,'all');

ZNBC=1;
Zh_Vertice=TransformData(AllSh_Vertice,phi,ZNBC);
Zl_Vertice=TransformData(AllSl_Vertice,phi,ZNBC);

ZhZl_Vertice=Zh_Vertice.*Zl_Vertice;
Mean_ZhZl=mean(ZhZl_Vertice,'all');

Mean_Zh=mean(Zh_Vertice,'all');
Mean_Zl=mean(Zl_Vertice,'all');

Zh2_Vertice=Zh_Vertice.^2;
Mean_Zh2=mean(Zh2_Vertice,'all');
Var_Zh=Mean_Zh2-(Mean_Zh)^2;

Zl2_Vertice=Zl_Vertice.^2;
Mean_Zl2=mean(Zl2_Vertice,'all');
Var_Zl=Mean_Zl2-(Mean_Zl)^2;

Zh=TransformData(AllSh,phi,ZNBC);
Zl=TransformData(AllSl,phi,ZNBC);

% corr(Zh, Zl+)      , L2 (Sh, Sl)
NoTrials=100;
L2_ShSlPlus=zeros(NoTrials,1);
L2_ShSlPlus_Way2=zeros(NoTrials,1);
corr_ShSlPlus=zeros(NoTrials,1);


for Trial=1:NoTrials
    Dl=MultiDataInput(Trial).Dl;    Dh=MultiDataInput(Trial).Dh;
    Yl=MultiDataInput(Trial).Yl;    Yh=MultiDataInput(Trial).Yh;
    [AllYlModified,ai_bi]=regress_aibi(Dl,Dh,Yl,Yh,AllYl);
    SlPlus=sum((AllYlModified-PhysData).^2,2);
    
    AllSlPlus_Vertice=SlPlus(Cellverticesindices);
    
    L2_ShSlPlus_Vertice=(AllSh_Vertice-AllSlPlus_Vertice).^2;
    L2_ShSlPlus(Trial,1)=mean(L2_ShSlPlus_Vertice,'all');
    
    ShSlPlus_Vertice=AllSh_Vertice.*AllSlPlus_Vertice;
    Mean_ShSlPlus=mean(ShSlPlus_Vertice,'all');
    
    SlPlus_Vertice=AllSlPlus_Vertice;
    Mean_SlPlus=mean(SlPlus_Vertice,'all');
    
    SlPlus2_Vertice=AllSlPlus_Vertice.^2;
    Mean_SlPlus2=mean(SlPlus2_Vertice,'all');
    Var_SlPlus=Mean_SlPlus2-(Mean_SlPlus)^2;
    
    corr_ShSlPlus(Trial,1)=(Mean_ShSlPlus-Mean_Sh*Mean_SlPlus)/ ( Var_Sh^0.5*Var_SlPlus^0.5  );
    
    
    ZlPlus_Vertice=TransformData(AllSlPlus_Vertice,phi,ZNBC);
    
    Zh_Vertice=TransformData(AllSh_Vertice,phi,ZNBC);
    
    
    ZhZlPlus_Vertice=Zh_Vertice.*ZlPlus_Vertice;
    Mean_ZhZlPlus=mean(ZhZlPlus_Vertice,'all');
    
    Mean_Zh=mean(Zh_Vertice,'all');
    Mean_ZlPlus=mean(ZlPlus_Vertice,'all');
    
    Zh2_Vertice=Zh_Vertice.^2;
    Mean_Zh2=mean(Zh2_Vertice,'all');
    Var_Zh=Mean_Zh2-(Mean_Zh)^2;
    
    ZlPlus2_Vertice=ZlPlus_Vertice.^2;
    Mean_ZlPlus2=mean(ZlPlus2_Vertice,'all');
    Var_ZlPlus=Mean_ZlPlus2-(Mean_ZlPlus)^2;
    
    corr_ZhZlPlus(Trial,1)=(Mean_ZhZlPlus-Mean_Zh*Mean_ZlPlus)/ ( Var_Zh^0.5*Var_ZlPlus^0.5  );
    
    NormalizedL2_ShSlPlus(Trial,1)=(L2_ShSlPlus(Trial,1)/Mean_Sh2)^0.5;
end

%%%%%%Conclusion in the Appendix I

%‖S_h (∙)-S_l (∙)‖_2/‖S_h (∙)‖_2
NormalizedL2_ShSl=(L2_ShSl/Mean_Sh2)^0.5%@@@A 

%The 0.5 and 0.9 quantiles of ‖S_h (∙)-S_l^+ (∙)‖_2/‖S_h (∙)‖_2
prctile(NormalizedL2_ShSlPlus,[ 50 90 ])%disp('0.5 0.9 quantile')
%= @@@ B  

%the median and 0.05 quantile of the correlation between g_φ (S_l^+ (X)) and g_φ (S_h (X)) over the 100 trials 
prctile(corr_ZhZlPlus,[ 50 5 ])
% = @@@C  

%the correlation between g_φ (S_h (X)) and g_φ (S_l (X))
corr_ZhZl=(Mean_ZhZl-Mean_Zh*Mean_Zl)/ ( Var_Zh^0.5*Var_Zl^0.5  )
% = @@@D 

function [YlTestModified,ai_bi,aicolumn_bicolumn]=regress_aibi(Dl,Dh,Yl,Yh,YlTest)
NTimePoints=size(Yl,2);
[~,idxinDl,idxinDh]=intersect(Dl,Dh,'rows','stable');
SameYl=Yl(idxinDl,:);
SameYh=Yh(idxinDh,:);
Ones2=ones(numel(idxinDh),1);
for kd=1:NTimePoints
    ai_bi(:,kd)=regress(SameYh(:,kd),[Ones2,SameYl(:,kd)]);
end
aicolumn_bicolumn=ai_bi';
aicolumn_bicolumn=aicolumn_bicolumn(:);
YlTestModified=ai_bi(1,:)+YlTest.*ai_bi(2,:);
end