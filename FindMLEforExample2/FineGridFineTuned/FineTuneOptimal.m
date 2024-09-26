
clear all
load Example2XMLEFineGridOkOptTop20.mat

lb = 0*ones(1,Dim);ub = 1*ones(1,Dim);
options=optimoptions('patternsearch','MaxIterations',10^6,'MeshTolerance',10^-4,'TolFun',10^-6,'MaxFunEvals',10^8);

SSHFun=@(x) sum((Simulator(x,2,Case)-PhysData).^2);
[XMLE2FineTuned,SSE_XMLE2FineTuned,exitflag] = patternsearch(SSHFun,XMLE,[],[],[],[],lb,ub,[],options)  ;

save Example2XMLEFineGridOkOptTop20FineTuned.mat
clearvars -except SSE_XMLE2FineTuned XMLE2FineTuned
save Example2TrueMLEFineTuned.mat