%for SVD
clear all
load DataforAppendixJ_SVD.mat

sortedAllSh=sort(AllSh); rangeSh=sortedAllSh(end)-sortedAllSh(1); RangeApproxShs=range(PrctileApproxShs,2);
disp('min(range of S_{h,1}), max(range of S_{h,1}), and range of S_h')
[min(RangeApproxShs) max(RangeApproxShs) rangeSh]

for i=1:100
    ApproxxstarMLTryTrial=ApproxxstarMLTry(:,:,i);
    DistanceBetweenEachPairof20fminconSolutions=pdist(ApproxxstarMLTryTrial);
    MaxDistanceBetweenEachPairof20fminconSolutions(i)=max(DistanceBetweenEachPairof20fminconSolutions);
    RangeofApproximateShOver20fminconSolutions(i)=range(ApproxShBestValTry(:,i));
    AnyTies(i)=sum(DistanceBetweenEachPairof20fminconSolutions==0);
end
if(sum(AnyTies)>0)
    disp('There are repeated solutions returned by fmincon in optimizing S_{h,1} in some trials.')
end
disp('Minimum and mean over the 100 trials of the maximum L2 distance between any two of the 20 solutions returned by fmincon in optimizing S_{h,1}:')
[min(MaxDistanceBetweenEachPairof20fminconSolutions) mean(MaxDistanceBetweenEachPairof20fminconSolutions)]
disp('Maximum over the 100 trials of the range of the S_{h,1} values at the 20 solutions returned by fmincon in optimizing S_{h,1}:')
[max(RangeofApproximateShOver20fminconSolutions)]

disp('Sequence of mean(S_{h,1}(xhatstarML)-min{S_{h,1}})')
mean((store5-store3(:,end)))
disp('Sequence of mean(S_h(xhatstarML)-min{S_h})')
mean((store6-store2(:,end-1)))
% disp('Sequence of mean(||xhatstarML-ApproxxstarML||_2), where ApproxxstarML=minimizer of S_{h,1}')
% mean(store7)
% disp('Sequence of mean(||xhatstarML-xstarML||_2))')
% mean(store8)
disp('Sequence of mean of fraction of points in 41^3 grid with smaller S_{h,1} value than S_{h,1}(xhatstarML):')
mean(store9/41^3)
disp('Sequence of mean of fraction of points in 41^3 grid with smaller S_h value than S_h(xhatstarML):')
mean(store10/41^3)

disp('For Figure I.2, S_{h,1}(xhatstarML)-min{S_{h,1}}, and S_h(xhatstarML)-min{S_h} are:')
[store4(82,end)-store3(82,end) store4(82,end-1)-store2(82,end-1)]
%A=[vecnorm(store4(:,1:3)-store3(:,1:3),2,2),vecnorm(store4(:,1:3)-store2(:,1:3),2,2)];
% disp('For Figure I.2, distance between xhatstarML and ApproxxstarML (minimizer of S_{h,1}), and distance between xhatstarML and xstarML are:')
% A(82,:)

disp('Above are some results for SVD method.')
%%
%for SVDAGP
load DataforAppendixJ_SVDAGP.mat

sortedAllSh=sort(AllSh); rangeSh=sortedAllSh(end)-sortedAllSh(1); RangeApproxShs=range(PrctileApproxShs,2);
disp('min(range of S_{h,2}), max(range of S_{h,2}), and range of S_h')
[min(RangeApproxShs) max(RangeApproxShs) rangeSh]

for i=1:100
    ApproxxstarMLTryTrial=ApproxxstarMLTry(:,:,i);
    DistanceBetweenEachPairof20fminconSolutions=pdist(ApproxxstarMLTryTrial);
    MaxDistanceBetweenEachPairof20fminconSolutions(i)=max(DistanceBetweenEachPairof20fminconSolutions);
    RangeofApproximateShOver20fminconSolutions(i)=range(ApproxShBestValTry(:,i));
    AnyTies(i)=sum(DistanceBetweenEachPairof20fminconSolutions==0);
end
if(sum(AnyTies)>0)
    disp('There are repeated solutions returned by fmincon in optimizing S_{h,2} in some trials.')
end
disp('Minimum and mean over the 100 trials of the maximum L2 distance between any two of the 20 solutions returned by fmincon in optimizing S_{h,2}:')
[min(MaxDistanceBetweenEachPairof20fminconSolutions) mean(MaxDistanceBetweenEachPairof20fminconSolutions)]
disp('Maximum over the 100 trials of the range of the S_{h,2} values at the 20 solutions returned by fmincon in optimizing S_{h,2}:')
[max(RangeofApproximateShOver20fminconSolutions)]

disp('Sequence of mean(S_{h,2}(xhatstarML)-min{S_{h,2}})')
mean((store5-store3(:,end)))
disp('Sequence of mean(S_h(xhatstarML)-min{S_h})')
mean((store6-store2(:,end-1)))
% disp('Sequence of mean(||xhatstarML-ApproxxstarML||_2), where ApproxxstarML=minimizer of S_{h,2}')
% mean(store7)
% disp('Sequence of mean(||xhatstarML-xstarML||_2))')
% mean(store8)
disp('Sequence of mean of fraction of points in 41^3 grid with smaller S_{h,2} value than S_{h,2}(xhatstarML):')
mean(store9/41^3)
disp('Sequence of mean of fraction of points in 41^3 grid with smaller S_h value than S_h(xhatstarML):')
mean(store10/41^3)

disp('For Figure I.2, S_{h,2}(xhatstarML)-min{S_{h,2}}, and S_h(xhatstarML)-min{S_h} are:')
[store4(82,end)-store3(82,end) store4(82,end-1)-store2(82,end-1)]
%B=[vecnorm(store4(:,1:3)-store3(:,1:3),2,2),vecnorm(store4(:,1:3)-store2(:,1:3),2,2)];
% disp('For Figure I.2, distance between xhatstarML and ApproxxstarML (minimizer of S_{h,2}), and distance between xhatstarML and xstarML are:')
% B(82,:)

disp('Above are some results for SVD-AGP method.')
%%
disp('Evidence of convergence for three calibration methods.')
load Example2.mat
nl=16; nh=8; n=nl+nh;
nh0=12;
for i=1:100
    DistLastTwoxhatstarMLsMBCAGP(i,:)=pdist(T_MBC_AGP{i}.xhatstarMLs((end-1):end,:));
    DistLastTwoxhatstarMLsSVD(i,:)=pdist(T_SVD{i}.xhatstarMLs((end-1):end,:));
    DistLastTwoxhatstarMLsSVDAGP(i,:)=pdist(T_SVD_AGP{i}.xhatstarMLs((end-2):2:end,:));
end
disp('Median, over the 100 trials, of the Euclidean distance between last two xhatstarML''s given by the MBC-AGP method in a trial:') 
median(DistLastTwoxhatstarMLsMBCAGP)
disp('Median, over the 100 trials, of the Euclidean distance between last two xhatstarML''s given by the SVD method in a trial:') 
median(DistLastTwoxhatstarMLsSVD)
disp('Median, over the 100 trials, of the Euclidean distance between last two xhatstarML''s given by the SVD-AGP method in a trial:') 
median(DistLastTwoxhatstarMLsSVDAGP)