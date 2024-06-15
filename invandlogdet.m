function [invCovMat, logdetCovMat,condCovMat]=invandlogdet(CovMat) 

    [Eigvec, Eigval]=svd(CovMat);
    Eigval=diag(Eigval);
    invCovMat=Eigvec*diag(1./Eigval)*Eigvec';
    logdetCovMat=sum(log(Eigval));
    condCovMat=max(Eigval)/min(Eigval);
end