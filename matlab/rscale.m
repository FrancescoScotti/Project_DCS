function[Nbar]=rscale(A,B,C,D,K)
    dim = size(A,1);
    Z = [zeros([1,dim]) 1];
    N = inv([A,B;C,D])*Z';
    Nx = N(1:dim);
    Nu = N(1+dim);
    Nbar=Nu + K*Nx;
end