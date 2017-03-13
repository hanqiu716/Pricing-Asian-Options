function [d_uncor,d_poscor]=KPA_brute_method(p)
% AIM: select effective dimensions by KPA in time-dep case, both R and K are not fixed

% output : 
% dd_uncor--effective dimension for uncorrelated case
%     dd_poscor--effective dimension for postice correlated case
%
%  input:
%    p--ratio
%


% call function time_dep_sigma first to get the covariance matrix A, A2
[A,A2]=time_dep_sigma;

% when rhol=0,rewrite A as tilde(A), named as S, S is a 250^2 * 100 matrix, 
%  when rhol=0.4, rewrite A2 as tilde(A2), named as S2
d=1;
for l=1:250
    for j=1:250
        S(d,:)=(A{j,l}(:))';
        S2(d,:)=(A2{j,l}(:))';
        d=d+1;
    end
end

% Do the SVD decompostion for S, just get the first column of U and V, and
% the first entry of diagonal matrix
% make the first column of U and V are positive
% r = sqrt(D)*U1; is a 250^2by 1 matrix
% k=sqrt(D)*V1;is a 250^2by 1 matrix
% reshape r and k to a 250 by 250 square matrix
% kron(R,K) is approximately to A

% rhol=0, uncorrelated
[U,D,V]=svds(S,1);
U1=-U;
V1=-V;
r=sqrt(D)*U1;
k=sqrt(D)*V1;
R=reshape(r,250,250);
K=reshape(k,10,10);

% effective dimension selection
eig_R=eig(R);
eig_K=eig(K);
eig_sigmaMN=kron(eig_R,eig_K);
denom=sum(eig_sigmaMN);
numer=0;
d=0;
eig_sigmaMNsort=sort(eig_sigmaMN,'descend');

    while  (numer/denom)<=p
        d=d+1;
        numer=numer+eig_sigmaMNsort(d);
    end

d_uncor=d;



% rhol= 0.4, positive correlated
%
[U2,D2,V2]=svds(S2,1);
U3=-U2;
V3=-V2;
r2=sqrt(D2)*U3;
k2=sqrt(D2)*V3;
R2=reshape(r2,250,250);
K2=reshape(k2,10,10);

% effective dimension selection
eig_R2=eig(R2);
eig_K2=eig(K2);
eig_sigmaMN2=kron(eig_R2,eig_K2);
denom2=sum(eig_sigmaMN2);
numer2=0;
d2=0;
eig_sigmaMN2sort=sort(eig_sigmaMN2,'descend');
while (numer2/denom2)<=p
    d2=d2+1;
     numer2=numer2+eig_sigmaMN2sort(d2);
end
d_poscor=d2;
    
end








