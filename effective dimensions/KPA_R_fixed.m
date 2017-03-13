function [dd_uncor,dd_poscor]=KPA_R_fixed(p)
% AIM: select effective dimensions by KPA in time-dep case, only R is fixed 
% output: 
%     dd_uncor--effective dimension for uncorrelated case
%     dd_poscor--effective dimension for postice correlated case
%
%  input:
%    p--ratio
% 

% call function time_dep_sigma first
[A,A2,t]=time_dep_sigma;

% Assume R is known we solve the least square problem: norm(Sigma-kron(R,K))
% RR is the covarianc matrix of a single brownian motion, it is constant
% and symmetric, and has a boomerang shape
% KK is the matrix that we want to compute
% We obtain a formula for optimum KK in Pitsisanis and Van Loan's paper

% when rhol= 0;
Sigma=cell2mat(A);


RR=zeros(250);
for f=1:250
    RR(f,f:end)=t(f);
    RR(f:end,f)=t(f);
end

KK=zeros(10);

btrace=trace(RR'*RR);
for pp=1:10
    for q=1:10
        Sigmahat=Sigma(pp:10:2500,q:10:2500);
        KK(pp,q)= trace(Sigmahat'*RR)/btrace;
    end
end

% effective dimension selection
% compute the eigenvalue of R and K, respectively
% since in this case A is approx equal to kron(R,K), 
% eig(A)is approxiamtely equal to kron(eig_R, eig_K)
% sort the eigenvalues of A into descending order, from largest to smallest
% let the sum of first d eigenvalues of A be the numerator,
% and sum(eig(A)) is denomenator,
% then compute d such that numer/denom not greater than ratio p
eig_RR=eig(RR);
eig_KK=eig(KK);
eig_sigma=kron(eig_RR,eig_KK);
denomm=sum(eig_sigma);
numerr=0;
d=0;
eig_sigmasort=sort(eig_sigma,'descend');
 
while  (numerr/denomm)<=p
       d=d+1;
        numerr=numerr+eig_sigmasort(d)  ; 
end
dd_uncor=d;
    


% when rhol= 0.4;
Sigma2=cell2mat(A2);

RR2=RR;

KK2=zeros(10);

btrace=trace(RR'*RR);
for pp=1:10
    for q=1:10
        Sigmahat2=Sigma2(pp:10:2500,q:10:2500);
        KK2(pp,q)= trace(Sigmahat2'*RR)/btrace;
    end
end

% effective dimension selection
eig_RR2=eig(RR2);
eig_KK2=eig(KK2);
eig_sigma2=kron(eig_RR2,eig_KK2);
denomm2=sum(eig_sigma2);
numerr2=0;
d2=0;
eig_sigmasort2=sort(eig_sigma2,'descend');
while  (numerr2/denomm2)<=p
       d2=d2+1;
        numerr2=numerr2+eig_sigmasort2(d2);    
       
end

dd_poscor=d2;
end 
