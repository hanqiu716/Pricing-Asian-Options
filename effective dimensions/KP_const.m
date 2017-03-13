function [d_uncor,d_poscor]=KP_const(p)
%AIM: select effective dimensions for covariance matrix A using Kronecker
%product constant volatilities case
% output : 
%     d_uncor--effective dimension for uncorrelated case
%     d_poscor--effective dimension for postice correlated case
%
%  input:
%    p--ratio

% call constant_sigma first to get the covariance R,Sigma,Sigma2
% Sigma-- uncorrelated case
% Sigma2 -- positive correlated case
[A,A2,t,Sigma,Sigma2,R]=constant_sigma;
 
 
% effective dimension selection
% compute the eigenvalues of R and Sigma respectively, 
% because A= kron(R,Sigma) in this case
% according to the property of Kronecker product, eig_SigmaMN= kron(eig(R),eig(Sigma))
% then we can get the eigenvalues of SigmaMN faster then computing the
% eig(SigmaMN) directly

% zero correlation
eig_R=eig(R);
eig_Sigma=eig(Sigma);
eig_SigmaMN=kron(eig_R,eig_Sigma);
denom=sum(eig_SigmaMN);
numer=0;
d=0;
eig_SigmaMNsort=sort(eig_SigmaMN,'descend');
 
while  (numer/denom)<=p
       d=d+1;
        numer=numer+eig_SigmaMNsort(d);
       
       
end
d_uncor=d;



% positive correlation
% effective dimension selection
eig_R2=eig(R);
eig_Sigma2=eig(Sigma2);
eig_SigmaMN2=kron(eig_R2,eig_Sigma2);
denom2=sum(eig_SigmaMN2);
numer2=0;
d2=0;
eig_SigmaMNsort2=sort(eig_SigmaMN2,'descend');
while  (numer2/denom2)<=p
        d2=d2+1; 
        numer2=numer2+eig_SigmaMNsort2(d2);
           
end

d_poscor=d2;
end
