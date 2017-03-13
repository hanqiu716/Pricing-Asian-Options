function [A,A2,t,Sigma,Sigma2,R]=constant_sigma
% Generate the Matrix A in constant case
% output:
%      A--convariance matrix in zero correlated case
%      A2--convariance matrix in positive correlated case
%      Sigma--covariance matrix of asset returns in zero correlated case
%      Sigma2--covariance matrix of asset returns in positive correlated case
%      R--autocovariance matrix of each Brownian motion


% A is a 2500-by-2500 matrix, 
% A= kron(R,B)
% R is the convariance matrix of each Brownian motion and has a boomerang form
% B is a 10-by-10 matrix, B=sigma_i*sigma_k

t=1/250:1/250:1;
R=zeros(250);
for f=1:250
    R(f,f:end)=t(f);
    R(f:end,f)=t(f);
end
sigma=zeros(10,1);
for n=1:10
  
         sigma(n)=0.1+(n-1)/9*0.4;
    
end

for i=1:10
    for k =1:10
        B(i,k)=sigma(i)*sigma(k);
    end
end
% when rho=0, uncorrelated case
   rho=eye(10);
    Sigma=rho.*B;
    A=kron(R,Sigma);
    
% when rho=0.4, positive correlated case
    rho2=0.4*ones(10)+0.6*eye(10);
    Sigma2=rho2.*B;
    A2=kron(R,Sigma2);
end
