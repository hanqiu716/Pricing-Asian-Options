function [a_uncor,a_pos]=MC_const(n)
%
% AIM: estimate  the asian options price with Monte Carlo Simulations in constant situation
%
% output: 
%     a_uncor= Asian options price in zero correlation case
%     a_pos= Asian options price in positive correlation case
%
%intput: 
%     n--number of simulations we need 


% run the fucntion constant_sigma to get the covariance matrix of Z

[A,A2,t]=constant_sigma;
K=100;
T=1;
r=0.04;
S=100;
M=10;
N=250;
sigma=zeros(10,1);
% w_ij represents the weight of each asset i at time j, sum (w)=1, i =1:M, j=1:N
W=rand(M,N);
ww=sum(sum(W),2);
w=W/ww;
miu=zeros(M*N,1);
for k =1:M*N
    k1=mod(k-1,M)+1;
    k2=floor((k-1)/M)+1;
    sigma(k1)=0.1+((k1-1)/9)*0.4;
    miu(k)=log(w(k1,k2)*S)+(r-sigma(k1)^2/2)*t(k2);
end

% approximate the integral of function f as the average of the function
% evaluated at a set of points u1,u2...un form pseudorandom sequence
% Since we are integrating over the MN-dimensional unit cube,
% each ui is a vector of M*N elements.
% repeat random sampling for n times

SUM=0;

% uncorrelated
for i=1:n
 u=rand(M*N,1);
 Z=norminv(u,zeros(M*N,1),sqrt(diag(A)));

 g=sum(exp(miu+Z));

 f=max((g-K),0);
 SUM=SUM+f;
 integral_const_uncor=SUM/n;
end

% postive correlated
for i=1:n
 u=rand(M*N,1);
 Z2=norminv(u,zeros(M*N,1),sqrt(diag(A2)));

 g2=sum(exp(miu+Z2));

 f2=max((g2-K),0);
 SUM=SUM+f2;
 integral_const_pos=SUM/n;
end

% evaluate asian options price
a_uncor=exp(r*(T-0))*integral_const_uncor;
a_pos=exp(r*(T-0))*integral_const_pos;

end