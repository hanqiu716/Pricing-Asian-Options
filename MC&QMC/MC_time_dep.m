function [a_uncor,a_pos]=MC_time_dep(n)
% AIM :estimate  the asian options price with Monte Carlo Simulations
% this is time-dependent volatilities situation
%

% output: 
%     a_uncor= Asian options price in zero correlation case
%     a_pos= Asian options price in positive correlation case
%
%intput: 
%     n--number of simulations we need 
%
%

% w_ij represents the weight of each asset i at time j, sum (w)=1, i =1:M, j=1:N
% run the fucntion time_dep_sigma to get the covariance matrix of Z

[A,A2,t]=time_dep_sigma;
Sigma=cell2mat(A);
Sigma2=cell2mat(A2);
K=100;
T=1;
r=0.04;
S=100;
M=10;
N=250;
W=rand(M,N);
ww=sum(sum(W),2);
w=W/ww;
miu=zeros(M*N,1);
for k =1:M*N
    k1=mod(k-1,M)+1;
    k2=floor((k-1)/M)+1;
    sigmahat(k1)=0.1+((k1-1)/9)*0.4-0.09;
    sigma=@(t)(sigmahat(k1)*exp(-t/1.5)+0.09).^2;
    miu(k)=log(w(k1,k2)*S)+r*t(k2)-integral(sigma,0,t(k2))/2;
end
% approximate the integral of a function f as the average of the function
% evaluated at a set of points u1,u2...un pseudorandom sequence
% Since we are integrating over the MN-dimensional unit cube,
% each ui is a vector of M*N elements.
% repeat random sampling for n times

% zero correlation case
SUM=0;
for i=1:n
 u=rand(M*N,1);
 Z=norminv(u,zeros(M*N,1),sqrt(diag(Sigma)));

 g=sum(exp(miu+Z));

 f=max((g-K),0);
 SUM=SUM+f;
 integral_time_uncor=SUM/n;
end

% postive correlation case
for i=1:n
 u=rand(M*N,1);
 Z2=norminv(u,zeros(M*N,1),sqrt(diag(Sigma2)));

 g2=sum(exp(miu+Z2));

 f2=max((g2-K),0);
 SUM=SUM+f2;
 integral_time_pos=SUM/n;
end

% evaluate asian optiosn price
a_uncor=exp(r*(T-0))*integral_time_uncor;
a_pos=exp(r*(T-0))*integral_time_pos;


end
