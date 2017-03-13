function [A,A2,t]=time_dep_sigma
% aim: build up the covariance matrix A in time-dependet volatilities case
% output:  
%               A, - covariance matrix for uncorrelated case; 
%               A2 - covariance matrix for positive correlated case; 

% A is a 250 by 250 block matrix, each block is 10 by 10
% t is the time grid has 250 equally spaced point in (0,1]
% sigmainf is the asymptotic volatility
% tau is the dacey constant£¬ equals to 1.5
% all block matrices in the first row and column of A are the same
% all  block matrices in the second row A{2,2:end} and column A{2:end,2} are
% the same
% ....
% A{i,i},A{i,i+!}...A{i,end}£¬A{i+1£¬i}...A{end,i} are the same, i=1,2...250
%
%

clear all
t=1/250:1/250:1;

sigmainf=0.09;
taui=1.5;
tauk=1.5;
tauik=taui*tauk/(taui+tauk);
B=zeros(10,10);
A=cell(250,250);
for m=1:250
    for i=1:10
        for k=1:10
            sigmai0=0.1+(i-1)/9*0.4;
            sigmak0=0.1+(k-1)/9*0.4;
            sigmahati0=sigmai0-sigmainf;
            sigmahatk0=sigmak0-sigmainf;
            B(i,k)=sigmahati0*sigmahatk0*tauik*(1-exp(-t(m)/tauik))+...
            sigmahati0*sigmainf*taui*(1-exp(-t(m)/taui))+...
            sigmahatk0*sigmainf*tauk*(1-exp(-t(m)/tauk))+...
            sigmainf*sigmainf*t(m);
        end
    end
    rho=eye(10);
    brho=rho.*B;
    rho2=0.4*ones(10)+0.6*eye(10);
    brho2=rho2.*B;
    for n=250:-1:m
        A{m,n}=brho;
        A{n,m}=brho;
        A2{m,n}=brho2;
        A2{n,m}=brho2;
    end
end
end







