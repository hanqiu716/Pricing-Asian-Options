function [d_uncor,d_poscor]=PCA_time_dep(p)
%AIM: select effective dimensions for covariance matrix A in time-dep case, using principal
%components analysis
% output: 
%     d_uncor--effective dimension for uncorrelated case
%     d_poscor--effective dimension for postice correlated case
% input:
%    p--ratio

% call time_dep_sigma first to get the covariance A and A2
[A,A2]=time_dep_sigma;


% convert cell matrix A and A2 into real matrix Sigma and Sigma2
Sigma=cell2mat(A);
Sigma2=cell2mat(A2);

 
% zero correlation
% effective dimension selection
eig_Sigma=eig(Sigma);
denom=sum(eig_Sigma);
numer=0;
d=0;
eig_Sigmasort=sort(eig_Sigma,'descend');
 
while  (numer/denom)<=p
       d=d+1;
        numer=numer+eig_Sigmasort(d);
       
end
d_uncor=d;



% positive correlation
% effective dimension selection
eig_Sigma2=eig(Sigma2);
denom2=sum(eig_Sigma2);
numer2=0;
d2=0;
eig_Sigmasort2=sort(eig_Sigma2,'descend');
while  (numer2/denom2)<=p
        d2=d2+1; 
        numer2=numer2+eig_Sigmasort2(d2);
           
end

d_poscor=d2;
end
