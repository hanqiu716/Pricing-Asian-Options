function [d_uncor,d_poscor]=PCA_const(p)
%AIM: select effective dimensions for covariance matrix A in constant case, using principal
%components analysis
% output : 
%     d_uncor--effective dimension for uncorrelated case
%     d_poscor--effective dimension for postice correlated case
%
%  input:
%    p--ratio

% call constant_sigma first to get the covariance A and A2
[A,A2]=constant_sigma;
 
 
% zero correlation
% effective dimension selection
eig_A=eig(A);
denom=sum(eig_A);
numer=0;
d=0;
eig_Asort=sort(eig_A,'descend');
 
while  (numer/denom)<=p
       d=d+1;
        numer=numer+eig_Asort(d);
       
       
end
d_uncor=d;



% positive correlation
% effective dimension selection

eig_A2=eig(A2);
denom2=sum(eig_A2);
numer2=0;
d2=0;
eig_Asort2=sort(eig_A2,'descend');
while  (numer2/denom2)<=p
        d2=d2+1; 
        numer2=numer2+eig_Asort2(d2);
           
end

d_poscor=d2;
end
