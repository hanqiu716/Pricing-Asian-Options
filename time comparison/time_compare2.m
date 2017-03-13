function time_compare2(M,N)
% aim: compare the time used by PCA and KPA when p =0.9:0.01:0.99 in
% time-dep volatilities case
% for each percentage p, we repeatedly run PCA and KPA for 10 times 
% and compute the average of time
% finally, plot time used by PCA and KPA, respetively

x=zeros(10,10);
y=zeros(10,10);
z=zeros(10,10);
j=0;
for p=0.9:0.01:0.99
j=j+1;
    for i=1:10
    tic;
    PCA_time_dep1(p,M,N);
    x(i,j)=toc;
    tic;
    KPA_brute_method1(p,M,N);
    y(i,j)=toc;
    tic;
    KPA_R_fixed1(p,M,N);
    z(i,j)=toc;
    end
end

t_PCA=sum(x)/10;
t_KPA1=sum(y)/10;
t_KPA2=sum(z)/10;


xx=0.9:0.01:0.99;
plot(xx,t_PCA,'ro')
hold on
plot(xx,t_KPA1,'b*')
hold on 
plot(xx,t_KPA2,'ko')
axis([0.9 0.99 0 0.5])
title('Time Comparison Between PCA and KPA in the Time-dependent Case');
ylabel('Computational Time');
xlabel('Ratio p')
legend('computational time for PCA','computationl time for KPA--full approximation','computationl time for KPA--R fixed');
end
