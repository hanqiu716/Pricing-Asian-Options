function time_compare(M,N)
% aim: compare the time used by PCA and KP when p =0.9:0.01:0.99 in
% constant volatilities case
% for each percentage p, we repeatedly run PCA and KP for 10 times 
% and compute the average of time
% finally, plot time used by PCA and KP, respetively

x=zeros(10,10);
y=zeros(10,10);
j=0;
for p=0.9:0.01:0.99
j=j+1;
    for i=1:10
    tic;
    KP_const1(p,M,N);
    x(i,j)=toc;
    tic;
    PCA_const1(p,M,N);
    y(i,j)=toc;
    end
end
t_KP=sum(x)/10;
t_PCA=sum(y)/10;


xx=0.9:0.01:0.99;
figure;
plot(xx,t_PCA,'ro')
hold on
plot(xx,t_KP,'b*')
axis([0.9 0.99 0 0.5])
title('Time Comparison Between PCA and KP in the Constant Case');
ylabel('Computational Time');
xlabel('percentage p')
legend('computationl time for PCA','computational time for KP');

end