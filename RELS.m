%递推增广最小二乘参数辨识RELS
clear all; 
close all;
clc

a=[1 -1.2 0.8]'; b=[1 0.6]'; c=[1 -1 0.2]'; %对象参数
d=3; 
na=length(a)-1; nb=length(b)-1; nc=length(c)-1; %na、nb、nc分别为为A、B、C阶次

L=1000; %仿真长度
uk=zeros(d+nb,1); %输入初值
yk=zeros(na,1); %输出初值
xik=zeros(nc,1); %噪声初值
xiek=zeros(nc,1); %噪声估计初值
u=randn(L,1); %输入采用白噪声序列
xi=sqrt(0.1)*randn(L,1); %白噪声序列

theta=[a(2:na+1);b;c(2:nc+1)]; %对象参数

thetae_1=zeros(na+nb+1+nc,1); %thetae初值
P=10^6*eye(na+nb+1+nc);
for k=1:L
    phi=[-yk;uk(d:d+nb);xik];
    y(k)=phi'*theta+xi(k); %采集输出数据
    
    phie=[-yk;uk(d:d+nb);xiek]; 
    
    %递推增广最小二乘法
    K=P*phie/(1+phie'*P*phie);
    thetae(:,k)=thetae_1+K*(y(k)-phie'*thetae_1);
    P=(eye(na+nb+1+nc)-K*phie')*P;
    
    xie=y(k)-phie'*thetae(:,k); %白噪声的估计值
    
    %更新数据
    thetae_1=thetae(:,k);
    
    for i=d+nb:-1:2
        uk(i)=uk(i-1);
    end
    uk(1)=u(k);
    
    for i=na:-1:2
        yk(i)=yk(i-1);
    end
    yk(1)=y(k);
    
    for i=nc:-1:2
        xik(i)=xik(i-1);
        xiek(i)=xiek(i-1);
    end
    xik(1)=xi(k);
    xiek(1)=xie;
end
figure('Name','RLS:a','NumberTitle','off')
plot([1:L],thetae(1:na,:));
xlabel('k'); ylabel('参数估计值');
legend('a1','a2'); axis([0 L -2 2]);
figure('Name','RLS:b','NumberTitle','off')
plot([1:L],thetae(na+1:na+nb+1,:));
xlabel('k'); ylabel('参数估计值');
legend('b0','b1'); axis([0 L 0 1.5]);
figure('Name','RLS:c','NumberTitle','off')
plot([1:L],thetae(na+nb+2:na+nb+nc+1,:));
xlabel('k'); ylabel('参数估计值');
legend('c1','c2'); axis([0 L -2 2]);
disp('参数估计值:');
disp(thetae_1);