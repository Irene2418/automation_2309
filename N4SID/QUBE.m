clc
clear
A_d=[0,1;-0.9799,1.9799];%2*2
B_d=[0;9.5694E-4];%2*1
C_d=[1,0];%1*2
D_d=0;%1*1
NN=10000;
%[A_d,B_d,C_d,D_d]=c2dm(A,B,C,D,0.1);%连续时间系统离散化
A_d %n*n
B_d %n*m
C_d %l*n
D_d %l*m

% N4SID算法
m = 1;%输入维数
l = 1;%输出维数
i = 50;  %＞n
j = 1000; %＞＞n

% 构建状态和输出矩阵
%u=1000*idinput(NN, 'rbs');
u=1000*wgn(NN,1,100);
%u = 10000*randn(1, NN);
x=zeros(2, NN+1);
y=zeros(1, NN);
for k=1:NN
    if k==1%time sequence 0
        x(:,k)=[1;1];
    end
    t(k)=k-1;
    %u(k)=100*sin(1*k);
    %    for s=2:(2*i*m)
    %       u(k)=u(k)+100*sin(s*k);
    %    end    
    x(:,k+1)=A_d*x(:,k)+B_d*u(k);
    y(k)=C_d*x(:, k);
end

%构造Hankel矩阵
Yp=zeros(i*l,j);
Up=zeros(i*m,j);
Yf=zeros(i*l,j);
Uf=zeros(i*m,j);
for k1=1:i*l
    for k2=1:j
        Yp(k1,k2)=y(k1+k2-1);
        Yf(k1,k2)=y(i+k1+k2-1);
    end
end
for k1=1:i*m
    for k2=1:j
        Up(k1,k2)=u(k1+k2-1);
        Uf(k1,k2)=u(i+k1+k2-1);
    end
end
Wp=[Up;Yp];

%根据定义求斜向投影Yf/UfWp
%UWP=pinv([Wp;Uf]*[Wp' Uf']);%斜向投影矩阵对应的大矩阵?
%Oi=Yf*[Wp' Uf']*UWP(:,1:i*m+i*l)*Wp;%Oi=Yf/UfWp
[Oid,Oilq]=obp(Uf,Wp,Yf);
[UU,SS,VV]=svd(Oilq);%奇异值分解，UU=[U1,U2]，SS=diag[S1 S2],VV=[V1^T,V2^T]^T
n=2;%人工核验SS主导奇异值的阶次
U1=UU(:,1:n);
S1=SS(1:n,1:n);
T=[1 2;0 1];%任意一个可逆矩阵?
VVT=VV';%V的转置?
V1T=VVT(1:n,:);%V1的转置?
Gammai=U1*sqrtm(S1)*T;
Xf=inv(T)*sqrtm(S1)*V1T;
Xip1jm1=Xf(:,2:j);%X_{i+1,j-1}
Xijm1=Xf(:,1:j-1);%X_{i,j-1}
Uijm1=zeros(m,j-1);%U_{i,j-1}
Yijm1=zeros(l,j-1);%Y_{i,j-1}
for k1=1:j-1
    Uijm1(:,k1)=u(i+k1);%U_{i,j-1}具体赋值??
    Yijm1(:,k1)=y(i+k1);%Y_{i,j-1}具体赋值??
end
hatABCD=[Xip1jm1;Yijm1]*[Xijm1;Uijm1]'*inv([Xijm1;Uijm1]*[Xijm1;Uijm1]');%利用最小二乘求解[A B;C D]矩阵
AT=hatABCD(1:n,1:n)%hatABCD的前n行n列为A矩阵
BT=hatABCD(1:n,n+1:n+m)
CT=hatABCD(n+1:n+l,1:n)
DT=hatABCD(n+1:n+l,n+1:n+m)
xT=zeros(2,NN+1);
for k=1:NN
    xT(:,k+1)=AT*xT(:,k)+BT*u(k);
    yT(k)=CT*xT(:,k)+DT*u(k);
    t(k)=k-1;
end
%绘制前200个数据
ltext={'y(k)','yT(k)'};
h(1,1)=plot(t(1:200),y(1:200),'-b','LineWidth',2);hold on;
h(2,1)=plot(t(1:200),yT(1:200),'--r','LineWidth',2);hold off;
legend(h,ltext);
xlabel('time sequence k')
ylabel('y(k) and yT(k)')

function [obp_ABCd,obp_ABClq]=obp(A,B,C)
obp_ABCd=A*[C' B']*pinv([C;B]*[C' B'])*[C;0*B];
[L11,L21,L22,L31,L32,L33,Q1T,Q2T,Q3T]=lqBCA(B,C,A);
obp_ABClq=L32*inv(L22)*[L21 L22]*[Q1T;Q2T];
end

function [L11,L21,L22,L31,L32,L33,Q1T,Q2T,Q3T]=lqBCA(B,C,A)
kB=size(B,1);
kC=size(C,1);
[Q,R]=qr([B;C;A]');
QT=Q';
L=R';%[B;C;A]=L*QT
L11=L(1:kB,1:kB);
L21=L(kB+1:kB+kC,1:kB);
L22=L(kB+1:kB+kC,kB+1:kB+kC);
sLr=size(L,1);
sLc=size(L,2);
L31=L(kB+kC+1:sLr,1:kB);
L32=L(kB+kC+1:sLr,kB+1:kB+kC);
L33=L(kB+kC+1:sLr,kB+kC+1:sLc);
Q1T=QT(1:kB,:);
Q2T=QT(kB+1:kB+kC,:);
sQTr=size(QT,1);
Q3T=QT(kB+kC+1:sQTr,:);
end