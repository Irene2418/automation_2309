clc
clear
A=[0,0,1,0 ; 0,0,0,1 ; 0,149.2751,-14.9287,4.9149;0,-261.6091,14.7551,-8.6136];%4*4
B=[0;0;49.7275;-49.1493];%4*1
C=[1,0,0,0;0,1,0,0];%2*4
D=[0;0];%2*1
NN=10000;
[A_d,B_d,C_d,D_d]=c2dm(A,B,C,D,0.1);%����ʱ��ϵͳ��ɢ��
A_d %n*n
B_d %n*m
C_d %l*n
D_d %l*m

% N4SID�㷨
m = 1;%����ά��
%n = 4;״̬ά��
l = 2;%���ά��
i = 6;  %��n
j = 3500; %����n

% ����״̬���������
%u=idinput(NN, 'rbs');

x=zeros(4, NN+1);
y=zeros(2, NN);
for k=1:NN
    if k==1%time sequence 0
        x(:,k)=[0;0;0;0];
    end
    t(k)=k-1;
    u(k)=sin(1*k);
        for s=2:(2*i*m)
           u(k)=u(k)+sin(s*k);
       end    
    x(:,k+1)=A_d*x(:,k)+B_d*u(k);
    y(:,k)=C_d*x(:, k);
end

%����Hankel����
Yp=zeros(i*l,j);
Up=zeros(i*m,j);
Yf=zeros(i*l,j);
Uf=zeros(i*m,j);
for k1=1:2:i*l
    for k2=1:j
        Yp(k1,k2)=y(1,(k1+1)/2+k2-1);
        Yp(k1+1,k2)=y(2,(k1+1)/2+k2-1);
        Yf(k1,k2)=y(1,i+(k1+1)/2+k2-1);
        Yf(k1+1,k2)=y(2,i+(k1+1)/2+k2-1);
    end
end
for k1=1:i*m
    for k2=1:j
        Up(k1,k2)=u(k1+k2-1);
        Uf(k1,k2)=u(i+k1+k2-1);
    end
end
Wp=[Up;Yp];

%���ݶ�����б��ͶӰYf/UfWp
UWP=pinv([Wp;Uf]*[Wp' Uf']);%б��ͶӰ�����Ӧ�Ĵ����?
Oi=Yf*[Wp' Uf']*UWP(:,1:i*m+i*l)*Wp;%Oi=Yf/UfWp

[UU,SS,VV]=svd(Oi);%����ֵ�ֽ⣬UU=[U1,U2]��SS=diag[S1 S2],VV=[V1^T,V2^T]^T
n=4;%�˹�����SS��������ֵ�Ľ״�
U1=UU(:,1:n);
S1=SS(1:n,1:n);
T=[1 2 3 4;0 1 0 0; 0 0 1 0; 0 0 0 1];%����һ���������?
VVT=VV';%V��ת��?
V1T=VVT(1:n,:);%V1��ת��?
Gammai=U1*sqrtm(S1)*T;
Xf=inv(T)*sqrtm(S1)*V1T;
Xip1jm1=Xf(:,2:j);%X_{i+1,j-1}
Xijm1=Xf(:,1:j-1);%X_{i,j-1}
Uijm1=zeros(m,j-1);%U_{i,j-1}
Yijm1=zeros(l,j-1);%Y_{i,j-1}
for k1=1:j-1
    Uijm1(:,k1)=u(i+k1);%U_{i,j-1}���帳ֵ??
    Yijm1(:,k1)=y(:,i+k1);%Y_{i,j-1}���帳ֵ??
end
hatABCD=[Xip1jm1;Yijm1]*[Xijm1;Uijm1]'*inv([Xijm1;Uijm1]*[Xijm1;Uijm1]');%������С�������[A B;C D]����
AT=hatABCD(1:n,1:n)%hatABCD��ǰn��n��ΪA����
BT=hatABCD(1:n,n+1:n+m)
CT=hatABCD(n+1:n+l,1:n)
DT=hatABCD(n+1:n+l,n+1:n+m)
xT=zeros(4,NN+1);
for k=1:NN
    xT(:,k+1)=AT*xT(:,k)+BT*u(k);
    yT(:,k)=CT*xT(:,k)+DT*u(k);
    t(k)=k-1;
end
%����ǰ200������
ltext={'y(k)','yT(k)'};
h(1,1)=plot(t(1:200),y(1:200),'-b','LineWidth',2);hold on;
h(2,1)=plot(t(1:200),yT(1:200),'--r','LineWidth',2);hold off;
legend(h,ltext);
xlabel('time sequence k')
ylabel('y(k) and yT(k)')



