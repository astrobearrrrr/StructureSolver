example

function example
M(1,1)=30/4;
M(1,2)=0;
M(2,1)=0;
M(2,2)=10/4;

K(1,1)=3*2;
K(1,2)=-1*2;
K(2,1)=-1*2;
K(2,2)=1*2;

C(1,1)=0;
C(1,2)=0;
C(2,1)=0;
C(2,2)=0;


tend=1.2;
dt=0.12;

for i=1:1:floor(tend/dt)+1
d1(1,i)=0;
d1(2,i)=0;
v1(1,i)=0;
v1(2,i)=0;
f(1,i)=0;
f(2,i)=1;
end

[d,v,a] = Newmark( K, M, C, f, d1, v1, dt, tend );
% d位移，v载荷，a加速度
fprintf('Newmark 节点位移 %12.8f %12.8f\n',d)

return
end

function [d,v,a] = Newmark( K, M, C, f, d1, v1, dt, tend )
% 利用Newmark 法计算结构的动力响应
% [d,v,a] = Newmark( K, M, C, f, d1, v1, dt, tend )
% 输入参数
% K ----- 刚度矩阵
% M ----- 质量矩阵
% C ----- 阻尼矩阵
% d1 ----- 初始位移
% v1 ----- 初始速度
% dt ----- 时间步长
% tend --- 结束时间
% 返回值
% d ----- 位移
% v ----- 速度
% a ----- 加速度

delta=0.5;
alpha=0.25;

[n,n] = size( K ) ;

alpha0 = 1/alpha/dt^2 ;
alpha1 = delta/alpha/dt ;
alpha2 = 1/alpha/dt ;
alpha3 = 1/2/alpha - 1 ;
alpha4 = delta/alpha - 1 ;
alpha5 = dt/2*(delta/alpha-2) ;
alpha6 = dt*(1-delta) ;
alpha7 = delta*dt ;

K1 = K + alpha0*M + alpha1*C ;
d = zeros( n, floor(tend/dt) + 1 ) ;
v = zeros( n, floor(tend/dt) + 1 ) ;
a = zeros( n, floor(tend/dt) + 1 ) ;
d(:,1) = d1(:,1) ;
v(:,1) = v1(:,1) ;
a(:,1) = M\(f(:,1)-K*d1(:,1)-C*v1(:,1)) ;

    outfile=fopen('out.txt','w');
    fprintf( outfile,'=========   Newmark method   =================\n' ) ;
    fprintf( outfile,'    Time point         Node 2       Node 3 \n' ) ;
time(1)=0;
for i=2:1:floor(tend/dt) + 1
t = (i-1)*dt ;
if t>0.5
  f(:,i)=0.0;
end
f1 = f(:,i) + M*(alpha0*d(:,i-1)+alpha2*v(:,i-1)+alpha3*a(:,i-1)) ...
+ C*(alpha1*d(:,i-1)+alpha4*v(:,i-1)+alpha5*a(:,i-1)) ;
d(:,i) = K1\f1 ;
a(:,i) = alpha0*(d(:,i)-d(:,i-1)) - alpha2*v(:,i-1) - alpha3*a(:,i-1) ;
v(:,i) = v(:,i-1) + alpha6*a(:,i-1) + alpha7*a(:,i) ;

fprintf(outfile, '%12.4f       %12.8f %12.8f\n', t,d(1,i),d(2,i)) ;
time(i)=t;
end

    figure ;
    plot(time,d(1,:),'-', time, d(2,:), ':') ;
    text(0.8,0.2,'Newmark 节点2') 
    text(0.4,0.4,'Newmark 节点3') 

return
end