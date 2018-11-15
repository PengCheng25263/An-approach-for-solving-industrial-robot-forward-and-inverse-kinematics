%project final.m
%������Ŀ�϶࣬����Ҫ���� ��Ҫ���еĲ��� �������������У�����ֱ�ӵ����С�

%��һ���������⣺
DR=pi/180
alpha=[0 -90*DR 0 -90*DR 90*DR -90*DR]
a=[0 312 1075 225 0 0]
d=[0 0 0 1280 0 0]
theta=[0 0 0 0 0 0]
%theta=[10*DR 20*DR 30*DR 10*DR -30*DR 10*DR]
%theta=[90*DR -35*DR 79*DR -80*DR 10*DR 120*DR]

A1=rotx(alpha(1))*transl(a(1),0,0)*rotz(theta(1))*transl(0,0,d(1))
A2=rotx(alpha(2))*transl(a(2),0,0)*rotz(theta(2))*transl(0,0,d(2))
A3=rotx(alpha(3))*transl(a(3),0,0)*rotz(theta(3))*transl(0,0,d(3))
A4=rotx(alpha(4))*transl(a(4),0,0)*rotz(theta(4))*transl(0,0,d(4))
A5=rotx(alpha(5))*transl(a(5),0,0)*rotz(theta(5))*transl(0,0,d(5))
A6=rotx(alpha(6))*transl(a(6),0,0)*rotz(theta(6))*transl(0,0,d(6))

%---------------------------------------------------------------------
%����syms pi �� clear pi �����������к������ϲ��ֵĴ���

%�����⣺
B1=A1*A2*A3   %T��03��
B2=A1*A2*A3*A4*A5*A6   %T��06��

%�����⣺
sigma=0   %0��ʾ�˶��ؽ�
L1=link([alpha(1),a(1),theta(1),d(1),sigma],'mod')
L2=link([alpha(2),a(2),theta(2),d(2),sigma],'mod')
L3=link([alpha(3),a(3),theta(3),d(3),sigma],'mod')
L4=link([alpha(4),a(4),theta(4),d(4),sigma],'mod')  
L5=link([alpha(5),a(5),theta(5),d(5),sigma],'mod')
L6=link([alpha(6),a(6),theta(6),d(6),sigma],'mod')
R=robot({L1,L2,L3,L4,L5,L6},'Rlan3R')    
q=[0 0 0 0 0 0] 
%q=theta
TR=fkine(R,q)
drivebot(R,q)

%�������ߡ����� ����⣬��֤��

%T1�ɽ⣬�������⣻T2Ϊ�ٽ�⣬��һ��,����ʱӦ��theta2,3���¸�ֵ��ȥ���������֣�T3�޽⣬���ڻ����˽Ӵ���Χ�ڡ�
%T1% T=[0.11013 0.52562 0.84356 1604.7;-0.96534 -0.1455 0.21668 926.49;0.23663 -0.83819 0.49138 1569.1;0 0 0 1]
%T2% T=[-0.73794 0.57972 0.34551 1499.3;-0.63372 -0.77128 -0.059391 0;0.23205 -0.26278 0.93654 2056.5;0 0 0 1]
%T3% T=[0.1658 -0.1736 -0.9708 2655;0.0292 0.9848 -0.1712 866.5;0.9857 0 0.1683 806.3;0 0 0 1]

%��֤��q = ikine(R, T)    T1�⣺q =[0.5236 -1.2217 -0.6196  0.7541 -0.3490 0.4362]
%����ikine���޷����

x=T(1,4)
y=T(2,4)
z=T(3,4)
Len1=312
Len2=1075
Len3=sqrt(1280^2+225^2)
Len=sqrt(x^2+y^2)-Len1
cmp1=sqrt(z^2+Len^2)
cmp2=Len2+Len3


%��theta1
theta1=atan2(y,x) 

%��theta3
if(cmp1>cmp2) 
    
    fprintf('out of reach');%T2�ڴ˴��޷�ͨ������CMP1��CMP2����С��ֱ������theta2 theta3����ʽ����ȡ����е�ʵ�����ֶ������¸�ֵ��
else
    theta3=acos((z^2+Len^2-Len2^2-Len3^2)/(2*Len2*Len3))-atan2(1280,225)
    fprintf('theta3=%d\n',theta3 )
end
%theta3_2=-(acos((z^2+Len^2-Len2^2-Len3^2)/(2*Len2*Len3))+atan2(1280,225))

%��theta2
beta=atan2(z,Len) 
phi=acos((Len^2+z^2+Len2^2-Len3^2)/(2*Len2*sqrt(Len^2+z^2)))
theta2=-(beta+phi)
%theta2_2=-(beta-phi)

%syms pi %�򻯾���,���ڹ۲��ϵ���������̡���ѡ�ã����ᵼ�»�����ģ���޷�������
%clear pi %���pi����������һ�������ⲿ�ֵĴ���ſɽ���������ģ�͡�
theta(1)=theta1
theta(2)=theta2
theta(3)=theta3
%��һ���
% theta(2)=theta2_2
% theta(3)=theta3_2

syms theta4 theta5 theta6
A1=rotx(alpha(1))*transl(a(1),0,0)*rotz(theta(1))*transl(0,0,d(1))
A2=rotx(alpha(2))*transl(a(2),0,0)*rotz(theta(2))*transl(0,0,d(2))
A3=rotx(alpha(3))*transl(a(3),0,0)*rotz(theta(3))*transl(0,0,d(3))
A4=rotx(alpha(4))*transl(a(4),0,0)*rotz(theta4)*transl(0,0,d(4))
A5=rotx(alpha(5))*transl(a(5),0,0)*rotz(theta5)*transl(0,0,d(5))
A6=rotx(alpha(6))*transl(a(6),0,0)*rotz(theta6)*transl(0,0,d(6))

X1=inv(A1*A2*A3)*T
vpa(X1,3) %Լ��������Ԫ��С��λ��Ϊ3
X2=A4*A5*A6
vpa(X2,2) %Լ��������Ԫ��С��λ��Ϊ3

%��theta4 theta5 theta6
theta5=-acos(X1(2,3))
%theta5=acos(X1(2,3))
theta4=acos(X1(3,3)/sin(theta5))
theta6=acos(X1(2,1)/sin(theta5))

theta(4)=theta4 
theta(5)=theta5
theta(6)=theta6

q= theta %���theta����
q=q/DR
drivebot(R,q)