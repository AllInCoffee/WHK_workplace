clear all
function f=ODEsystem(t,y,m,u)
f = zeros(2, 1);
%n_kerr=0.55;
tao=18.5;
sigma_FCD=7.2;
gamma_TPA=0.11;
gamma_FCA=0.2;
%%%

x = 0.55;
k=150;
k_l=50;
K = k * abs(x);
K_T=(k+k_l)*abs(x);
Delta_th= -sqrt(3)/2.0 * K_T;

%%%
f(1)=(-1i.*m.*Delta_th-1i.*(x.*abs(y(1)).^2 - (y(2)+sigma_FCD.*y(2)^0.8) )-(K_T./2+gamma_TPA.*abs(y(1)).^2+gamma_FCA.*y(2))).*y(1) -sqrt(K)*u;
f(2)=abs(y(1)).^4 - y(2)./tao;
end

tSpan1=[0:1:50];
tSpan2=[0:0.01:2];
u=20;
m=1.3;
figure(1)
[Time,Y]=ode45(@(t,y) ODEsystem(t,y,m,u), tSpan1, [1,3]); %  
plot(Time,Y(:,2),'linewidth',1.5);
hold on
plot(Time,real(Y(:,1)),'linewidth',1.5);
hold on
plot(Time,imag(Y(:,1)),'linewidth',1.5);
hold on
legend('n','a\_real','a\_imag')

m_range = 0.7:0.3:2;
figure(2)
for i = 1:length(m_range)
    m=m_range(i);
    txt1 = ['n , m = ',num2str(m)];
    [Time,Y]=ode45(@(t,y) ODEsystem(t,y,m,u), tSpan1, [-0.2628 - 0.8714i,13]); % 2 stays stable, 
    plot(Time,Y(:,2),'linewidth',1.5,'DisplayName',txt1);
    hold on
    legend()  
end
title('Development of n Over Time')
figure(3)
for j = 1:length(m_range)
    m=m_range(j);
    txt2 = ['|a| , m = ',num2str(m)];
    [Time,Y]=ode45(@(t,y) ODEsystem(t,y,m,u), tSpan2, [-0.2628 - 0.8714i,13]); % 2 stays stable, 
    plot(Time,abs(Y(:,1)),'linewidth',1.5,'DisplayName',txt2);
    hold on
    legend()  
end
title('Development of |a| Over Time')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % find the stationary points 
% options = optimoptions('fsolve','Display','iter');
% y0=[10,3]; %randomly chosen
% tao=18.5;
% sigma_FCD=7.2;
% gamma_TPA=0.11;
% gamma_FCA=0.2;
% %%%
% 
% x = 0.55;
% k=150;
% k_l=50;
% K = k * abs(x);
% K_T=(k+k_l)*abs(x);
% %Delta_th= -sqrt(3)/2.0 * K_T;
% 
% u=20;%%
% m=1.3;%% old 1.3 Delta_th
% Delta_th= -sqrt(3)/2.0 * K_T;
% m.*Delta_th
% F=@(y) [(-1i.*m.*Delta_th-1i.*(x.*abs(y(1)).^2-(y(2)+sigma_FCD.*y(2).^0.8))-(K_T./2+gamma_TPA.*abs(y(1)).^2+gamma_FCA.*y(2))).*y(1)-sqrt(K).*u;abs(y(1)).^4 - y(2)./tao];
% 
% [Ys,fval]=fsolve(F, y0,options)
% % calculate N0
% N0=tao.*abs(Ys(1)).^4
% % 
% 
% % stationary points:
% %N0=12.6979 (calculated through tao.*E.^2), -0.2628 - 0.8714i (fsolve)
% 
% %% if we insert the stationary points as y0, the value stay near the stationary points.
% %% which proves that the results are correct. 

%%%%%%%%%%%%%%%%%
% calculate Delta_th_new
% from dz/dt=0 we get : 
sigma_FCD=7.2;
gamma_TPA=0.11;
gamma_FCA=0.2;
x = 0.55;
k=150;
k_l=50;
K = k * abs(x);
K_T=(k+k_l)*abs(x);
syms E Delta n
P=((K_T./2 + gamma_TPA.*E + gamma_FCA.* n).^2 + (Delta + x.*E-n-sigma_FCD.*n.^0.8).^2).*E./K
DP=diff(P,E)

%(2*((11*E)/100 + n/5 + 55)^2)/165 + (2*E*((11*Delta)/10 + (1573*E)/2500 - (132*n)/125 - (198*n^(4/5))/25 + 121/10))/165 + (2*(Delta + (11*E)/20 - n - (36*n^(4/5))/5)^2)/165
 





