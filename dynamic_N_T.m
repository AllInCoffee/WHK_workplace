clear all
function dy = fun(t,y)
n_kerr=0.55;
tao=18.5;
tao_theta=185;
sigma_FCD=7.2;
xi_T=0.074;
eta_lin=0.4;
eta_c=1;
alpha_TPA=0.11;
gamma_FCA=0.2;
%%%

P=0.5;
delta=-3;
%%%
dy=[(1i.*delta-1i.*(n_kerr.*abs(y(1)).^2 - (y(2)+sigma_FCD.*y(2)^0.8) + y(3))...
    -(1+alpha_TPA.*abs(y(1)).^2+gamma_FCA.*y(2))).*y(1) + sqrt(P);
    abs(y(1)).^4 - y(2)./tao;
    xi_T.*abs(y(1)).^2.*(eta_c.*eta_lin + 2.*alpha_TPA.*abs(y(1))^2 + 2.*gamma_FCA.*y(2)) - y(3)./tao_theta]  ;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % % find the stationary points 
% options = optimoptions('fsolve','Display','iter');
% y0=[0.5,1,1]; %randomly chosen
% n_kerr=0.55;
% tao=18.5;
% tao_theta=185;
% sigma_FCD=7.2;
% xi_T=0.074;
% eta_lin=0.4;
% eta_c=1;
% alpha_TPA=0.11;
% gamma_FCA=0.2;
% %%%
% 
% P=0.5;
% delta=-3;
% F=@(y) [(1i.*delta-1i.*(n_kerr.*abs(y(1)).^2 - (y(2)+sigma_FCD.*y(2)^0.8) + y(3))...
%     -(1+alpha_TPA.*abs(y(1)).^2+gamma_FCA.*y(2))).*y(1) + sqrt(P);
%     abs(y(1)).^4 - y(2)./tao;
%     xi_T.*abs(y(1)).^2.*(eta_c.*eta_lin + 2.*alpha_TPA.*abs(y(1))^2 + 2.*gamma_FCA.*y(2)) - y(3)./tao_theta]  ;
% [Ys,fval]=fsolve(F, y0,options)
% % % calculate N0
% % N0=tao.*abs(Ys(1)).^4
% % % % calculate T0
% % T0=tao_theta.*xi_T.*abs(Ys(:,1))^2.*(eta_lin.*eta_c+2.*alpha_TPA.*abs(Ys(:,1))^2+2.*gamma_FCA.*tao.*abs(Ys(:,1))^4)
% % abs(Ys(2))
% % abs(Ys(3))

% stationary points: A=0.4089 + 0.0359i 
%N0=0.5254 (calculated through tao.*E.^2), 0.5083 + 0.1135i (fsolve)
%T0=1.4933 (calculated), 1.4492+0.2647i(fsolve)
%% if we insert the stationary points as y0, the value stay near the stationary points.
%% which proves that the results are correct. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% plot the development of a_real and a_imag, n, T
tSpan=[0 50];
[Time,Y]=ode45(@(t,y) fun(t,y), tSpan, [0.4+0.04i,0.5,1.5]); % 2 stays stable, 
figure(1)
plot(Time,real(Y(:,1)),'linewidth',1.5);
hold on
plot(Time,imag(Y(:,1)),'linewidth',1.5);
hold on
plot(Time,abs(Y(:,1)),'linewidth',1.5);
hold on
plot(Time,Y(:,2),'linewidth',1.5);
hold on
plot(Time,Y(:,3),'linewidth',1.5)
xlabel('t')
legend('a\_real','a\_imag','|a|','n','T','Location','best')
figure(2)
plot(real(Y(:,1)),imag(Y(:,1)))
xlabel('a\_real')
ylabel('a\_imag')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% with changing N, we only need the first two equations of the coupled systems
%% in paper Multibistability and self-pulsation in nonlinear high- Q silicon microring resonators considering
%% thermo-optical effect (2013)

function dy = funN(t,y)
n_kerr=0.55;
tao=18.5;
tao_theta=185;
sigma_FCD=7.2;
xi_T=0.074;
eta_lin=0.4;
eta_c=1;
alpha_TPA=0.11;
gamma_FCA=0.2;
%%%

P=0.5;
delta=-3;
T=2; % fix the T 
%%%
dy=[(1i.*delta-1i.*(n_kerr.*abs(y(1)).^2 - (y(2)+sigma_FCD.*y(2)^0.8) + T)...
    -(1+alpha_TPA.*abs(y(1)).^2+gamma_FCA.*y(2))).*y(1) + sqrt(P);
    abs(y(1)).^4 - y(2)./tao]  ;
end


figure(3)
for n=0:1:2
    for a_real=-2:1:2
        for a_imag=-2:1:2
            [ts,ys]= ode45(@(t,y) funN(t,y), tSpan, [a_real+1i.*a_imag,n]);
            %ys(:,1)
            % plot(ys(:,2),real(ys(:,1)),'b');hold on
            % plot(ys(:,2),imag(ys(:,1)),'r');hold on
            %plot(ys(:,2),abs(ys(:,1)),'k');hold on
            plot3(real(ys(:,1)),imag(ys(:,1)),ys(:,2));hold on
        end
    end
end
xlabel('n')
ylabel('a')
%legend('a\_real','a\_imag','n','Location','best')
xlabel('a\_real')
ylabel('a\_imag')
zlabel('n')

%% a versus n with direction:
n_kerr=0.55;
tao=18.5;
tao_theta=185;
sigma_FCD=7.2;
xi_T=0.074;
eta_lin=0.4;
eta_c=1;
alpha_TPA=0.11;
gamma_FCA=0.2;
%%%

P=0.5;
delta=-3;
T=2; % fix the T
[X1,X2,Y] = meshgrid(-2:0.5:2,-2:0.5:2,0:0.5:2);
diff_a=(1i.*delta-1i.*(n_kerr.*(X1.^2+X2.^2) - (Y+sigma_FCD.*Y.^0.8) + T)...
    -(1+alpha_TPA.*(X1.^2+X2.^2)+gamma_FCA.*Y)).*(X1+1i.*X2) + sqrt(P);
[DX1,DX2,DY]=gradient(diff_a,0.5,0.5,0.5);

quiver3(X1,X2,Y,DX1,DX2,DY, 'k')
% 
% %%%%%%%%%%%%%%%%%%%%%%%%% with fixed n and chaning T %%%%%%%%%%%%%%%%%%%%%%%%
function dy = funT(t,y)
n_kerr=0.55;
tao=18.5;
tao_theta=185;
sigma_FCD=7.2;
xi_T=0.074;
eta_lin=0.4;
eta_c=1;
alpha_TPA=0.11;
gamma_FCA=0.2;
%%%

P=0.5;
delta=-3;
n=0.52; % fix the N 
%%%
dy=[(1i.*delta-1i.*(n_kerr.*abs(y(1)).^2 - (n+sigma_FCD.*n^0.8) + y(2))...
    -(1+alpha_TPA.*abs(y(1)).^2+gamma_FCA.*n)).*y(1) + sqrt(P);
    xi_T.*abs(y(1)).^2.*(eta_c.*eta_lin + 2.*alpha_TPA.*abs(y(1))^2 + 2.*gamma_FCA.*n) - y(2)./tao_theta]  ;
end


figure(4)
for T=0:0.5:3
    for a_real=-2:0.5:2
        for a_imag=-2:0.5:2
            [ts,ys]= ode45(@(t,y) funT(t,y), tSpan, [a_real+1i.*a_imag,T]);
            %ys(:,1)
            % plot(ys(:,2),real(ys(:,1)),'b');hold on
            % plot(ys(:,2),imag(ys(:,1)),'r');hold on
            %plot(ys(:,2),abs(ys(:,1)),'k');hold on
            plot3(real(ys(:,1)),imag(ys(:,1)),ys(:,2));hold on
        end
    end
end

%legend('a\_real','a\_imag','n','Location','best')
xlabel('a\_real')
ylabel('a\_imag')
zlabel('T')

%% a versus T with direction:
n_kerr=0.55;
tao=18.5;
tao_theta=185;
sigma_FCD=7.2;
xi_T=0.074;
eta_lin=0.4;
eta_c=1;
alpha_TPA=0.11;
gamma_FCA=0.2;
%%%

P=0.5;
delta=-3;
n=0.52; % fix the n
[X1,X2,Z] = meshgrid(-2:0.5:2,-2:0.5:2,0:0.5:3);
diff_a=(1i.*delta-1i.*(n_kerr.*(X1.^2+X2.^2) - (n+sigma_FCD.*n.^0.8) + Z)...
    -(1+alpha_TPA.*(X1.^2+X2.^2)+gamma_FCA.*n)).*(X1+1i.*X2) + sqrt(P);
[DX1,DX2,DZ]=gradient(diff_a,0.5,0.5,0.5);

quiver3(X1,X2,Z,DX1,DX2,DZ, 'k')