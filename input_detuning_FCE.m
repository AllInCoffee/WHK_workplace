clear all

tao=18.5;
tao_theta=185;
sigma_FCD=7.2;
xi_T=0.074;
eta_lin=0.4;
eta_c=1;
gamma_TPA=0.11;
gamma_FCA=0.2;
chi = 0.55;
k_l=50 * abs(chi);
k = 150 * abs(chi);
k_T=k+k_l;

%tao_theta.*xi_T.*E.*(eta_lin.*eta_c + 2.*alpha_TPA.*E + 2.*gamma_FCA.*tao.*E.^2)).^2 ...

P_range=1:160:1000;
t=size(P_range,2);
for i=1:t
    P=P_range(i);
    h=@(delta,E) k*P-E.*( (delta+chi.*E - tao.*E.^2 - sigma_FCD.*tao.^0.8.*E.^1.6).^2 + ...
    (k_T./2+gamma_TPA.*E+gamma_FCA.*tao.*E.^2).^2 );
    fimplicit(h,[-700 700 0 5]);hold on
end

xlabel('\Delta')
ylabel('E')
title('Single Ring Resonator with Kerr and FCE')
