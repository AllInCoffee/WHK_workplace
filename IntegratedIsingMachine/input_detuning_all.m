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


figure()
P_range=1:160:1000;
t=size(P_range,2);
for i=1:t
    P=P_range(i);
    h=@(delta,E) P-E.*( (delta+chi.*E).^2 + 1 );
    fimplicit(h,[-200 200 0 30]);hold on
end

xlabel('\Delta')
ylabel('E')
title('Single Ring Resonator with Only Kerr ')

figure()
P_range=1:16000:100000;
t=size(P_range,2);
for i=1:t
    P=P_range(i);
    h=@(delta,E) P-E.*( (delta+chi.*E).^2 + 1 );
    fimplicit(h,[-200 200 0 300]);hold on
end

xlabel('\Delta')
ylabel('E')
title('Single Ring Resonator with Only Kerr ')


figure()
P_range=1:160:1000;
t=size(P_range,2);
for i=1:t
    P=P_range(i);
    h=@(delta,E) P-E.*( (delta+chi.*E - tao.*E.^2 - sigma_FCD.*tao.^0.8.*E.^1.6).^2 + ...
    (1+gamma_TPA.*E+gamma_FCA.*tao.*E.^2).^2 );
    fimplicit(h,[-200 200 0 3]);hold on
end

xlabel('\Delta')
ylabel('E')
title('Single Ring Resonator with Kerr and FCE')


figure()
P_range=1:160:1000;
t=size(P_range,2);
for i=1:t
    P=P_range(i);
    h=@(delta,E) P-E.*( (delta+chi.*E - tao.*E.^2 - sigma_FCD.*tao.^0.8.*E.^1.6...
        +tao_theta.*xi_T.*E.*(eta_lin.*eta_c + 2.*gamma_TPA.*E + 2.*gamma_FCA.*tao.*E.^2)).^2 + ...
    (1+gamma_TPA.*E+gamma_FCA.*tao.*E.^2).^2 );
    fimplicit(h,[-200 200 0 3]);hold on
end

xlabel('\Delta')
ylabel('E')
title('Single Ring Resonator with Kerr, FCE and TOE')