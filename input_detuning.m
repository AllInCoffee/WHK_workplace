clear all
n_kerr=0.55;
tao=18.5;
tao_theta=185;
sigma_FCD=7.2;
xi_T=0.074;
eta_lin=0.4;
eta_c=1;
alpha_TPA=0.11;
gamma_FCA=0.2;

delta=-3;

f= @(x,y) x-y.*( (-delta+n_kerr.*y - (tao.*y.^2 + sigma_FCD.*tao.^0.8.*y.^1.6) + ...
    tao_theta.*xi_T.*y.*(eta_lin.*eta_c + 2.*alpha_TPA.*y + 2.*gamma_FCA.*tao.*y.^2)).^2 ...
    + (1+alpha_TPA.*y+gamma_FCA.*tao.*y.^2).^2 );
%%Method1 fimplicit
fimplicit(f,[0 80 0 1], 'bx') 
xlabel('P')
ylabel('E')

figure()
E_range=0:0.1:1;
s = size(E_range,2);
for i=1:s
    E=E_range(i);
    g=@(delta,P) P-E.*( (-delta+n_kerr.*E - (tao.*E.^2 + sigma_FCD.*tao.^0.8.*E.^1.6) + ...
    tao_theta.*xi_T.*E.*(eta_lin.*eta_c + 2.*alpha_TPA.*E + 2.*gamma_FCA.*tao.*E.^2)).^2 ...
    + (1+alpha_TPA.*E+gamma_FCA.*tao.*E.^2).^2 );
    fimplicit(g,[-30 30 0 100]);hold on
end
xlabel('\delta')
ylabel('P')

figure()
P_range=1:10:100;
t=size(P_range,2);
for i=1:t
    P=P_range(i);
    h=@(delta,E) P-E.*( (-delta+n_kerr.*E - (tao.*E.^2 + sigma_FCD.*tao.^0.8.*E.^1.6) + ...
    tao_theta.*xi_T.*E.*(eta_lin.*eta_c + 2.*alpha_TPA.*E + 2.*gamma_FCA.*tao.*E.^2)).^2 ...
    + (1+alpha_TPA.*E+gamma_FCA.*tao.*E.^2).^2 );
    fimplicit(h,[-40 40 0 1.5]);hold on
end
xlabel('\delta')
ylabel('E')

% example for multiple solutions
% syms Z
% eqn=0==4*Z^4-11*Z^2+7-3i;
% Root=double(vpasolve(eqn))


% syms E A
% %delta=-3;
% %P=10;
% delta1=-3:1:3;
% P1=0.1:20:100;
% for i=1:size(delta1,2)
%     delta=delta1(i);
%     for j = 1:size(P1,2)
%         P=P1(j);
% 
%         % eqn=0+0i==A.*(1i.*delta - 1i.*(n_kerr.*abs(A)^4 - (tao.*abs(A).^4+sigma_FCD.*(tao.*abs(A).^4).^0.8) ...
%         %     + tao_theta.*xi_T.*abs(A).^2.*(eta_lin.*eta_c + 2.*alpha_TPA.*abs(A).^2 + 2.*gamma_FCA.*tao.*abs(A).^4)) ...
%         %     - (1+alpha_TPA.*abs(A).^2 + gamma_FCA.*n_kerr.*abs(A)^4)) + sqrt(P);
%         eqn1=0== P-E.*( (-delta+n_kerr.*E - (tao.*E.^2 + sigma_FCD.*tao.^0.8.*E.^1.6) + ...
%              tao_theta.*xi_T.*E.*(eta_lin.*eta_c + 2.*alpha_TPA.*E + 2.*gamma_FCA.*tao.*E.^2)).^2 ...
%              + (1+alpha_TPA.*E+gamma_FCA.*tao.*E.^2).^2 );
%         root=double(vpasolve(eqn1,E));
%         % xlabel('$A$', 'Interpreter','latex')
%         % legend()
%         syms re_A im_A A
%         k=@(re_A,im_A) root - im_A.^2+re_A.^2;
%         %fimplicit(k,[-5 5 -5 5])
% 
%         temp=sqrt(root);
% 
%         re_A_range=linspace(-temp, temp,10);
%         for i=1:size(re_A_range,2)
%             re_A=re_A_range(i);
%             m11= - 1i.*delta - 1i.*2.*n_kerr.*abs(A).^2+ 1i.*(tao.*abs(A).^4+sigma_FCD.*(tao.*abs(A).^4).^0.8) ...
%                 -1i.*tao_theta.*xi_T.*abs(A).^2.*(eta_lin.*eta_c + 2.*alpha_TPA.*abs(A).^2 + 2.*gamma_FCA.*tao.*abs(A).^4)...
%                 -1-gamma_FCA.*tao.*abs(A).^4-2.*alpha_TPA.*abs(A).^2;
%             im_A=sqrt(root-re_A.^2);
% 
%             M11=subs(m11,A,re_A + im_A*1i);
% 
%             %fplot(real(M11),[-temp temp],'x')  %% exceeding this range, im_A will not be real
% 
% 
%             m12=  -1i.*n_kerr.*A.^2 - alpha_TPA.*A.^2;
%             m13=  1i.*A.*(1+0.8.*sigma_FCD.*(tao.*abs(A).^4).^-0.2) - gamma_FCA.*A;
%             m14=  -1i.*A;
% 
%             m21=  conj(m12);
%             m22=  conj(m11);
%             m23=  conj(m13);
%             m24=  conj(m14);
% 
%             m31=  2.*abs(A).^2.*conj(A);
%             m32=  conj(m31);
%             m33=  -1./tao;
%             m34=0;
% 
%             m41=  xi_T.*(eta_lin.*eta_c + 4.*alpha_TPA.*abs(A).^2 + 2.*gamma_FCA.*tao.*abs(A).^4).*conj(A);
%             m42=  conj(m41);
%             m43=  2.*xi_T.*gamma_FCA.*abs(A).^2;
%             m44=  -1./tao_theta;
% 
%             M12=subs(m12,A,re_A + im_A*1i);
%             M13=subs(m13,A,re_A + im_A*1i);
%             M14=subs(m14,A,re_A + im_A*1i);
% 
%             M21=subs(m21,A,re_A + im_A*1i);
%             M22=subs(m22,A,re_A + im_A*1i);
%             M23=subs(m23,A,re_A + im_A*1i);
%             M24=subs(m24,A,re_A + im_A*1i);
% 
%             M31=subs(m31,A,re_A + im_A*1i);
%             M32=subs(m32,A,re_A + im_A*1i);
%             M33=subs(m33,A,re_A + im_A*1i);
%             M34=subs(m34,A,re_A + im_A*1i);
% 
%             M41=subs(m41,A,re_A + im_A*1i);
%             M42=subs(m42,A,re_A + im_A*1i);
%             M43=subs(m43,A,re_A + im_A*1i);
%             M44=subs(m44,A,re_A + im_A*1i);
%             % 
%             M= double([M11 M12 M13 M14;
%                 M21 M22 M23 M24;
%                 M31 M32 M33 M34;
%                 M41 M42 M43 M44]);
%             % 
%             lambda=eig(M);  % method 2: B = sort(A,'ComparisonMethod','real')
%             sel= real(lambda);
%             res=sort(sel,'descend'); 
%             if res(1)>0
%                 if res(2)<0
%                     %plot(delta,P,'r');hold on % BI region
%                     semilogy(delta,P,'ro');hold on
%                 else
%                     %plot(delta,P,'b');hold on % SP
%                     semilogy(delta,P,'bo');hold on
% 
%                 end
%                 % to set according to P and E range
%                 ylim([0,100]);
%                 xlim([-3,3]);
%             else
%                 %disp('stable')
% 
%             end
% 
% 
%             %tf = isequal(lambda(1), lambda(2),lambda(3),lambda(4))
% 
% 
%         end
%     end
% end
% 