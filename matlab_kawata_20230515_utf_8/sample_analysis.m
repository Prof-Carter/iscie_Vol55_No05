% sample_analysis.m
% last modified: 2023/05/15 by Masakatsu KAWATA

Acl = A  + B2*K_opt;
Ccl = C1 + D12*K_opt;

% /////////////////////////////////////////////////////////////////////////
% 極配置仕様の確認：Acl = A + B2*K の固有値
% /////////////////////////////////////////////////////////////////////////
% ---------------------------------------------------
% 円領域の描画
theta = 0:0.01:2*pi;
x_circle = r*cos(theta) + c;
y_circle = r*sin(theta);
figure(1);
plot(x_circle,y_circle,'LineWidth',1.5)
hold on
% ---------------------------------------------------
% A + B2*K の固有値の描画
eigen = eig(Acl)
figure(1)
plot(real(eigen),imag(eigen),'x','LineWidth',1.5,'MarkerSize',10)
hold off
grid on
axis('square')
xlim([-15 5])
ylim([-10 10])
set(gca,'FontSize',16,'FontName','Arial')
xlabel('${\rm{Re}}[\lambda]$','FontSize',18,'Interpreter','latex')
ylabel('${\rm{Im}}[\lambda]$','FontSize',18,'Interpreter','latex')

% /////////////////////////////////////////////////////////////////////////
% H-infinity 制御仕様の確認
% /////////////////////////////////////////////////////////////////////////
% ---------------------------------------------------
% Gzw(jw) の特異値の計算
w   = logspace(-2,2,500);
sys = ss(Acl,B1,Ccl,D11);
sv  = sigma(sys,w);
% ---------------------------------------------------
% 特異値の描画
figure(2)
semilogx(w,sv,[1e-2 1e2],gamma_opt*[1 1],'--','LineWidth',1.5)
grid on
set(gca,'FontSize',16,'FontName','Arial')
xlabel('$\omega$ [rad/s]','FontSize',18,'Interpreter','latex')
ylabel('Magnitude','FontSize',18,'Interpreter','latex')
set(gca,'XTick',[1e-2 1e-1 1e0 1e1 1e2])
legend({'$\max\sigma[{G}_{zw}(j\omega)]$','$\gamma$'},...
       'Location','southwest')
set(legend,'FontSize',18,'Interpreter','latex')

% ---------------------------------------------------
sys1 = ss(Acl,B1,-[Cp 0],1);    %%% sys1: w --> e
sys2 = ss(Acl,B1,-[Cp 0],0);    %%% sys2: w --> theta1
sv1 = sigma(sys1,w);            %%% sys1 の特異値の計算
sv2 = sigma(sys2,w);            %%% sys2 の特異値の計算
figure(3)
semilogx(w,20*log10(sv1),w,20*log10(sv2),'LineWidth',1.5)
hold on

sys3 = gamma_opt*tf([1 0],bs);          %%% sys3: gamma/Ws(s)
sys4 = gamma_opt*tf(1,[bt2 bt1 bt0]);   %%% sys4: gamma/Wt(s)
sv3 = sigma(sys3,w);                    %%% sys3 の特異値の計算
sv4 = sigma(sys4,w);                    %%% sys4 の特異値の計算
figure(3)
semilogx(w,20*log10(sv3),'--',w,20*log10(sv4),'--','LineWidth',1.5)
hold off

grid on
ylim([-100 20])
set(gca,'FontSize',16,'FontName','Arial')
xlabel('$\omega$ [rad/s]','FontSize',18,'Interpreter','latex')
ylabel('Gain [dB]','FontSize',18,'Interpreter','latex')
set(gca,'XTick',[1e-2 1e-1 1e0 1e1 1e2])
legend({'$w \rightarrow e$',...
        '$w \rightarrow \theta_{1}$',...
        '$\gamma/{W}_{\rm{s}}(s)$','$\gamma/{W}_{\rm{t}}(s)$'},...
       'Location','southwest')
set(legend,'FontSize',18,'Interpreter','latex')

% /////////////////////////////////////////////////////////////////////////
figure(1); movegui('northwest')
figure(2); movegui('north')
figure(3); movegui('northeast')

figure(1); exportgraphics(gcf,'eigen.pdf','ContentType','vector')
figure(2); exportgraphics(gcf,'gain1.pdf','ContentType','vector')
figure(3); exportgraphics(gcf,'gain2.pdf','ContentType','vector')