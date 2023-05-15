% plant.m
% last modified: 2023/05/15 by Masakatsu KAWATA

clear
format compact

% /////////////////////////////////////////////////////////////////////////
% アクロボットのパラメータ
% /////////////////////////////////////////////////////////////////////////
g   = 9.81;         %% 重力加速度
m1  = 7.02e-001;    %% リンク 1 の質量
m2  = 1.40e-001;    %% リンク 2 の質量
% --------------------------------------
L1  = 1.23e-001;    %% リンク 1 の全長
L2  = 2.35e-001;    %% リンク 2 の全長
l1  = 1.25e-002;    %% リンク 1 の軸から重心までの長さ
l2  = 1.60e-001;    %% リンク 2 の軸から重心までの長さ
% --------------------------------------
Jg1 = 5.98e-003;    %% リンク 1 の重心まわりの慣性モーメント
c1  = 4.79e-003;    %% リンク 1 の粘性摩擦係数
Jg2 = 1.09e-003;    %% リンク 2 の重心まわりの慣性モーメント
% --------------------------------------
a2s = 1.46e+001;    %% リンクおよびモータ，モータドライバの特性により定まる定数
c2s = 9.46e+001;
km  = 2.90e+002;
% --------------------------------------
a1 = Jg1 + m1*l1^2 + m2*L1^2;
a2 = Jg2 + m2*l2^2;
a3 = m2*L1*l2;
a4 = m1*g*l1 + m2*g*L1;
a5 = m2*g*l2;

% /////////////////////////////////////////////////////////////////////////
% 重み
% /////////////////////////////////////////////////////////////////////////
% Wt(p) = bt2*s^2 + bt1*s + bt0
bt0 = 1/10;
bt1 = 2/10;
bt2 = 1/10;

% Ws(p) = bs/s
bs  = 1;

% /////////////////////////////////////////////////////////////////////////
% 一般可制御対象 S の係数行列の定義
% /////////////////////////////////////////////////////////////////////////
Ep = [ 1     0      0         0
       0     1      0         0
       0     0  a1+a2-2*a3  a2-a3
       0     0    a2-a3      a2s  ];
Ap = [ 0     0      1         0
       0     0      0         1
     a4-a5  -a5    -c1        0
      -a5   -a5     0       -c2s ];
Bp = [ 0
       0
       0
       km ];
% --------------------------------------
Cp1 = [ 1  0 ];
Cp  = [ Cp1  zeros(1,2) ];
% --------------------------------------
Ct1  = [ bt0  0 ];
Ct2  = [ bt1  0 ];
Ct   = [ Ct1  Ct2 ];
% --------------------------------------
Ctd2 = [ bt2  0 ];
Ctd  = [ zeros(1,2)  Ctd2 ];
% --------------------------------------

Ed   = [ Ep          zeros(4,1)
         zeros(1,4)  1          ]; 
Ad   = [ Ap          zeros(4,1)
        -Cp          zeros(1,1) ];
B1d  = [ zeros(4,1)
         1          ];
B2d  = [ Bp
          0  ];
C11d = [ Ct          0
         zeros(1,4)  bs ]; 
C12d = [ Ctd         0
         zeros(1,4)  0 ]; 

% --------------------------------------
A   = inv(Ed)*Ad;
B1  = inv(Ed)*B1d;
B2  = inv(Ed)*B2d;

C1  = C11d + C12d*A;
D11 =        C12d*B1;
D12 =        C12d*B2;

% /////////////////////////////////////////////////////////////////////////
% 次数の定義
% /////////////////////////////////////////////////////////////////////////
n = 5;      % x の次数
q = 1;      % w の次数
m = 1;      % u の次数
p = 2;      % z の次数

% /////////////////////////////////////////////////////////////////////////
% 円領域の指定：中心 (c,0), 半径 r
% /////////////////////////////////////////////////////////////////////////
c = -6;
r =  8;
