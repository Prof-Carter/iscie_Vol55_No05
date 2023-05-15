% sample_yalmip_sedumi.m
% last modified: 2023/05/15 by Masakatsu KAWATA

plant                                           % アクロボットに対してシステム行列などを定義した M ファイルの実行
% -------------------------------------------
gamma = sdpvar(1,1);                            % 決定変数 gamma（ガンマ）：スカラ
X     = sdpvar(n,n,'sy');                       % 決定変数 X：n×n の対称行列
Z     = sdpvar(n,m,'f');                        % 決定変数 Z：n×m の長方行列
% -------------------------------------------
ep = 1e-5;                                      % 十分小さな正数 <=== 加筆
% -------------------------------------------
LMI = [];                                       % LMI の記述の初期化
% -------------------------------------------
AX = A*X + B2*Z';
CX = C1*X + D12*Z';
% -------------------------------------------
M1 = [ r*X      AX-c*X
       AX'-c*X  r*X    ];
LMI = [LMI, M1 >= eps*eye(length(M1))];         % M1 ≧ eps*I (> 0)
% -------------------------------------------
M2 = [ AX+AX' B1            CX'
       B1'   -gamma*eye(q)  D11'
       CX     D11          -gamma*eye(p) ];
LMI = [LMI, M2 <= -eps*eye(length(M2))];        % M2 ≦ -eps*I (< 0)
% -------------------------------------------
opt = sdpsettings; opt.solver = 'sedumi';       % ソルバとして SeDuMi を利用
optimize(LMI,gamma,opt)                         % 目的関数を E = gamma とした凸最適化問題を解く
%%% もしくは上 2 行の代わりに
%%% optimize(LMI,gamma,sdpsettings('solver','sedumi')) 
% -------------------------------------------
gamma_opt = value(gamma)                        % 得られた gamma の最適解 gamma_opt
X_opt = value(X)                                % 得られた X の最適解 X_opt
Z_opt = value(Z)                                % 得られた Z の最適解 Z_opt
% -------------------------------------------
K_opt = Z_opt'*inv(X_opt)                       % コントローラゲイン K_opt
