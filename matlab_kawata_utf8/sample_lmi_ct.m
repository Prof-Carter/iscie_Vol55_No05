% sample_lmi_ct.m
% last modified: 2023/05/15 by Masakatsu KAWATA

plant                                           % アクロボットに対してシステム行列などを定義した M ファイルの実行
% -------------------------------------------
setlmis([])                                     % LMI の記述の初期化
gamma = lmivar(1, [1 0]);                       % 決定変数 gamma（ガンマ）：スカラ ===> gammma を X, Z より前に定義
X     = lmivar(1, [n 1]);                       % 決定変数 X：n×n の対称行列
Z     = lmivar(2, [n m]);                       % 決定変数 Z：n×m の長方行列
% -------------------------------------------
lmi1 = newlmi;                                  % 現在記述されている LMI に識別子 lmi1 の LMI を追加（-lmi1：正定）
lmiterm([-lmi1 1 1  X     ],  r,   1);          % (1,1) 要素に  r*X *1 を追加
lmiterm([-lmi1 1 2  X     ],  A,   1);          % (1,2) 要素に  A*X *1 を追加
lmiterm([-lmi1 1 2 -Z     ],  B2,  1);          % (1,2) 要素に B2*Z'*1 を追加（lmiterm の -Z は Z' を意味する）
lmiterm([-lmi1 1 2  X     ], -c,   1);          % (1,2) 要素に -c*X *1 を追加
lmiterm([-lmi1 2 2  X     ],  r,   1);          % (2,2) 要素に  r*X *1 を追加
% -------------------------------------------
lmi2 = newlmi;                                  % 現在記述されている LMI に識別子 lmi2 の LMI を追加（lmi2：負定）
lmiterm([ lmi2 1 1  X     ],  A,   1, 's');     % (1,1) 要素に   A*X *1 + ( A*X *1)' を追加
lmiterm([ lmi2 1 1 -Z     ],  B2,  1, 's');     % (1,1) 要素に  B2*Z'*1 + (B2*Z'*1)' を追加
lmiterm([ lmi2 1 2  0     ],  B1);              % (1,2) 要素に  B1 を追加・・・定数行例のみ
lmiterm([ lmi2 3 1  X     ],  C1,  1);          % (3,1) 要素に  C1*X *1 を追加
lmiterm([ lmi2 3 1 -Z     ],  D12, 1);          % (3,1) 要素に D12*Z'*1 を追加
lmiterm([ lmi2 3 2  0     ],  D11);             % (3,2) 要素に D11 を追加・・・定数行例のみ
lmiterm([ lmi2 2 2  gamma ], -1,   1);          % (2,2) 要素に -1*gamma*1 を追加
lmiterm([ lmi2 3 3  gamma ], -1,   1);          % (3,3) 要素に  1*gamma*1 を追加
% -------------------------------------------
lmisys = getlmis;                               % 定義された LMI を取得
% -------------------------------------------
cobj = zeros(1,decnbr(lmisys));                 % 決定変数の要素数を取得し，その長さの零ベクトル cobj を生成
cobj(1) = 1;                                    % cobj の 1 番目の要素を 1 に変更（決定変数の 1 番目の要素は gamma）
% -------------------------------------------
[cost,xopt] = mincx(lmisys,cobj);               % 目的関数を E = gamma とした凸最適化問題を解く
% -------------------------------------------
gamma_opt  = dec2mat(lmisys,xopt,gamma)         % 得られた gamma の最適解 gamma_opt
X_opt = dec2mat(lmisys,xopt,X)                  % 得られた X の最適解 X_opt
Z_opt = dec2mat(lmisys,xopt,Z)                  % 得られた Z の最適解 Z_opt
% -------------------------------------------
K_opt = Z_opt'*inv(X_opt)                       % コントローラゲイン K_opt