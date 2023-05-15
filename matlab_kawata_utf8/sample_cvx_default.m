% sample_cvx_default.m
% last modified: 2023/05/15 by Masakatsu KAWATA

plant                                           % アクロボットに対してシステム行列などを定義した M ファイルの実行
% -------------------------------------------
cvx_begin sdp                                   % SDP（半正定値計画問題；LMI を制約とした凸最適化問題）の開始
    variable gamma(1,1)                         % 決定変数 gamma（ガンマ）：スカラ
    variable X(n,n) symmetric                   % 決定変数 X：n×n の対称行列
    variable Z(n,m)                             % 決定変数 Z：n×m の長方行列
    % ---------------------------------------
    minimize(gamma)                             % 目的関数 E = gamma を最小化することを宣言
    % ---------------------------------------
    AX = A*X + B2*Z'; 
    CX = C1*X + D12*Z'; 
    % ---------------------------------------
    M1 = [ r*X      AX-c*X 
           AX'-c*X  r*X    ];
    M1 > 0;                                     % M1 > 0
    % ---------------------------------------
    M2 = [ AX+AX' B1            CX'
           B1'   -gamma*eye(q)  D11'
           CX     D11          -gamma*eye(p) ];
    M2 < 0;                                     % M2 < 0
cvx_end                                         % SDP の終了
% -------------------------------------------
K = Z'*inv(X)                                   % コントローラゲイン K
