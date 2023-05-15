% sample_yalmip_sedumi.m
% last modified: 2023/05/15 by Masakatsu KAWATA

plant                                           % �A�N���{�b�g�ɑ΂��ăV�X�e���s��Ȃǂ��`���� M �t�@�C���̎��s
% -------------------------------------------
gamma = sdpvar(1,1);                            % ����ϐ� gamma�i�K���}�j�F�X�J��
X     = sdpvar(n,n,'sy');                       % ����ϐ� X�Fn�~n �̑Ώ̍s��
Z     = sdpvar(n,m,'f');                        % ����ϐ� Z�Fn�~m �̒����s��
% -------------------------------------------
ep = 1e-5;                                      % �\�������Ȑ��� <=== ���M
% -------------------------------------------
LMI = [];                                       % LMI �̋L�q�̏�����
% -------------------------------------------
AX = A*X + B2*Z';
CX = C1*X + D12*Z';
% -------------------------------------------
M1 = [ r*X      AX-c*X
       AX'-c*X  r*X    ];
LMI = [LMI, M1 >= eps*eye(length(M1))];         % M1 �� eps*I (> 0)
% -------------------------------------------
M2 = [ AX+AX' B1            CX'
       B1'   -gamma*eye(q)  D11'
       CX     D11          -gamma*eye(p) ];
LMI = [LMI, M2 <= -eps*eye(length(M2))];        % M2 �� -eps*I (< 0)
% -------------------------------------------
opt = sdpsettings; opt.solver = 'sedumi';       % �\���o�Ƃ��� SeDuMi �𗘗p
optimize(LMI,gamma,opt)                         % �ړI�֐��� E = gamma �Ƃ����ʍœK����������
%%% �������͏� 2 �s�̑����
%%% optimize(LMI,gamma,sdpsettings('solver','sedumi')) 
% -------------------------------------------
gamma_opt = value(gamma)                        % ����ꂽ gamma �̍œK�� gamma_opt
X_opt = value(X)                                % ����ꂽ X �̍œK�� X_opt
Z_opt = value(Z)                                % ����ꂽ Z �̍œK�� Z_opt
% -------------------------------------------
K_opt = Z_opt'*inv(X_opt)                       % �R���g���[���Q�C�� K_opt
