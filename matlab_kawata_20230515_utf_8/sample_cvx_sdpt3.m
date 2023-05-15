% sample_cvx_sdpt3.m
% last modified: 2023/05/15 by Masakatsu KAWATA

plant                                           % �A�N���{�b�g�ɑ΂��ăV�X�e���s��Ȃǂ��`���� M �t�@�C���̎��s
% -------------------------------------------
cvx_solver('sdpt3');                            % �\���o�Ƃ��� SDPT3 ���g�p
% -------------------------------------------
cvx_begin sdp                                   % SDP�i������l�v����GLMI �𐧖�Ƃ����ʍœK�����j�̊J�n
    variable gamma(1,1)                         % ����ϐ� gamma�i�K���}�j�F�X�J��
    variable X(n,n) symmetric                   % ����ϐ� X�Fn�~n �̑Ώ̍s��
    variable Z(n,m)                             % ����ϐ� Z�Fn�~m �̒����s��
    % ---------------------------------------
    minimize(gamma)                             % �ړI�֐� E = gamma ���ŏ������邱�Ƃ�錾
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
cvx_end                                         % SDP �̏I��
% -------------------------------------------
K = Z'*inv(X)                                   % �R���g���[���Q�C�� K
