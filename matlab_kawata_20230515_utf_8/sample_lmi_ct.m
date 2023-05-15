% sample_lmi_ct.m
% last modified: 2023/05/15 by Masakatsu KAWATA

plant                                           % �A�N���{�b�g�ɑ΂��ăV�X�e���s��Ȃǂ��`���� M �t�@�C���̎��s
% -------------------------------------------
setlmis([])                                     % LMI �̋L�q�̏�����
gamma = lmivar(1, [1 0]);                       % ����ϐ� gamma�i�K���}�j�F�X�J�� ===> gammma �� X, Z ���O�ɒ�`
X     = lmivar(1, [n 1]);                       % ����ϐ� X�Fn�~n �̑Ώ̍s��
Z     = lmivar(2, [n m]);                       % ����ϐ� Z�Fn�~m �̒����s��
% -------------------------------------------
lmi1 = newlmi;                                  % ���݋L�q����Ă��� LMI �Ɏ��ʎq lmi1 �� LMI ��ǉ��i-lmi1�F����j
lmiterm([-lmi1 1 1  X     ],  r,   1);          % (1,1) �v�f��  r*X *1 ��ǉ�
lmiterm([-lmi1 1 2  X     ],  A,   1);          % (1,2) �v�f��  A*X *1 ��ǉ�
lmiterm([-lmi1 1 2 -Z     ],  B2,  1);          % (1,2) �v�f�� B2*Z'*1 ��ǉ��ilmiterm �� -Z �� Z' ���Ӗ�����j
lmiterm([-lmi1 1 2  X     ], -c,   1);          % (1,2) �v�f�� -c*X *1 ��ǉ�
lmiterm([-lmi1 2 2  X     ],  r,   1);          % (2,2) �v�f��  r*X *1 ��ǉ�
% -------------------------------------------
lmi2 = newlmi;                                  % ���݋L�q����Ă��� LMI �Ɏ��ʎq lmi2 �� LMI ��ǉ��ilmi2�F����j
lmiterm([ lmi2 1 1  X     ],  A,   1, 's');     % (1,1) �v�f��   A*X *1 + ( A*X *1)' ��ǉ�
lmiterm([ lmi2 1 1 -Z     ],  B2,  1, 's');     % (1,1) �v�f��  B2*Z'*1 + (B2*Z'*1)' ��ǉ�
lmiterm([ lmi2 1 2  0     ],  B1);              % (1,2) �v�f��  B1 ��ǉ��E�E�E�萔�s��̂�
lmiterm([ lmi2 3 1  X     ],  C1,  1);          % (3,1) �v�f��  C1*X *1 ��ǉ�
lmiterm([ lmi2 3 1 -Z     ],  D12, 1);          % (3,1) �v�f�� D12*Z'*1 ��ǉ�
lmiterm([ lmi2 3 2  0     ],  D11);             % (3,2) �v�f�� D11 ��ǉ��E�E�E�萔�s��̂�
lmiterm([ lmi2 2 2  gamma ], -1,   1);          % (2,2) �v�f�� -1*gamma*1 ��ǉ�
lmiterm([ lmi2 3 3  gamma ], -1,   1);          % (3,3) �v�f��  1*gamma*1 ��ǉ�
% -------------------------------------------
lmisys = getlmis;                               % ��`���ꂽ LMI ���擾
% -------------------------------------------
cobj = zeros(1,decnbr(lmisys));                 % ����ϐ��̗v�f�����擾���C���̒����̗�x�N�g�� cobj �𐶐�
cobj(1) = 1;                                    % cobj �� 1 �Ԗڂ̗v�f�� 1 �ɕύX�i����ϐ��� 1 �Ԗڂ̗v�f�� gamma�j
% -------------------------------------------
[cost,xopt] = mincx(lmisys,cobj);               % �ړI�֐��� E = gamma �Ƃ����ʍœK����������
% -------------------------------------------
gamma_opt  = dec2mat(lmisys,xopt,gamma)         % ����ꂽ gamma �̍œK�� gamma_opt
X_opt = dec2mat(lmisys,xopt,X)                  % ����ꂽ X �̍œK�� X_opt
Z_opt = dec2mat(lmisys,xopt,Z)                  % ����ꂽ Z �̍œK�� Z_opt
% -------------------------------------------
K_opt = Z_opt'*inv(X_opt)                       % �R���g���[���Q�C�� K_opt