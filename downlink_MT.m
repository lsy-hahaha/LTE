function [allocation_RB]=downlink_MT(num_rb,num_ue,CQI,num_symbol,num_subcarrier)       %������·����ظɱ��㷨����OFDMA��Դ����䷽����
%%num_rb��ϵͳ��Դ�飬num_ue��ϵͳ���û���Ŀ��num_symbol��һ����֡��OFDM��������num_subcarriler��ÿ����֡Я�������ز���
%%CQI���û���ӦRB��CQI���� CQI��6*20��
%%allocation_RB�Ƿ������ ��һ��ֵ1�����RB�ѷ��䣬0δ���䣬�ڶ��д����RBʹ�õ�CQIֵ�������д����RB��������ĸ��û�
%%������Ҫ�����ÿ���û���ÿ����Դ���������������ʣ�Ȼ������㷨���ȼ���ʽ�����ÿ���û���ÿ����Դ���ϵ����ȼ���Ȼ�����RB��ʼ����ֱ�����е�RBȫ���������
allocation_RB=zeros(num_rb,4);      % ��¼RB����λ�ü�����
for rb_count = 1:1:num_rb           % ����ÿһ��RB
    K_argmax = zeros(1,num_ue);     % ���ȼ�����
    CQI_ue = zeros(1,num_ue);       % ���ڴ洢�û���CQI
    %%��ʼ����ÿһ���û���CQI����K��ѡ��K���
    for ue_n= 1:1:num_ue            % ����ÿһ���û� 
        CQI_ = CQI(rb_count, ue_n); % �û���CQIֵ������ÿһ����Դ����˵���������û�����õ�                              
        CQI_ue(ue_n) = CQI_;        % �����û���CQIֵ
        [MCS_para,Qm,~]=modulation_CQI_mapping(CQI_);% �����û���CQI������������ʵĹ���
        bit_rate = (num_symbol-3)*num_subcarrier*Qm*MCS_para/1024;% PPT���û���������Դ���Ͽ��Ի�õ�������������
        K_argmax(ue_n) = bit_rate; % ����ue�����ȼ���ֵԽ�����ȼ�Խ��
    %ѭ����ֱ�����������һ����Դ���ϵ����е��û���Kֵ��Ȼ���ٱ�����Դ��  
    end
    [rate, ue_choose] = max(K_argmax, [], 2); % ��¼�ڸ���Դ�������ȼ���ߵ��û����ʺͱ��
                                    % ����RB
    allocation_RB(rb_count, 1) = 1; % ��������һ��Ϊ1��ʾ�ѷ���
	allocation_RB(rb_count, 2) = CQI_ue(ue_choose); % �������ڶ��б�ʾ��RB��������ȼ���ߵ��û�ʹ�õ�CQIֵ
    allocation_RB(rb_count, 3) = ue_choose; % ��¼��RB����������ȼ���ߵ��û�
    allocation_RB(rb_count, 4) = rate; % �ֵ���Դ��������
end
