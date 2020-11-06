function [allocation_RB]=downlink_PF(num_rb,num_ue,CQI,windows_throuput_ue,num_symbol,num_subcarrier)       %������·������ƽ�㷨����OFDMA��Դ����䷽����
%%num_rb��ϵͳ��Դ�飬num_ue��ϵͳ���û���Ŀ��num_symbol��һ����֡��OFDM��������num_subcarriler��ÿ����֡Я�������ز���
%%CQI���û���ӦRB��CQI���� CQI��6*20��past_throuput_ue��ÿ���û�ʱ�䴰�ڵ�ƽ��������
%%allocation_RB�Ƿ������ ��һ��ֵ1�����RB�ѷ��䣬0δ���䣬�ڶ��д����RBʹ�õ�CQIֵ�������д����RB��������ĸ��û�
%% ������Ҫ�����ÿ���û���ÿ����Դ���������������ʣ�Ȼ������㷨���ȼ���ʽ�����ÿ���û���ÿ����Դ���ϵ����ȼ���Ȼ�����RB��ʼ����ֱ�����е�RBȫ���������
allocation_RB=zeros(num_rb,4);        %������󣬵�һ��ֵ1�����RB�ѷ��䣬0δ���䣬�ڶ��д����RBʹ�õ�CQIֵ�������д����RB��������ĸ��û�,�����д���ֵ����RB��Ĵ�������      
for rb_count = 1:1:num_rb
    K_argmax = zeros(1,num_ue);       % �洢һ��RB��ÿ���û��ϵ�Kֵ
    bit_rate_array = zeros(1,num_ue); % �洢��RB��ÿ���û��ı�������
    CQI_ue = zeros(1,num_ue);         % �洢��RB��ÿ���û���CQIֵ
    for ue_n= 1:1:num_ue              % ����ÿһ���û�
        CQI_ = CQI(rb_count, ue_n);   % �û���CQIֵ������ÿһ����Դ����˵��K��õ�
        CQI_ue(ue_n) = CQI_;          % �����û���CQIֵ
        [MCS_para,Qm,~]=modulation_CQI_mapping(CQI_);% �����û���CQI������������ʵĹ���
        bit_rate = (num_symbol-3)*num_subcarrier*Qm*MCS_para/1024; % PPT���û���������Դ���Ͽ��Ի�õ�������������
        bit_rate_array(ue_n)=bit_rate; % �洢��������
        if(windows_throuput_ue(ue_n)==0)
            windows_throuput_ue(ue_n)=0.00001;
        end       
        K_argmax(ue_n) = bit_rate/windows_throuput_ue(ue_n); % ���ݱ�����ƽ����������/ʱ�䴰��ƽ��������������ue�����ȼ���ֵԽ�����ȼ�Խ�� 
    %ѭ����ֱ�����������һ����Դ���ϵ����е��û���Kֵ��Ȼ���ٱ�����Դ�� 
    end
    ue_all=find(K_argmax==max( K_argmax));
    Sizeall=size(ue_all,2);
    one_ue=round(rand*(Sizeall-1)+1);
    ue_choose=ue_all(1,one_ue);
    %[~, ue_choose] = max(K_argmax);% ���������ƽ���ȼ���ߵ��û�Kֵ
    allocation_RB(rb_count, 1) = 1;       % ��������һ��Ϊ1��ʾ�ѷ���
    allocation_RB(rb_count, 2) = CQI_ue(ue_choose); % �������ڶ��б�ʾ��RB��������ȼ���ߵ��û�ʹ�õ�CQIֵ
    allocation_RB(rb_count, 3) = ue_choose; % ��¼��RB����������ȼ���ߵ��û�
    allocation_RB(rb_count, 4) = bit_rate_array(ue_choose); % �ֵ���Դ��������
end

