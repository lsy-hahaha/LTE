function [allocation_RB,allocation_begin_ue_next]=downlink_RR(num_rb,num_ue,CQI,allocation_begin_ue,num_symbol,num_subcarrier)       %LTE������·��ѯ�㷨����OFDMA��Դ����䷽����
%%num_rb��ϵͳ��Դ�飬num_ue��ϵͳ���û���Ŀ��num_symbol��һ����֡��OFDM��������num_subcarriler��ÿ����֡Я�������ز���
%%CQI���û���ӦRB��CQI����  CQI=zeros(num_prb,num_ue),allocation_begin_ue���ڵ�ǰTTI��Ҫ��һ�����ܵ��ȵ��û����
%%allocation_RB�Ƿ������ ��һ��ֵ1�����RB�ѷ��䣬0δ���䣬�ڶ��д����RBʹ�õ�CQIֵ�������д����RB��������ĸ��û�
%%allocation_begin_ue_next����һ��TTI��Ҫ��һ�����ܵ��ȵ��û����
%% �ӵ�һ��RB��ʼ���䣬��һ��RB������ڵ�ǰTTI��Ҫ��һ�����ܵ��ȵ��û���Ȼ��ڶ���RB�������һ����Ҫ���ܵ��ȵ��û���ֱ������RB������ɣ���ʱ�ٷ�����һ��TTI��Ҫ���ܵ��ȵ��û����
    %��Դ��ķ���������һ��1����ֳ���0����δ�ֳ����ڶ��д���CQI��ֵ�������зָ��ĸ��û�,�����д���ֵ����RB��Ĵ�������
    allocation_RB=zeros(num_rb,4);    % ��¼RB����λ�ü�����
    %% ȷ����ʼλ��
    begin_ue = mod(allocation_begin_ue*num_rb, num_ue);   % ��ʼ����RB���û����
    %% ����
    tmp = 0;     % ��¼����ʣ��RB�ı��
	for rb_count = 1:1:num_rb         % ����ÿһ��RB
        ue_n = rb_count + begin_ue; % ʵ�ʷ�����û����
		if ue_n > num_ue            % �����õ��û������䵽��RB��ֹͣ����
			tmp = rb_count;           % ��¼����ʣ�µ�RB������ʣ��RB���н�һ������
        break                           
        end                           
        allocation_RB(rb_count, 1) = 1;              % ��������һ��ֵΪ1��ʾ�ѷ���
	    CQI_ = CQI(rb_count, ue_n);                % �û���CQIֵ
		allocation_RB(rb_count, 2) = CQI_;           % �������ڶ��б�ʾ��RBʹ�õ�CQIֵ
		allocation_RB(rb_count, 3) = ue_n;         % ��¼��RB��������ĸ��û�
		[MCS_para,Qm,~]=modulation_CQI_mapping(CQI_);% ���ú����������û���CQI������������ʵĹ���
		rate = (num_symbol-3)*num_subcarrier*Qm*MCS_para/1024; % PPT���û���������Դ���Ͽ��Ի�õ�������������
		allocation_RB(rb_count, 4) = rate;           % �������ĵ����б�ʾ������Դ��������
    end
    %%��ʼ����ʣ���RB
    if tmp
        for rb_count = tmp:1:num_rb       % ����ʣ�����Դ��
            ue_n = rb_count+1-tmp;      % +1��ʾ�ӵ�һ���û���ʼ���䣬��¼ʵʱ�����û����
            allocation_RB(rb_count, 1) = 1;% ��������һ��Ϊ1��ʾ�ѷ��䣻�ظ�֮ǰ�ķ�������¼����
            CQI_ = CQI(rb_count, ue_n);  % �û���CQIֵ
            allocation_RB(rb_count, 2) = CQI_;   % �������ڶ��б�ʾ��RBʹ�õ�CQIֵ
            allocation_RB(rb_count, 3) = ue_n; % ��¼��RB��������ĸ��û�
            [MCS_para,Qm,~]=modulation_CQI_mapping(CQI_); % �����û���CQI������������ʵĹ���
            rate = (num_symbol-3)*num_subcarrier*Qm*MCS_para/1024; % PPT���û���������Դ���Ͽ��Ի�õ�������������
            allocation_RB(rb_count, 4) = rate;            % ����
        end       
    end
    allocation_begin_ue_next = allocation_begin_ue + 1;   % ��¼��һ��TTI��ʼ������û�λ��
end