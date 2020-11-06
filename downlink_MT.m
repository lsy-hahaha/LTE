function [allocation_RB]=downlink_MT(num_rb,num_ue,CQI,num_symbol,num_subcarrier)       %下行链路最大载干比算法（即OFDMA资源块分配方案）
%%num_rb是系统资源块，num_ue是系统中用户数目，num_symbol是一个子帧的OFDM符号数，num_subcarriler是每个子帧携带的子载波数
%%CQI是用户对应RB的CQI矩阵 CQI（6*20）
%%allocation_RB是分配矩阵 第一列值1代表该RB已分配，0未分配，第二列代表该RB使用的CQI值，第三列代表该RB分配给了哪个用户
%%首先需要计算出每个用户在每个资源块上理论数据速率，然后根据算法优先级公式计算出每个用户在每个资源块上的优先级，然后基于RB开始分配直至所有的RB全部分配完成
allocation_RB=zeros(num_rb,4);      % 记录RB分配位置及速率
for rb_count = 1:1:num_rb           % 遍历每一个RB
    K_argmax = zeros(1,num_ue);     % 优先级矩阵
    CQI_ue = zeros(1,num_ue);       % 用于存储用户的CQI
    %%开始根据每一个用户的CQI计算K，选择K大的
    for ue_n= 1:1:num_ue            % 遍历每一个用户 
        CQI_ = CQI(rb_count, ue_n); % 用户的CQI值，对于每一个资源块来说遍历所有用户找最好的                              
        CQI_ue(ue_n) = CQI_;        % 保存用户的CQI值
        [MCS_para,Qm,~]=modulation_CQI_mapping(CQI_);% 根据用户的CQI完成求数据速率的过程
        bit_rate = (num_symbol-3)*num_subcarrier*Qm*MCS_para/1024;% PPT：用户在物理资源块上可以获得的理论数据速率
        K_argmax(ue_n) = bit_rate; % 保存ue的优先级，值越大优先级越高
    %循环，直到计算出在这一个资源块上的所有的用户的K值，然后再遍历资源块  
    end
    [rate, ue_choose] = max(K_argmax, [], 2); % 记录在该资源块上优先级最高的用户速率和编号
                                    % 分配RB
    allocation_RB(rb_count, 1) = 1; % 分配矩阵第一列为1表示已分配
	allocation_RB(rb_count, 2) = CQI_ue(ue_choose); % 分配矩阵第二列表示该RB分配给优先级最高的用户使用的CQI值
    allocation_RB(rb_count, 3) = ue_choose; % 记录该RB分配给了优先级最高的用户
    allocation_RB(rb_count, 4) = rate; % 分到资源块后的速率
end
