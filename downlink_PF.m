function [allocation_RB]=downlink_PF(num_rb,num_ue,CQI,windows_throuput_ue,num_symbol,num_subcarrier)       %下行链路比例公平算法（即OFDMA资源块分配方案）
%%num_rb是系统资源块，num_ue是系统中用户数目，num_symbol是一个子帧的OFDM符号数，num_subcarriler是每个子帧携带的子载波数
%%CQI是用户对应RB的CQI矩阵 CQI（6*20）past_throuput_ue是每个用户时间窗内的平均吞吐量
%%allocation_RB是分配矩阵 第一列值1代表该RB已分配，0未分配，第二列代表该RB使用的CQI值，第三列代表该RB分配给了哪个用户
%% 首先需要计算出每个用户在每个资源块上理论数据速率，然后根据算法优先级公式计算出每个用户在每个资源块上的优先级，然后基于RB开始分配直至所有的RB全部分配完成
allocation_RB=zeros(num_rb,4);        %分配矩阵，第一列值1代表该RB已分配，0未分配，第二列代表该RB使用的CQI值，第三列代表该RB分配给了哪个用户,第四列代表分到这个RB后的传送速率      
for rb_count = 1:1:num_rb
    K_argmax = zeros(1,num_ue);       % 存储一个RB在每个用户上的K值
    bit_rate_array = zeros(1,num_ue); % 存储该RB下每个用户的比特速率
    CQI_ue = zeros(1,num_ue);         % 存储该RB下每个用户的CQI值
    for ue_n= 1:1:num_ue              % 遍历每一个用户
        CQI_ = CQI(rb_count, ue_n);   % 用户的CQI值，对于每一个资源块来说找K最好的
        CQI_ue(ue_n) = CQI_;          % 保存用户的CQI值
        [MCS_para,Qm,~]=modulation_CQI_mapping(CQI_);% 根据用户的CQI完成求数据速率的过程
        bit_rate = (num_symbol-3)*num_subcarrier*Qm*MCS_para/1024; % PPT：用户在物理资源块上可以获得的理论数据速率
        bit_rate_array(ue_n)=bit_rate; % 存储比特速率
        if(windows_throuput_ue(ue_n)==0)
            windows_throuput_ue(ue_n)=0.00001;
        end       
        K_argmax(ue_n) = bit_rate/windows_throuput_ue(ue_n); % 根据比例公平，比特速率/时间窗内平均吞吐量，保存ue的优先级，值越大优先级越高 
    %循环，直到计算出在这一个资源块上的所有的用户的K值，然后再遍历资源块 
    end
    ue_all=find(K_argmax==max( K_argmax));
    Sizeall=size(ue_all,2);
    one_ue=round(rand*(Sizeall-1)+1);
    ue_choose=ue_all(1,one_ue);
    %[~, ue_choose] = max(K_argmax);% 输出比例公平优先级最高的用户K值
    allocation_RB(rb_count, 1) = 1;       % 分配矩阵第一列为1表示已分配
    allocation_RB(rb_count, 2) = CQI_ue(ue_choose); % 分配矩阵第二列表示该RB分配给优先级最高的用户使用的CQI值
    allocation_RB(rb_count, 3) = ue_choose; % 记录该RB分配给了优先级最高的用户
    allocation_RB(rb_count, 4) = bit_rate_array(ue_choose); % 分到资源块后的速率
end

