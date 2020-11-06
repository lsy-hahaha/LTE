function [allocation_RB,allocation_begin_ue_next]=downlink_RR(num_rb,num_ue,CQI,allocation_begin_ue,num_symbol,num_subcarrier)       %LTE下行链路轮询算法（即OFDMA资源块分配方案）
%%num_rb是系统资源块，num_ue是系统中用户数目，num_symbol是一个子帧的OFDM符号数，num_subcarriler是每个子帧携带的子载波数
%%CQI是用户对应RB的CQI矩阵  CQI=zeros(num_prb,num_ue),allocation_begin_ue是在当前TTI需要第一个接受调度的用户编号
%%allocation_RB是分配矩阵 第一列值1代表该RB已分配，0未分配，第二列代表该RB使用的CQI值，第三列代表该RB分配给了哪个用户
%%allocation_begin_ue_next是下一个TTI需要第一个接受调度的用户编号
%% 从第一个RB开始分配，第一个RB分配给在当前TTI需要第一个接受调度的用户，然后第二个RB分配给下一个需要接受调度的用户，直至所有RB分配完成，此时再返回下一个TTI需要接受调度的用户编号
    %资源块的分配结果，第一列1代表分出，0代表未分出，第二列代表CQI的值，第三列分给哪个用户,第四列代表分到这个RB后的传送速率
    allocation_RB=zeros(num_rb,4);    % 记录RB分配位置及速率
    %% 确定起始位置
    begin_ue = mod(allocation_begin_ue*num_rb, num_ue);   % 开始分配RB的用户编号
    %% 分配
    tmp = 0;     % 记录分配剩下RB的编号
	for rb_count = 1:1:num_rb         % 遍历每一个RB
        ue_n = rb_count + begin_ue; % 实际分配的用户编号
		if ue_n > num_ue            % 当所用的用户都分配到了RB，停止分配
			tmp = rb_count;           % 记录分配剩下的RB，并将剩下RB进行进一步分配
        break                           
        end                           
        allocation_RB(rb_count, 1) = 1;              % 分配矩阵第一列值为1表示已分配
	    CQI_ = CQI(rb_count, ue_n);                % 用户的CQI值
		allocation_RB(rb_count, 2) = CQI_;           % 分配矩阵第二列表示该RB使用的CQI值
		allocation_RB(rb_count, 3) = ue_n;         % 记录该RB分配给了哪个用户
		[MCS_para,Qm,~]=modulation_CQI_mapping(CQI_);% 调用函数，根据用户的CQI完成求数据速率的过程
		rate = (num_symbol-3)*num_subcarrier*Qm*MCS_para/1024; % PPT：用户在物理资源块上可以获得的理论数据速率
		allocation_RB(rb_count, 4) = rate;           % 分配矩阵的第四列表示分配资源块后的速率
    end
    %%开始分配剩余的RB
    if tmp
        for rb_count = tmp:1:num_rb       % 遍历剩余的资源块
            ue_n = rb_count+1-tmp;      % +1表示从第一个用户开始分配，记录实时分配用户编号
            allocation_RB(rb_count, 1) = 1;% 分配矩阵第一列为1表示已分配；重复之前的分配矩阵记录过程
            CQI_ = CQI(rb_count, ue_n);  % 用户的CQI值
            allocation_RB(rb_count, 2) = CQI_;   % 分配矩阵第二列表示该RB使用的CQI值
            allocation_RB(rb_count, 3) = ue_n; % 记录该RB分配给了哪个用户
            [MCS_para,Qm,~]=modulation_CQI_mapping(CQI_); % 根据用户的CQI完成求数据速率的过程
            rate = (num_symbol-3)*num_subcarrier*Qm*MCS_para/1024; % PPT：用户在物理资源块上可以获得的理论数据速率
            allocation_RB(rb_count, 4) = rate;            % 速率
        end       
    end
    allocation_begin_ue_next = allocation_begin_ue + 1;   % 记录下一个TTI开始分配的用户位置
end