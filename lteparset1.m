function LTEconfig= lteparset1 %%系统仿真参数设置
LTEconfig.SubFrame.Number=50;%%系统仿真子帧数
LTEconfig.System.Bandwidth=5e6;%%系统总带宽(Hz)
LTEconfig.bandwidth=180e3;%%系统物理资源块的带宽(Hz)
LTEconfig.PRB.deta_bandwidth=15e3;%%系统物理资源块的间隔带宽(Hz)
LTEconfig.PRB.Number=25;%%%系统物理资源块数
LTEconfig.Subcarrier.bandwidth=15e3;%系统子载波带宽(Hz)
LTEconfig.Num_Subcarrier_eachPRB=12;
LTEconfig.Subcarrier.Number=12*LTEconfig.PRB.Number;%系统子载波数
LTEconfig.SpeedOfLight=3e8;%%光速(m/s)
LTEconfig.inter_bts_distance=500;%%小区边长
LTEconfig.number_of_sectors=3;%六边形小区被划分为3个扇区
LTEconfig.eNB_max_antenna_gain=14;%基站的最大天线增益
LTEconfig.data_res=20;%网格边长
LTEconfig.antenna_azimuth_offset= 60;%天线方向的初始偏置角度 
LTEconfig.SubFrame.Num_Symbol=14; %%%一个子帧的OFDMA符号数
LTEconfig.SubFrame.Duration=1e-3;%%%一个子帧的周期(s)
%% 系统中UE的参数设计 %%
LTEconfig.UE.Num_TX=1;%%%UE的发送天线数目
LTEconfig.UE.Speed=2/3.6;%%%UE的移动速度
LTEconfig.UE.Number=10;%%系统中UE的数目
LTEconfig.UE.thermal_noise_density = -174;%温度为290K时，热噪声功率谱密度（dBm/Hz）
LTEconfig.UE.receiver_noise_figure = 9;    % Receiver noise figure in dB
LTEconfig.BS.Transmit_Power=40;%%%基站发射功率(W)
LTEconfig.BS.Transmit_Power_dB=16;%单位16dB
%data_res=LTEconfig.data_res;
%% 布置中心小区 六边形 %%
cell_side_length=LTEconfig.inter_bts_distance;
eNodeB = LTE_init_create_eNodeB(cell_side_length);
Neighborhood = zeros(6,2);%%六个顶点的位置坐标
Neighborhood(1,:) = eNodeB(1).pos;
Neighborhood(2,:) = eNodeB(2).pos;
Neighborhood(3,:) = eNodeB(3).pos;
Neighborhood(4,:) = eNodeB(6).pos;
Neighborhood(5,:) = eNodeB(7).pos;
Neighborhood(6,:) = eNodeB(4).pos;
LTEconfig.BS.Position_X=eNodeB(5).pos(1);
LTEconfig.BS.Position_Y=eNodeB(5).pos(2);
%% 布置中心小区用户 %%
LTEconfig.roi_x = [min(Neighborhood(:,1)),max(Neighborhood(:,1))];%%顶点横坐标的最小值和最大值 roi_x（1）是最小值，roi_x(2)是最大值
LTEconfig.roi_y = [min(Neighborhood(:,2)),max(Neighborhood(:,2))];%%顶点纵坐标的最大、小值
[roi_max_pixel roi_pixel_exact] = LTE_common_pos_to_pixel([LTEconfig.roi_x(2) LTEconfig.roi_y(2)],[LTEconfig.roi_x(1) LTEconfig.roi_y(1)],LTEconfig.data_res);%%横向纵向对应的网格数
roi_height_pixel = roi_max_pixel(2);%%纵向
roi_width_pixel  = roi_max_pixel(1);%%横向
 pos_grid_pixel = zeros(roi_height_pixel*roi_width_pixel,2);%%网格坐标
 pos_grid_pixel(:,1) = reshape(repmat(1:roi_width_pixel,roi_height_pixel,1),1,roi_width_pixel*roi_height_pixel);
 pos_grid_pixel(:,2) = repmat(1:roi_height_pixel,1,roi_width_pixel);
 pos_grid_meter(:,1) = pos_grid_pixel(:,1)*LTEconfig.data_res-LTEconfig.roi_x(2);
 pos_grid_meter(:,2) = pos_grid_pixel(:,2)*LTEconfig.data_res-LTEconfig.roi_y(2);
Neighborhood = [Neighborhood;Neighborhood(1,:)];
 in = inpolygon(pos_grid_meter(:,1),pos_grid_meter(:,2),Neighborhood(:,1),Neighborhood(:,2));%%规范区域：小区内
 us_pos = zeros(sum(in),2);%%小区内用户数
 b = 1;
 for i = 1 : length(in)
     if in(i) > 0
         us_pos(b,:) = pos_grid_meter(i,:);%%小区内用户坐标
         b = b + 1;
     end     
 end
figure(1);
plot(Neighborhood(:,1),Neighborhood(:,2),us_pos(:,1),us_pos(:,2),'.r');
a=randperm(b-2);%随机顺序
ue_rand_pos=us_pos(a(1:LTEconfig.UE.Number),:);
distance=sqrt((ue_rand_pos(:,1)-LTEconfig.BS.Position_X).^2+(ue_rand_pos(:,2)-LTEconfig.BS.Position_Y).^2);
LTEconfig.UE.Number_Final=length(ue_rand_pos(:,1));
LTEconfig.UE.Position_X=ue_rand_pos(:,1);
LTEconfig.UE.Position_Y=ue_rand_pos(:,2); 
figure(2);
plot(Neighborhood(:,1),Neighborhood(:,2),LTEconfig.UE.Position_X,LTEconfig.UE.Position_Y,'.r');
LTEconfig.UE.pixel(:,1) = round((ue_rand_pos(:,1) + LTEconfig.roi_x(2))./LTEconfig.data_res);%%无单位制
LTEconfig.UE.pixel(:,2) = round((ue_rand_pos(:,2) + LTEconfig.roi_y(2))./LTEconfig.data_res);
%% 布置七小区
eNodeBs_7cell=LTE_init_create_eNodeB_7cell;
tx_pos = zeros(6,2);
tx_pos(1,:)= eNodeBs_7cell(1).pos;
tx_pos(2,:)= eNodeBs_7cell(2).pos;
tx_pos(3,:)= eNodeBs_7cell(3).pos;
tx_pos(4,:)= eNodeBs_7cell(6).pos;
tx_pos(5,:)= eNodeBs_7cell(7).pos;
tx_pos(6,:)= eNodeBs_7cell(4).pos;
%% 计算邻小区基站对中心小区用户的天线增益 %%
angle=zeros(LTEconfig.number_of_sectors,LTEconfig.UE.Number_Final,6);
antenna = antennas.TS36942Antenna(LTEconfig.eNB_max_antenna_gain);
axis=zeros(3,2,6);%各基站天线参考坐标
antenna_gain=zeros(LTEconfig.UE.Number_Final,6);%用户到各基站的天线增益
% for s_ = 1:LTEconfig.number_of_sectors
%     azimuth = wrapTo360(LTEconfig.antenna_azimuth_offset+ 120*(s_-1));%小区天线角
% end
xy=[1,sqrt(3);-1,0;1,-sqrt(3)];
for i=1:6
    axis(:,1,i)=xy(:,1)+tx_pos(i,1);
    axis(:,2,i)=xy(:,2)+tx_pos(i,2);
end
for k=1:6
    for i=1:LTEconfig.number_of_sectors
        for j=1:LTEconfig.UE.Number_Final
    %根据acosd(dot([x1-x2,y1-y2],[x3-x2,y3-y2])/(norm([x1-x2,y1-y2])*norm([x3-x2,y3-y2])))
        angle(i,j,k)=acosd(dot([axis(i,1,k)-tx_pos(k,1),axis(i,2,k)-tx_pos(k,2)],[LTEconfig.UE.Position_X(j)-tx_pos(k,1),LTEconfig.UE.Position_Y(j)-tx_pos(k,2)])/(norm([axis(i,1,k)-tx_pos(k,1),axis(i,2,k)-tx_pos(k,2)])*norm([LTEconfig.UE.Position_X(j)-tx_pos(k,1),LTEconfig.UE.Position_Y(j)-tx_pos(k,2)])));
        end
    end
end
for k=1:6
    angle_min=zeros(1,LTEconfig.UE.Number_Final);
    for i=1:LTEconfig.UE.Number_Final
        angle_min(i)=min(angle(:,i,k));
    end
    antenna_gain(:,k)=antenna.gain(angle_min);
end
%% 计算用户收到基站的干扰，中心小区中继受到邻小区基站的干扰
distance_to_enbs=zeros(LTEconfig.UE.Number_Final,6);%用户到六个邻小区基站的距离
passloss_to_enbs=zeros(LTEconfig.UE.Number_Final,6);
interference_to_enbs=zeros(LTEconfig.UE.Number_Final,6);
for i=1:6
    distance_to_enbs(:,i)=sqrt((LTEconfig.UE.Position_X-tx_pos(i,1)).^2+(LTEconfig.UE.Position_Y-tx_pos(i,2)).^2);
    passloss_to_enbs(:,i)=LTE_pathloss_enb_to_ue(distance_to_enbs(:,i));%用户到邻小区基站的路径损耗
    interference_to_enbs(:,i)=10.^((LTEconfig.BS.Transmit_Power_dB+antenna_gain(:,i)-passloss_to_enbs(:,i))/10);%用户受到基站的干扰   
end
LTEconfig.ue_interference_enbs=zeros(LTEconfig.UE.Number_Final,1);
for ii=1:LTEconfig.UE.Number_Final
    LTEconfig.ue_interference_enbs(ii)=sum(interference_to_enbs(ii,:));
end
end

