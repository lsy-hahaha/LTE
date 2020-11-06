function LTEconfig_add=lteparset2
LTEconfig_add.SubFrame.Number=50;%%ϵͳ������֡��
LTEconfig_add.System.Bandwidth=5e6;%%ϵͳ�ܴ���(Hz)
LTEconfig_add.bandwidth=180e3;%%ϵͳ������Դ��Ĵ���(Hz)
LTEconfig_add.PRB.deta_bandwidth=15e3;%%ϵͳ������Դ��ļ������(Hz)
LTEconfig_add.PRB.Number=25;%%%ϵͳ������Դ����
LTEconfig_add.Subcarrier.bandwidth=15e3;%ϵͳ���ز�����(Hz)
LTEconfig_add.Num_Subcarrier_eachPRB=12;
LTEconfig_add.Subcarrier.Number=12*LTEconfig_add.PRB.Number;%ϵͳ���ز���
LTEconfig_add.SpeedOfLight=3e8;%%����(m/s)
LTEconfig_add.inter_bts_distance=500;%%С���߳�
LTEconfig_add.number_of_sectors=3;%������С��������Ϊ3������
LTEconfig_add.eNB_max_antenna_gain=14;%��վ�������������
LTEconfig_add.data_res=20;%����߳�
LTEconfig_add.antenna_azimuth_offset= 60;%���߷���ĳ�ʼƫ�ýǶ� 
LTEconfig_add.SubFrame.Num_Symbol=14; %%%һ����֡��OFDMA������
LTEconfig_add.SubFrame.Duration=1e-3;%%%һ����֡������(s)

%% ϵͳ��UE�Ĳ������ %%
LTEconfig_add.UE.Num_TX=1;%%%UE�ķ���������Ŀ
LTEconfig_add.UE.Speed=2/3.6;%%%UE���ƶ��ٶ�
LTEconfig_add.UE.Number=10;%%ϵͳ��UE����Ŀ
LTEconfig_add.UE.thermal_noise_density = -174;%�¶�Ϊ290Kʱ���������������ܶȣ�dBm/Hz��
LTEconfig_add.UE.receiver_noise_figure = 9;    % Receiver noise figure in dB
LTEconfig_add.BS.Transmit_Power=40;%%%��վ���书��(W)
LTEconfig_add.BS.Transmit_Power_dB=16;%��λ16dB
%data_res=LTEconfig_add.data_res;
%% ��������С�� ������ %%
cell_side_length=LTEconfig_add.inter_bts_distance;
eNodeB = LTE_init_create_eNodeB(cell_side_length);
Neighborhood = zeros(6,2);%%���������λ������
Neighborhood(1,:) = eNodeB(1).pos;
Neighborhood(2,:) = eNodeB(2).pos;
Neighborhood(3,:) = eNodeB(3).pos;
Neighborhood(4,:) = eNodeB(6).pos;
Neighborhood(5,:) = eNodeB(7).pos;
Neighborhood(6,:) = eNodeB(4).pos;
LTEconfig_add.BS.Position_X=eNodeB(5).pos(1);
LTEconfig_add.BS.Position_Y=eNodeB(5).pos(2);
%% ��������С���û� %%
LTEconfig_add.roi_x = [min(Neighborhood(:,1)),max(Neighborhood(:,1))];%%������������Сֵ�����ֵ roi_x��1������Сֵ��roi_x(2)�����ֵ
LTEconfig_add.roi_y = [min(Neighborhood(:,2)),max(Neighborhood(:,2))];%%��������������Сֵ
[roi_max_pixel roi_pixel_exact] = LTE_common_pos_to_pixel([LTEconfig_add.roi_x(2) LTEconfig_add.roi_y(2)],[LTEconfig_add.roi_x(1) LTEconfig_add.roi_y(1)],LTEconfig_add.data_res);%%���������Ӧ��������
roi_height_pixel = roi_max_pixel(2);%%����
roi_width_pixel  = roi_max_pixel(1);%%����
 pos_grid_pixel = zeros(roi_height_pixel*roi_width_pixel,2);%%��������
 pos_grid_pixel(:,1) = reshape(repmat(1:roi_width_pixel,roi_height_pixel,1),1,roi_width_pixel*roi_height_pixel);
 pos_grid_pixel(:,2) = repmat(1:roi_height_pixel,1,roi_width_pixel);
 pos_grid_meter(:,1) = pos_grid_pixel(:,1)*LTEconfig_add.data_res-LTEconfig_add.roi_x(2);
 pos_grid_meter(:,2) = pos_grid_pixel(:,2)*LTEconfig_add.data_res-LTEconfig_add.roi_y(2);
Neighborhood = [Neighborhood;Neighborhood(1,:)];
 in = inpolygon(pos_grid_meter(:,1),pos_grid_meter(:,2),Neighborhood(:,1),Neighborhood(:,2));%%�淶����С����
 us_pos = zeros(sum(in),2);%%С�����û���
 b = 1;
 for i = 1 : length(in)
     if in(i) > 0
         us_pos(b,:) = pos_grid_meter(i,:);%%С�����û�����
         b = b + 1;
     end     
 end

a=randperm(b-2);%���˳��
ue_rand_pos=us_pos(a(1:LTEconfig_add.UE.Number),:);
distance=sqrt((ue_rand_pos(:,1)-LTEconfig_add.BS.Position_X).^2+(ue_rand_pos(:,2)-LTEconfig_add.BS.Position_Y).^2);
LTEconfig_add.UE.Number_Final=length(ue_rand_pos(:,1));
LTEconfig_add.UE.Position_X=ue_rand_pos(:,1);
LTEconfig_add.UE.Position_Y=ue_rand_pos(:,2); 
figure(3);
plot(Neighborhood(:,1),Neighborhood(:,2),LTEconfig_add.UE.Position_X,LTEconfig_add.UE.Position_Y,'.r');
hold on;
LTEconfig_add.UE.pixel(:,1) = round((ue_rand_pos(:,1) + LTEconfig_add.roi_x(2))./LTEconfig_add.data_res);%%�޵�λ��
LTEconfig_add.UE.pixel(:,2) = round((ue_rand_pos(:,2) + LTEconfig_add.roi_y(2))./LTEconfig_add.data_res);
%% ������С��
eNodeBs_7cell=LTE_init_create_eNodeB_7cell;
tx_pos = zeros(6,2);
tx_pos(1,:)= eNodeBs_7cell(1).pos;
tx_pos(2,:)= eNodeBs_7cell(2).pos;
tx_pos(3,:)= eNodeBs_7cell(3).pos;
tx_pos(4,:)= eNodeBs_7cell(6).pos;
tx_pos(5,:)= eNodeBs_7cell(7).pos;
tx_pos(6,:)= eNodeBs_7cell(4).pos;
%% ������С����վ������С���û����������� %%
angle=zeros(LTEconfig_add.number_of_sectors,LTEconfig_add.UE.Number_Final,6);
antenna = antennas.TS36942Antenna(LTEconfig_add.eNB_max_antenna_gain);
axis=zeros(3,2,6);%����վ���߲ο�����
antenna_gain=zeros(LTEconfig_add.UE.Number_Final,6);%�û�������վ����������
% for s_ = 1:LTEconfig_add.number_of_sectors
%     azimuth = wrapTo360(LTEconfig_add.antenna_azimuth_offset+ 120*(s_-1));%С�����߽�
% end
xy=[1,sqrt(3);-1,0;1,-sqrt(3)];
for i=1:6
    axis(:,1,i)=xy(:,1)+tx_pos(i,1);
    axis(:,2,i)=xy(:,2)+tx_pos(i,2);
end
for k=1:6
    for i=1:LTEconfig_add.number_of_sectors
        for j=1:LTEconfig_add.UE.Number_Final
    %����acosd(dot([x1-x2,y1-y2],[x3-x2,y3-y2])/(norm([x1-x2,y1-y2])*norm([x3-x2,y3-y2])))
        angle(i,j,k)=acosd(dot([axis(i,1,k)-tx_pos(k,1),axis(i,2,k)-tx_pos(k,2)],[LTEconfig_add.UE.Position_X(j)-tx_pos(k,1),LTEconfig_add.UE.Position_Y(j)-tx_pos(k,2)])/(norm([axis(i,1,k)-tx_pos(k,1),axis(i,2,k)-tx_pos(k,2)])*norm([LTEconfig_add.UE.Position_X(j)-tx_pos(k,1),LTEconfig_add.UE.Position_Y(j)-tx_pos(k,2)])));
        end
    end
end
for k=1:6
    angle_min=zeros(1,LTEconfig_add.UE.Number_Final);
    for i=1:LTEconfig_add.UE.Number_Final
        angle_min(i)=min(angle(:,i,k));
    end
    antenna_gain(:,k)=antenna.gain(angle_min);
end
%% �����û��յ���վ�ĸ��ţ�����С���м��ܵ���С����վ�ĸ���
distance_to_enbs=zeros(LTEconfig_add.UE.Number_Final,6);%�û���������С����վ�ľ���
passloss_to_enbs=zeros(LTEconfig_add.UE.Number_Final,6);
interference_to_enbs=zeros(LTEconfig_add.UE.Number_Final,6);
for i=1:6
    distance_to_enbs(:,i)=sqrt((LTEconfig_add.UE.Position_X-tx_pos(i,1)).^2+(LTEconfig_add.UE.Position_Y-tx_pos(i,2)).^2);
    passloss_to_enbs(:,i)=LTE_pathloss_enb_to_ue(distance_to_enbs(:,i));%�û�����С����վ��·�����
    interference_to_enbs(:,i)=10.^((LTEconfig_add.BS.Transmit_Power_dB+antenna_gain(:,i)-passloss_to_enbs(:,i))/10);%�û��ܵ���վ�ĸ���   
end
LTEconfig_add.ue_interference_enbs=zeros(LTEconfig_add.UE.Number_Final,1);
for ii=1:LTEconfig_add.UE.Number_Final
    LTEconfig_add.ue_interference_enbs(ii)=sum(interference_to_enbs(ii,:));
end
end
