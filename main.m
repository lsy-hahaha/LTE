clc;
clear all;
num_simulate=20;%%仿真的次数
num_ue_set=9;
num_algorithm=3;
sum_throuput_all=zeros(num_algorithm,num_ue_set,num_simulate);
Jain_all=zeros(num_algorithm,num_ue_set,num_simulate);
sum_throuput_average=zeros(num_algorithm,num_ue_set);
Jain_average=zeros(num_algorithm,num_ue_set);
for iiii=1:1:num_simulate
close all;
 clear UE.distance_to_enb;
 clear UE.pathloss_to_enb;
 clear sf;
 clear ue_sf;
 clear UE.snr_to_enb;

%%  产生固定用户的信道  %%%%
LTEconfig=lteparset1;
num_prb=LTEconfig.PRB.Number;
num_ue=LTEconfig.UE.Number_Final;%%%用户数目
num_subframe=LTEconfig.SubFrame.Number;
num_symbol=LTEconfig.SubFrame.Num_Symbol;
num_subcarrier=LTEconfig.Num_Subcarrier_eachPRB;
%% 计算发送给用户的天线增益
angle=zeros(LTEconfig.number_of_sectors,num_ue);
angle_min=zeros(1,num_ue);
antenna = antennas.TS36942Antenna(LTEconfig.eNB_max_antenna_gain);
for s_ = 1:LTEconfig.number_of_sectors
    azimuth = wrapTo360(LTEconfig.antenna_azimuth_offset+ 120*(s_-1));%小区天线角
end
%xy=[sqrt(3),1;-sqrt(3),1;0,-1];%对应30度偏角
xy=[1,sqrt(3);-1,0;1,-sqrt(3)];%对应60度
for i=1:LTEconfig.number_of_sectors
    for j=1:num_ue
%根据acosd(dot([x1-x2,y1-y2],[x3-x2,y3-y2])/(norm([x1-x2,y1-y2])*norm([x3-x2,y3-y2])))
    angle(i,j)=acosd(dot([xy(i,1)-LTEconfig.BS.Position_X,xy(i,2)-LTEconfig.BS.Position_Y],[LTEconfig.UE.Position_X(j)-LTEconfig.BS.Position_X,LTEconfig.UE.Position_Y(j)-LTEconfig.BS.Position_Y])/(norm([xy(i,1)-LTEconfig.BS.Position_X,xy(i,2)-LTEconfig.BS.Position_Y])*norm([LTEconfig.UE.Position_X(j)-LTEconfig.BS.Position_X,LTEconfig.UE.Position_Y(j)-LTEconfig.BS.Position_Y])));
    end
end
for i=1:num_ue
    [angle_min(i),sector]=min(angle(:,i));
end
distance.sector_antenna_gain=antenna.gain(angle_min);
%% 计算路径损耗
UE.distance_to_enb=sqrt((LTEconfig.UE.Position_X-LTEconfig.BS.Position_X).^2+(LTEconfig.UE.Position_Y-LTEconfig.BS.Position_Y).^2);
UE.pathloss_to_enb= LTE_pathloss_enb_to_ue(UE.distance_to_enb);%num_ue*1
%% 计算阴影衰落
sf=shadowfading(LTEconfig.roi_x,LTEconfig.roi_y,LTEconfig.data_res);%每个网格对应的阴影衰落 
for i=1:num_ue
    ue_sf(i,1)=sf(LTEconfig.UE.pixel(i,1),LTEconfig.UE.pixel(i,2));%num_ue个用户的阴影衰落 num_ue*1
end
LTEconfig.UE.receiver_noise_figure = 9;    % Receiver noise figure in dB
thermal_noise_W=10^(LTEconfig.UE.thermal_noise_density/10) / 1000 * LTEconfig.System.Bandwidth* 10^(LTEconfig.UE.receiver_noise_figure/10);%系统中的热噪声
UE.snr_to_enb=(10.^((LTEconfig.BS.Transmit_Power_dB+distance.sector_antenna_gain'-UE.pathloss_to_enb-ue_sf)/10))./(thermal_noise_W+ LTEconfig.ue_interference_enbs);%num_ue*1
UE.snr_to_enb_db=10*log10(UE.snr_to_enb);

%% 产生随时间变化的小尺度衰落信道
 SimType=8;   % 设置仿真类型
ltepar=lteparset(SimType,num_ue); % 配置仿真参数
ltepar.Nrb=num_prb; %系统PRB个数
ltepar.nUE=num_ue; %用户数目
 ltepar.Ntot=300;   %占有子载波数
ltepar.Bandwidth=5e6; %系统带宽
ltepar.Nfft=512; %FFT点数
ltepar.Fs=7.68e6; %采样频率
 h=cell(1,ltepar.nUE,num_subframe);
 sc_snr=cell(1,ltepar.nUE);%%元细胞数组
 prb_snr=zeros(num_prb,ltepar.nUE);
 prb_snr_TTI=zeros(num_prb,ltepar.nUE,num_subframe);
subframe_i=1;
while subframe_i<=num_subframe
     chan=cell(1,ltepar.nUE);
    for k=1:ltepar.nUE
      [chan{k}, h{1,k,subframe_i}]=channelmodel(ltepar.ChanModPar.PDP_dB,ltepar.Fs,ltepar.CarrierFrequency*1e9,...
      ltepar.UEpar.speed,ltepar.ChanModPar.FadeType,ltepar.ChanModPar.CorrType,ltepar.ChanModPar.M,ltepar.psi(:,:,:,:,k),ltepar.theta(:,:,:,:,k),...
      ltepar.phi(:,:,:,:,k),subframe_i*1e-3,ltepar.UEpar.nTX,ltepar.BSpar.nRX,ltepar.ChanModPar.TxCorrMatrix,ltepar.ChanModPar.RxCorrMatrix);
      sc_snr{1,k}=calculate_SC_SNR(h{1,k,subframe_i},UE.snr_to_enb_db(k)); 
    end
    prb_snr=calculate_PRB_SNR(sc_snr,num_prb);
    prb_snr_TTI(:,:,subframe_i)=prb_snr;
    subframe_i=subframe_i+1;
end

%% 产生增加用户的信道%%%%
for jjj=1:1:num_ue_set
    clear UE.distance_to_enb;
    clear UE.pathloss_to_enb;
    clear sf;
    clear ue_sf;
    clear UE.snr_to_enb;
LTEconfig_add=lteparset2;
num_ue=LTEconfig_add.UE.Number_Final;%%%用户数目
%% 计算发送给用户的天线增益
angle=zeros(LTEconfig_add.number_of_sectors,num_ue);
angle_min=zeros(1,num_ue);
antenna = antennas.TS36942Antenna(LTEconfig_add.eNB_max_antenna_gain);
for s_ = 1:LTEconfig_add.number_of_sectors
    azimuth = wrapTo360(LTEconfig_add.antenna_azimuth_offset+ 120*(s_-1));%小区天线角
end
%xy=[sqrt(3),1;-sqrt(3),1;0,-1];%对应30度偏角
xy=[1,sqrt(3);-1,0;1,-sqrt(3)];%对应60度
for i=1:LTEconfig_add.number_of_sectors
    for j=1:num_ue
%根据acosd(dot([x1-x2,y1-y2],[x3-x2,y3-y2])/(norm([x1-x2,y1-y2])*norm([x3-x2,y3-y2])))
    angle(i,j)=acosd(dot([xy(i,1)-LTEconfig_add.BS.Position_X,xy(i,2)-LTEconfig_add.BS.Position_Y],[LTEconfig_add.UE.Position_X(j)-LTEconfig_add.BS.Position_X,LTEconfig_add.UE.Position_Y(j)-LTEconfig_add.BS.Position_Y])/(norm([xy(i,1)-LTEconfig_add.BS.Position_X,xy(i,2)-LTEconfig_add.BS.Position_Y])*norm([LTEconfig_add.UE.Position_X(j)-LTEconfig_add.BS.Position_X,LTEconfig_add.UE.Position_Y(j)-LTEconfig_add.BS.Position_Y])));
    end
end
for i=1:num_ue
    [angle_min(i),sector]=min(angle(:,i));
end
distance.sector_antenna_gain=antenna.gain(angle_min);
%% 计算路径损耗
UE.distance_to_enb=sqrt((LTEconfig_add.UE.Position_X-LTEconfig_add.BS.Position_X).^2+(LTEconfig_add.UE.Position_Y-LTEconfig_add.BS.Position_Y).^2);
UE.pathloss_to_enb= LTE_pathloss_enb_to_ue(UE.distance_to_enb);%num_ue*1
%% 计算阴影衰落
sf=shadowfading(LTEconfig_add.roi_x,LTEconfig_add.roi_y,LTEconfig_add.data_res);%每个网格对应的阴影衰落 
for i=1:num_ue
    ue_sf(i,1)=sf(LTEconfig_add.UE.pixel(i,1),LTEconfig_add.UE.pixel(i,2));%num_ue个用户的阴影衰落 num_ue*1
end
LTEconfig_add.UE.receiver_noise_figure = 9;    % Receiver noise figure in dB
thermal_noise_W=10^(LTEconfig_add.UE.thermal_noise_density/10) / 1000 * LTEconfig_add.System.Bandwidth* 10^(LTEconfig_add.UE.receiver_noise_figure/10);%系统中的热噪声
UE.snr_to_enb=(10.^((LTEconfig_add.BS.Transmit_Power_dB+distance.sector_antenna_gain'-UE.pathloss_to_enb-ue_sf)/10))./(thermal_noise_W+ LTEconfig_add.ue_interference_enbs);%num_ue*1
UE.snr_to_enb_db=10*log10(UE.snr_to_enb);

%% 产生随时间变化的小尺度衰落信道
 SimType=8;   % 设置仿真类型
ltepar=lteparset(SimType,num_ue); % 配置仿真参数
ltepar.Nrb=num_prb; %系统PRB个数
ltepar.nUE=num_ue; %用户数目
 ltepar.Ntot=300;   %占有子载波数
ltepar.Bandwidth=5e6; %系统带宽
ltepar.Nfft=512; %FFT点数
ltepar.Fs=7.68e6; %采样频率
 h=cell(1,ltepar.nUE,num_subframe);
 sc_snr=cell(1,ltepar.nUE);%%元细胞数组
 prb_snr=zeros(num_prb,ltepar.nUE);
 prb_snr_TTI_add=zeros(num_prb,ltepar.nUE,num_subframe);
subframe_i=1;
while subframe_i<=num_subframe
     chan=cell(1,ltepar.nUE);
    for k=1:ltepar.nUE
      [chan{k}, h{1,k,subframe_i}]=channelmodel(ltepar.ChanModPar.PDP_dB,ltepar.Fs,ltepar.CarrierFrequency*1e9,...
      ltepar.UEpar.speed,ltepar.ChanModPar.FadeType,ltepar.ChanModPar.CorrType,ltepar.ChanModPar.M,ltepar.psi(:,:,:,:,k),ltepar.theta(:,:,:,:,k),...
      ltepar.phi(:,:,:,:,k),subframe_i*1e-3,ltepar.UEpar.nTX,ltepar.BSpar.nRX,ltepar.ChanModPar.TxCorrMatrix,ltepar.ChanModPar.RxCorrMatrix);
      sc_snr{1,k}=calculate_SC_SNR(h{1,k,subframe_i},UE.snr_to_enb_db(k)); 
    end
    prb_snr=calculate_PRB_SNR(sc_snr,num_prb);
    prb_snr_TTI_add(:,:,subframe_i)=prb_snr;
    subframe_i=subframe_i+1;
end

%% 原有的用户信道与新增加的用户信道合在一起
LTEconfig.UE.Number_Final= LTEconfig.UE.Number_Final+LTEconfig_add.UE.Number_Final;
LTEconfig.UE.Position_X=cat(1,LTEconfig.UE.Position_X,LTEconfig_add.UE.Position_X);
LTEconfig.UE.Position_Y=cat(1,LTEconfig.UE.Position_Y,LTEconfig_add.UE.Position_Y);
LTEconfig.UE.pixel=cat(1,LTEconfig.UE.pixel,LTEconfig_add.UE.pixel);
LTEconfig.ue_interference_enbs=cat(1,LTEconfig.ue_interference_enbs,LTEconfig_add.ue_interference_enbs);
prb_snr_TTI=cat(2,prb_snr_TTI,prb_snr_TTI_add);
num_ue=LTEconfig.UE.Number_Final;%%%用户数目

%% 变量定义初始化
sum_throuput_ue_perTTI=zeros(num_ue,num_subframe,num_algorithm);%每个TTI内每个用户的吞吐量
sum_throuput_ue_until_TTI=zeros(num_ue,num_algorithm);%所有用户到当前TTI的累积吞吐量
sum_throuput_total_all_TTI_RT=zeros(num_algorithm,1);%所有用户在所有TTI内的吞吐量
sum_throuput_system_kbps_RT=zeros(num_algorithm,1);%系统吞吐量
Jain=zeros(num_algorithm,1);
for iii=1:1:3
        if iii==1
            allocation_type='RR';
        elseif iii==2
            allocation_type='MT';
        elseif  iii==3
            allocation_type='PF';
        end
       allocation_begin_ue=zeros(1,num_subframe)+1;%分配用户初始化（仅用于RR算法）
      CQI=zeros(num_prb,num_ue); %CQI矩阵置0
      windows_band=5;%窗长
      windows_throuput_ue=zeros(1,num_ue); %每个用户的窗吞吐量
      past_throuput_ue=zeros(1,num_ue);
      for TTI_ID=1:1:num_subframe
             CQI=SINR_mapping_CQI(prb_snr_TTI(:,:,1:TTI_ID),TTI_ID,num_ue,num_prb);


             %%  调度算法
             if(strcmp(allocation_type,'RR'))%选择三种算法    strcmp是比较算法
                [allocation_RB,allocation_begin_ue(1,TTI_ID+1)]=downlink_RR(num_prb,num_ue,CQI(:,:),allocation_begin_ue(1,TTI_ID),num_symbol,num_subcarrier); 
            elseif (strcmp(allocation_type,'PF'))
                [allocation_RB]=downlink_PF(num_prb,num_ue,CQI(:,:),windows_throuput_ue,num_symbol,num_subcarrier);   
                
            elseif (strcmp(allocation_type,'MT'))
                [allocation_RB]=downlink_MT(num_prb,num_ue,CQI(:,:),num_symbol,num_subcarrier);
            else error('NOT support this schedule way!');
             end
           
            %% 计算每个TTI内每个用户的吞吐量
            for jj=1:1:num_ue
                 ue_prb=find(allocation_RB(:,3)==jj);
                 num_ue_prb=length(ue_prb);
                 for i=1:1:num_ue_prb
                   sum_throuput_ue_perTTI(jj,TTI_ID,iii)=sum_throuput_ue_perTTI(jj,TTI_ID,iii)+allocation_RB(ue_prb(i),4);
                 end
            end
           
            %% 计算每个用户到当前TTI的累积吞吐量
            for jj=1:1:num_ue
               sum_throuput_ue_until_TTI(jj,iii)=sum_throuput_ue_until_TTI(jj,iii)+sum_throuput_ue_perTTI(jj,TTI_ID,iii);
            end
        
       
            %% 计算用户的窗吞吐量
             for jj=1:1:num_ue    
                if (TTI_ID<windows_band+1)
                    past_throuput_ue(jj)=sum_throuput_ue_until_TTI(jj,iii); %当前TTI不足窗长，当前TTI之前的吞吐量作为该用户的窗吞吐量
                elseif (TTI_ID>windows_band)
                    past_throuput_ue(jj)=past_throuput_ue(jj)+sum_throuput_ue_perTTI(jj,TTI_ID,iii)-sum_throuput_ue_perTTI(jj,TTI_ID-windows_band,iii);%计算TTI大于窗长以后的每个用户的窗吞吐量
                end
             windows_throuput_ue(jj)=(1-1/windows_band)*past_throuput_ue(jj)+(1/windows_band)*sum_throuput_ue_perTTI(jj,TTI_ID,iii);
             end   
       
end
         %% 计算实时业务用户在当前分配方式下的系统吞吐量
         sum_throuput_system_kbps_RT(iii,1)=sum(sum_throuput_ue_until_TTI(:,iii))/num_subframe; 
         %%  计算所有实时业务用户Jain指数
        Jain(iii,1)=(sum(sum_throuput_ue_until_TTI(:,iii)))^2/num_ue/sum(sum_throuput_ue_until_TTI(:,iii).^2);
end

if num_ue==20
   save 20UE_sum_throuput_system_kbps_RT.txt sum_throuput_system_kbps_RT -ascii;
   save 20UE_Jain.txt Jain -ascii;
elseif num_ue==30
   save 30UE_sum_throuput_system_kbps_RT.txt sum_throuput_system_kbps_RT -ascii;
   save 30UE_Jain.txt Jain -ascii;   
elseif num_ue==40
   save 40UE_sum_throuput_system_kbps_RT.txt sum_throuput_system_kbps_RT -ascii;
   save 40UE_Jain.txt Jain -ascii;    
elseif num_ue==50
   save 50UE_sum_throuput_system_kbps_RT.txt sum_throuput_system_kbps_RT -ascii;
   save 50UE_Jain.txt Jain -ascii;   
elseif num_ue==60
   save 60UE_sum_throuput_system_kbps_RT.txt sum_throuput_system_kbps_RT -ascii;
   save 60UE_Jain.txt Jain -ascii;   
elseif num_ue==70
   save 70UE_sum_throuput_system_kbps_RT.txt sum_throuput_system_kbps_RT -ascii;
   save 70UE_Jain.txt Jain -ascii;   
elseif num_ue==80
   save 80UE_sum_throuput_system_kbps_RT.txt sum_throuput_system_kbps_RT -ascii;
   save 80UE_Jain.txt Jain -ascii; 
elseif num_ue==90
   save 90UE_sum_throuput_system_kbps_RT.txt sum_throuput_system_kbps_RT -ascii;
   save 90UE_Jain.txt Jain -ascii; 
elseif num_ue==100
   save 100UE_sum_throuput_system_kbps_RT.txt sum_throuput_system_kbps_RT -ascii;
   save 100UE_Jain.txt Jain -ascii; 
end
end
sum_throuput_all(:,1,iiii)=load('20UE_sum_throuput_system_kbps_RT.txt');
Jain_all(:,1,iiii)=load('20UE_Jain.txt');

sum_throuput_all(:,2,iiii)=load('30UE_sum_throuput_system_kbps_RT.txt');
Jain_all(:,2,iiii)=load('30UE_Jain.txt');

sum_throuput_all(:,3,iiii)=load('40UE_sum_throuput_system_kbps_RT.txt');
Jain_all(:,3,iiii)=load('40UE_Jain.txt');

sum_throuput_all(:,4,iiii)=load('50UE_sum_throuput_system_kbps_RT.txt');
Jain_all(:,4,iiii)=load('50UE_Jain.txt');

sum_throuput_all(:,5,iiii)=load('60UE_sum_throuput_system_kbps_RT.txt');
Jain_all(:,5,iiii)=load('60UE_Jain.txt');

sum_throuput_all(:,6,iiii)=load('70UE_sum_throuput_system_kbps_RT.txt');
Jain_all(:,6,iiii)=load('70UE_Jain.txt');

sum_throuput_all(:,7,iiii)=load('80UE_sum_throuput_system_kbps_RT.txt');
Jain_all(:,7,iiii)=load('80UE_Jain.txt');

sum_throuput_all(:,8,iiii)=load('90UE_sum_throuput_system_kbps_RT.txt');
Jain_all(:,8,iiii)=load('90UE_Jain.txt');

sum_throuput_all(:,9,iiii)=load('100UE_sum_throuput_system_kbps_RT.txt');
Jain_all(:,9,iiii)=load('100UE_Jain.txt');


end
for iii=1:1:num_algorithm
    for jjj=1:1:num_ue_set
        sum_throuput_average(iii,jjj)=mean(sum_throuput_all(iii,jjj,:))/1000;
        Jain_average(iii,jjj)=mean(Jain_all(iii,jjj,:));
    end
end

figure(4);
t=20:10:100;
plot(t,sum_throuput_average(1,:),'-s',t,sum_throuput_average(2,:),'-p',t,sum_throuput_average(3,:),'-o','linewidth',2);
legend('RR','MT','PF');
xlabel('系统接入非实时业务用户数（个）');
ylabel('系统吞吐量（Mbps）');
grid on;
hold on;



figure(5);
t=20:10:100;
plot(t,Jain_average(1,:),'-s',t,Jain_average(2,:),'-p',t,Jain_average(3,:),'-o','linewidth',2);
legend('RR','MT','PF');
xlabel('系统接入非实时业务用户数（个）');
ylabel('Jains指数');
grid on;
hold on;

