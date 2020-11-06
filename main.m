clc;
clear all;
num_simulate=20;%%����Ĵ���
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

%%  �����̶��û����ŵ�  %%%%
LTEconfig=lteparset1;
num_prb=LTEconfig.PRB.Number;
num_ue=LTEconfig.UE.Number_Final;%%%�û���Ŀ
num_subframe=LTEconfig.SubFrame.Number;
num_symbol=LTEconfig.SubFrame.Num_Symbol;
num_subcarrier=LTEconfig.Num_Subcarrier_eachPRB;
%% ���㷢�͸��û�����������
angle=zeros(LTEconfig.number_of_sectors,num_ue);
angle_min=zeros(1,num_ue);
antenna = antennas.TS36942Antenna(LTEconfig.eNB_max_antenna_gain);
for s_ = 1:LTEconfig.number_of_sectors
    azimuth = wrapTo360(LTEconfig.antenna_azimuth_offset+ 120*(s_-1));%С�����߽�
end
%xy=[sqrt(3),1;-sqrt(3),1;0,-1];%��Ӧ30��ƫ��
xy=[1,sqrt(3);-1,0;1,-sqrt(3)];%��Ӧ60��
for i=1:LTEconfig.number_of_sectors
    for j=1:num_ue
%����acosd(dot([x1-x2,y1-y2],[x3-x2,y3-y2])/(norm([x1-x2,y1-y2])*norm([x3-x2,y3-y2])))
    angle(i,j)=acosd(dot([xy(i,1)-LTEconfig.BS.Position_X,xy(i,2)-LTEconfig.BS.Position_Y],[LTEconfig.UE.Position_X(j)-LTEconfig.BS.Position_X,LTEconfig.UE.Position_Y(j)-LTEconfig.BS.Position_Y])/(norm([xy(i,1)-LTEconfig.BS.Position_X,xy(i,2)-LTEconfig.BS.Position_Y])*norm([LTEconfig.UE.Position_X(j)-LTEconfig.BS.Position_X,LTEconfig.UE.Position_Y(j)-LTEconfig.BS.Position_Y])));
    end
end
for i=1:num_ue
    [angle_min(i),sector]=min(angle(:,i));
end
distance.sector_antenna_gain=antenna.gain(angle_min);
%% ����·�����
UE.distance_to_enb=sqrt((LTEconfig.UE.Position_X-LTEconfig.BS.Position_X).^2+(LTEconfig.UE.Position_Y-LTEconfig.BS.Position_Y).^2);
UE.pathloss_to_enb= LTE_pathloss_enb_to_ue(UE.distance_to_enb);%num_ue*1
%% ������Ӱ˥��
sf=shadowfading(LTEconfig.roi_x,LTEconfig.roi_y,LTEconfig.data_res);%ÿ�������Ӧ����Ӱ˥�� 
for i=1:num_ue
    ue_sf(i,1)=sf(LTEconfig.UE.pixel(i,1),LTEconfig.UE.pixel(i,2));%num_ue���û�����Ӱ˥�� num_ue*1
end
LTEconfig.UE.receiver_noise_figure = 9;    % Receiver noise figure in dB
thermal_noise_W=10^(LTEconfig.UE.thermal_noise_density/10) / 1000 * LTEconfig.System.Bandwidth* 10^(LTEconfig.UE.receiver_noise_figure/10);%ϵͳ�е�������
UE.snr_to_enb=(10.^((LTEconfig.BS.Transmit_Power_dB+distance.sector_antenna_gain'-UE.pathloss_to_enb-ue_sf)/10))./(thermal_noise_W+ LTEconfig.ue_interference_enbs);%num_ue*1
UE.snr_to_enb_db=10*log10(UE.snr_to_enb);

%% ������ʱ��仯��С�߶�˥���ŵ�
 SimType=8;   % ���÷�������
ltepar=lteparset(SimType,num_ue); % ���÷������
ltepar.Nrb=num_prb; %ϵͳPRB����
ltepar.nUE=num_ue; %�û���Ŀ
 ltepar.Ntot=300;   %ռ�����ز���
ltepar.Bandwidth=5e6; %ϵͳ����
ltepar.Nfft=512; %FFT����
ltepar.Fs=7.68e6; %����Ƶ��
 h=cell(1,ltepar.nUE,num_subframe);
 sc_snr=cell(1,ltepar.nUE);%%Ԫϸ������
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

%% ���������û����ŵ�%%%%
for jjj=1:1:num_ue_set
    clear UE.distance_to_enb;
    clear UE.pathloss_to_enb;
    clear sf;
    clear ue_sf;
    clear UE.snr_to_enb;
LTEconfig_add=lteparset2;
num_ue=LTEconfig_add.UE.Number_Final;%%%�û���Ŀ
%% ���㷢�͸��û�����������
angle=zeros(LTEconfig_add.number_of_sectors,num_ue);
angle_min=zeros(1,num_ue);
antenna = antennas.TS36942Antenna(LTEconfig_add.eNB_max_antenna_gain);
for s_ = 1:LTEconfig_add.number_of_sectors
    azimuth = wrapTo360(LTEconfig_add.antenna_azimuth_offset+ 120*(s_-1));%С�����߽�
end
%xy=[sqrt(3),1;-sqrt(3),1;0,-1];%��Ӧ30��ƫ��
xy=[1,sqrt(3);-1,0;1,-sqrt(3)];%��Ӧ60��
for i=1:LTEconfig_add.number_of_sectors
    for j=1:num_ue
%����acosd(dot([x1-x2,y1-y2],[x3-x2,y3-y2])/(norm([x1-x2,y1-y2])*norm([x3-x2,y3-y2])))
    angle(i,j)=acosd(dot([xy(i,1)-LTEconfig_add.BS.Position_X,xy(i,2)-LTEconfig_add.BS.Position_Y],[LTEconfig_add.UE.Position_X(j)-LTEconfig_add.BS.Position_X,LTEconfig_add.UE.Position_Y(j)-LTEconfig_add.BS.Position_Y])/(norm([xy(i,1)-LTEconfig_add.BS.Position_X,xy(i,2)-LTEconfig_add.BS.Position_Y])*norm([LTEconfig_add.UE.Position_X(j)-LTEconfig_add.BS.Position_X,LTEconfig_add.UE.Position_Y(j)-LTEconfig_add.BS.Position_Y])));
    end
end
for i=1:num_ue
    [angle_min(i),sector]=min(angle(:,i));
end
distance.sector_antenna_gain=antenna.gain(angle_min);
%% ����·�����
UE.distance_to_enb=sqrt((LTEconfig_add.UE.Position_X-LTEconfig_add.BS.Position_X).^2+(LTEconfig_add.UE.Position_Y-LTEconfig_add.BS.Position_Y).^2);
UE.pathloss_to_enb= LTE_pathloss_enb_to_ue(UE.distance_to_enb);%num_ue*1
%% ������Ӱ˥��
sf=shadowfading(LTEconfig_add.roi_x,LTEconfig_add.roi_y,LTEconfig_add.data_res);%ÿ�������Ӧ����Ӱ˥�� 
for i=1:num_ue
    ue_sf(i,1)=sf(LTEconfig_add.UE.pixel(i,1),LTEconfig_add.UE.pixel(i,2));%num_ue���û�����Ӱ˥�� num_ue*1
end
LTEconfig_add.UE.receiver_noise_figure = 9;    % Receiver noise figure in dB
thermal_noise_W=10^(LTEconfig_add.UE.thermal_noise_density/10) / 1000 * LTEconfig_add.System.Bandwidth* 10^(LTEconfig_add.UE.receiver_noise_figure/10);%ϵͳ�е�������
UE.snr_to_enb=(10.^((LTEconfig_add.BS.Transmit_Power_dB+distance.sector_antenna_gain'-UE.pathloss_to_enb-ue_sf)/10))./(thermal_noise_W+ LTEconfig_add.ue_interference_enbs);%num_ue*1
UE.snr_to_enb_db=10*log10(UE.snr_to_enb);

%% ������ʱ��仯��С�߶�˥���ŵ�
 SimType=8;   % ���÷�������
ltepar=lteparset(SimType,num_ue); % ���÷������
ltepar.Nrb=num_prb; %ϵͳPRB����
ltepar.nUE=num_ue; %�û���Ŀ
 ltepar.Ntot=300;   %ռ�����ز���
ltepar.Bandwidth=5e6; %ϵͳ����
ltepar.Nfft=512; %FFT����
ltepar.Fs=7.68e6; %����Ƶ��
 h=cell(1,ltepar.nUE,num_subframe);
 sc_snr=cell(1,ltepar.nUE);%%Ԫϸ������
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

%% ԭ�е��û��ŵ��������ӵ��û��ŵ�����һ��
LTEconfig.UE.Number_Final= LTEconfig.UE.Number_Final+LTEconfig_add.UE.Number_Final;
LTEconfig.UE.Position_X=cat(1,LTEconfig.UE.Position_X,LTEconfig_add.UE.Position_X);
LTEconfig.UE.Position_Y=cat(1,LTEconfig.UE.Position_Y,LTEconfig_add.UE.Position_Y);
LTEconfig.UE.pixel=cat(1,LTEconfig.UE.pixel,LTEconfig_add.UE.pixel);
LTEconfig.ue_interference_enbs=cat(1,LTEconfig.ue_interference_enbs,LTEconfig_add.ue_interference_enbs);
prb_snr_TTI=cat(2,prb_snr_TTI,prb_snr_TTI_add);
num_ue=LTEconfig.UE.Number_Final;%%%�û���Ŀ

%% ���������ʼ��
sum_throuput_ue_perTTI=zeros(num_ue,num_subframe,num_algorithm);%ÿ��TTI��ÿ���û���������
sum_throuput_ue_until_TTI=zeros(num_ue,num_algorithm);%�����û�����ǰTTI���ۻ�������
sum_throuput_total_all_TTI_RT=zeros(num_algorithm,1);%�����û�������TTI�ڵ�������
sum_throuput_system_kbps_RT=zeros(num_algorithm,1);%ϵͳ������
Jain=zeros(num_algorithm,1);
for iii=1:1:3
        if iii==1
            allocation_type='RR';
        elseif iii==2
            allocation_type='MT';
        elseif  iii==3
            allocation_type='PF';
        end
       allocation_begin_ue=zeros(1,num_subframe)+1;%�����û���ʼ����������RR�㷨��
      CQI=zeros(num_prb,num_ue); %CQI������0
      windows_band=5;%����
      windows_throuput_ue=zeros(1,num_ue); %ÿ���û��Ĵ�������
      past_throuput_ue=zeros(1,num_ue);
      for TTI_ID=1:1:num_subframe
             CQI=SINR_mapping_CQI(prb_snr_TTI(:,:,1:TTI_ID),TTI_ID,num_ue,num_prb);


             %%  �����㷨
             if(strcmp(allocation_type,'RR'))%ѡ�������㷨    strcmp�ǱȽ��㷨
                [allocation_RB,allocation_begin_ue(1,TTI_ID+1)]=downlink_RR(num_prb,num_ue,CQI(:,:),allocation_begin_ue(1,TTI_ID),num_symbol,num_subcarrier); 
            elseif (strcmp(allocation_type,'PF'))
                [allocation_RB]=downlink_PF(num_prb,num_ue,CQI(:,:),windows_throuput_ue,num_symbol,num_subcarrier);   
                
            elseif (strcmp(allocation_type,'MT'))
                [allocation_RB]=downlink_MT(num_prb,num_ue,CQI(:,:),num_symbol,num_subcarrier);
            else error('NOT support this schedule way!');
             end
           
            %% ����ÿ��TTI��ÿ���û���������
            for jj=1:1:num_ue
                 ue_prb=find(allocation_RB(:,3)==jj);
                 num_ue_prb=length(ue_prb);
                 for i=1:1:num_ue_prb
                   sum_throuput_ue_perTTI(jj,TTI_ID,iii)=sum_throuput_ue_perTTI(jj,TTI_ID,iii)+allocation_RB(ue_prb(i),4);
                 end
            end
           
            %% ����ÿ���û�����ǰTTI���ۻ�������
            for jj=1:1:num_ue
               sum_throuput_ue_until_TTI(jj,iii)=sum_throuput_ue_until_TTI(jj,iii)+sum_throuput_ue_perTTI(jj,TTI_ID,iii);
            end
        
       
            %% �����û��Ĵ�������
             for jj=1:1:num_ue    
                if (TTI_ID<windows_band+1)
                    past_throuput_ue(jj)=sum_throuput_ue_until_TTI(jj,iii); %��ǰTTI���㴰������ǰTTI֮ǰ����������Ϊ���û��Ĵ�������
                elseif (TTI_ID>windows_band)
                    past_throuput_ue(jj)=past_throuput_ue(jj)+sum_throuput_ue_perTTI(jj,TTI_ID,iii)-sum_throuput_ue_perTTI(jj,TTI_ID-windows_band,iii);%����TTI���ڴ����Ժ��ÿ���û��Ĵ�������
                end
             windows_throuput_ue(jj)=(1-1/windows_band)*past_throuput_ue(jj)+(1/windows_band)*sum_throuput_ue_perTTI(jj,TTI_ID,iii);
             end   
       
end
         %% ����ʵʱҵ���û��ڵ�ǰ���䷽ʽ�µ�ϵͳ������
         sum_throuput_system_kbps_RT(iii,1)=sum(sum_throuput_ue_until_TTI(:,iii))/num_subframe; 
         %%  ��������ʵʱҵ���û�Jainָ��
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
xlabel('ϵͳ�����ʵʱҵ���û���������');
ylabel('ϵͳ��������Mbps��');
grid on;
hold on;



figure(5);
t=20:10:100;
plot(t,Jain_average(1,:),'-s',t,Jain_average(2,:),'-p',t,Jain_average(3,:),'-o','linewidth',2);
legend('RR','MT','PF');
xlabel('ϵͳ�����ʵʱҵ���û���������');
ylabel('Jainsָ��');
grid on;
hold on;

