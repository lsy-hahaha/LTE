function [chan h] = channelmodel(pdp,fs,central_frequency,user_speed,fadetype,corrtype,M,psi,theta,phi,time_i,nTX,nRX,TxCorrMatrix,RxCorrMatrix)
%% 产生衰落信道，发射信号滤波， 被SIM_SISO_TEST调用
% % 作者： 赵龙海
% 2012 哈工大通信所
% 输入： data   需要滤波的SC_FDMA信号
%        pdp    功率延迟分布 
%        fs     系统采样频率 ： 15k x Nfft
%        central_frequency   系统载频（Hz）
%        user_speed   用户速度（m/s）
% 输出： chan   产生的信道对象
%         h    1xL  用于理想信道估计下的均衡
%        outdata    滤波后的信号
% filter the data
% References: 1.Yahong Rosa Zheng, et al. " SImulation Models With Correct Statistical Properties for 
%               Rayleigh Fading Channels", IEEE Trans. On Communications, Vol.51, No.6,  June 2003
%             2.Jean Philippe Kermoal ,et al. " A Stochastic MIMO Radio Channel Model With Experimental
%               Validation", IEEE Journal on Selected Areas in Communications, Vol.20, No.6, August 2002
relative_power_dB=pdp(1,:);
delays_in_samples=round(pdp(2,:).*fs);
distinct_tap_samples=unique(delays_in_samples);
num_faders=length(distinct_tap_samples);

for i=1:num_faders
equal_delay_samples=distinct_tap_samples(i)==delays_in_samples;
if(length(relative_power_dB(equal_delay_samples))==1)
    power(i)=relative_power_dB(equal_delay_samples);
else
    power(i)=sum(10.^(relative_power_dB(equal_delay_samples)/10));%合并相同采样时刻的多径功率
    power(i)=10*log10(power(i));
end
end
speed=user_speed*1000/3600;
fd=speed*central_frequency/3e8;
if strcmp(fadetype,'Fast Fading')
    chan=mimochan(nTX,nRX,1/fs,fd,distinct_tap_samples/fs,power);
    chan.StoreHistory=1;
    chan.ResetBeforeFiltering=0;
    chan.StorePathGains=1;
    chan.TxCorrelationMatrix=TxCorrMatrix;
    chan.RxCorrelationMatrix=RxCorrMatrix;
    outdata=filter(chan,data');
    channel=chan.PathGains(500,:,:,:);
    h=zeros(nRX,nTX,distinct_tap_samples(end)+1);
    for rx_i=1:nRX
        for tx_i=1:nTX
               for k=1:size(channel,2)
                   index=distinct_tap_samples(k);
                   h(rx_i,tx_i,index+1)=channel(1,k,tx_i,rx_i);
               end
        end
    end            
elseif strcmp(fadetype,'Block Fading')   % 'Block Fading' 一个TTI内不变
    power=10.^(power/10);
 %   sum(power)
    path_num=length(power);
    X=zeros(nRX,nTX,path_num);
    chan=[];
    if strcmp(corrtype,'correlated')
        xc=zeros(nRX,nTX,path_num,M);
        xs=zeros(nRX,nTX,path_num,M);
        Xc=zeros(nRX,nTX,path_num);
        Xs=zeros(nRX,nTX,path_num);
       
        for rx_i=1:nRX
            for tx_i=1:nTX
                for path_i=1:path_num
                    for subpath_i=1:M
                     
                      xc(rx_i,tx_i,path_i,subpath_i)=cos(psi(rx_i,tx_i,path_i,subpath_i))*cos(2*pi*fd*time_i*cos((2*pi*subpath_i-pi+theta(rx_i,tx_i,path_i,subpath_i))/(4*M))+phi(rx_i,tx_i,path_i,subpath_i));
                      xs(rx_i,tx_i,path_i,subpath_i)=sin(psi(rx_i,tx_i,path_i,subpath_i))*cos(2*pi*fd*time_i*cos((2*pi*subpath_i-pi+theta(rx_i,tx_i,path_i,subpath_i))/(4*M))+phi(rx_i,tx_i,path_i,subpath_i));
                    end
                  Xc(rx_i,tx_i,path_i)=[2/sqrt(M)]*sum(xc(rx_i,tx_i,path_i,:));
              
                  Xs(rx_i,tx_i,path_i)=[2/sqrt(M)]*sum(xs(rx_i,tx_i,path_i,:));
                 
                  X(rx_i,tx_i,path_i)=sqrt(power(path_i)/2)*(Xc(rx_i,tx_i,path_i)+1i*Xs(rx_i,tx_i,path_i));
                end
            end
        end
                        
    elseif strcmp(corrtype,'independent')
        for rx_i=1:nRX
            for tx_i=1:nTX
                for path_i=1:path_num
                    X(rx_i,tx_i,path_i)=sqrt(power(path_i)/2)*(randn+1i*randn);
                end
            end
        end      
    end
      if ~isequal(TxCorrMatrix,eye(nTX))||~isequal(RxCorrMatrix,eye(nRX))
            disp('antennas are correlated!');
            Rmimo=kron(TxCorrMatrix,RxCorrMatrix);
            C=chol(Rmimo,'lower');
            A=zeros(nRX*nTX,path_num);
            CorrX=zeros(nRX,nTX,path_num);
            for path_i=1:path_num
                A(:,path_i)=C*reshape(X(:,path_i),nRX*nTX,1);
                CorrX(:,:,path_i)=reshape(A(:,path_i),nRX,nTX);
            end
      else
          CorrX=X;%非相关
      %    disp('asdfasfd');
     %     CorrX
      end 
    h=zeros(nRX,nTX,distinct_tap_samples(end)+1);
    for rx_i=1:nRX
        for tx_i=1:nTX
                for k=1:path_num
                     index=distinct_tap_samples(k);
                     h(rx_i,tx_i,index+1)=CorrX(rx_i,tx_i,k);
                end
        end
    end
  %  h
   % CorrX
%   sum(abs(h).^2)
% for rx_i=1:nRX
%     tmp=zeros(1,length(data));
% 
%     for tx_i=1:nTX
%    %     squeeze(h(rx_i,tx_i,:))
%     tmp2=filter(squeeze(h(rx_i,tx_i,:)),1,data(tx_i,:));
% 
%     tmp=tmp2+tmp;
%     end
%     outdata(rx_i,:)=tmp;
% %    chan=[chan;h];
% end

else
    error('no such channel type!');
end

end

