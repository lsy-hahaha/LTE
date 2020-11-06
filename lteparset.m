function LtePar=lteparset(simtype,num_ue)
%% 仿真参数配置，被SIM_SISO_TEST调用
%  作者： 赵龙海
% 2012 哈工大通信所
% simtype=1   SUSISO-- Single User Single input and Single output
% simtype=2   MUSISO-- MultiUser Single input and Single output
% simtype=3   SUMIMO-- Single User Single input and Single output
% simtype=4   MUMIMO-- MultiUser  Single input and Single output
% 
% if length(varargin)>1
%     error('No such functionality yet. Try ''ltepar=lteparset( simtype )'' instead.')
% end

LtePar.release='r10';  % LTE version
if simtype==1
    LtePar.SimType='SUSISO'; % simulation type, SUSISO, MUSISO, SUMIMO, MUMIMO
elseif simtype==2
    LtePar.SimType='MUSISO';
elseif simtype==3
    LtePar.SimType='SUSISO_1_4M';     % 
elseif simtype==4
    LtePar.SimType='Fading_Test1x1';   % 1x1 Fading 
elseif simtype==5
    LtePar.SimType='SU_MIMO2x2';       % Open Loop Spatial Multiplexing
elseif simtype==6
    LtePar.SimType='SU_SIMO1x2'   ;   % Rx Diversity
elseif simtype==7
    LtePar.SimType='SU_MIMO4x4';  %  
elseif simtype==8
    LtePar.SimType='MU_SISO1x1';
else
    error('not implemented yet');
end

LtePar.PRB=180e3;     % Physical Resource Block Bandwidth (Hz)
LtePar.SubcarrierSpacing=15e3;  % Subcarrier spacing (Hz)
LtePar.Nsc=12;    % Number of subcarriers of one PRB

%LtePar.Nrb=50;  % Number of PRB of system
LtePar.CarrierFrequency=2;     % GHz
LtePar.SpeedOfLight=3e8;


%% SimType Configuration
%% SUSISO
if strcmp(LtePar.SimType,'SUSISO')
 LtePar.Bandwidth=10e6;  % system bandwidth (Hz)  10MHz   
LtePar.UEpar.mode=1;   %  trasmission model: 1, single antenna; 2, receiver diversity; 3, OLSM; 4, CLSM; 5 MUMIMO;
LtePar.UEpar.nTX=1;   % number of antennas
LtePar.UEpar.speed=1/3.6;  % user speed
LtePar.UEpar.nCodewords=1;
LtePar.UEpar.nLayers=1;

LtePar.BSpar.channel_estimation_method='PERFECT';   %  'PERFECT', 'LS', 'MMSE'
LtePar.BSpar.receiver='MMSE';  % 'MMSE' 'ZF'
LtePar.BSpar.nRX=1;

LtePar.nUE=1;
LtePar.nBS=1;

LtePar.ChanModPar.ChanType='AWGN';   % 'AWGN', 'PDP', 'WINNER2'
LtePar.Scheduler='Round Robin';    % Round Robin , PF , max CI
%% MUSISO
elseif strcmp(LtePar.SimType,'MUSISO')
    LtePar.Bandwidth=10e6;  % system bandwidth (Hz)  10MHz
LtePar.UEpar.mode=1;
LtePar.UEpar.nTX=1;
LtePar.UEpar.speed=1/3.6;
LtePar.UEpar.nCodewords=1;
LtePar.UEpar.nLayers=1;

LtePar.BSpar.channel_estimation_method='PERFECT';
LtePar.BSpar.receiver='MMSE';  % 'MMSE' 'ZF'
LtePar.BSpar.nRX=1;
% 
 LtePar.nUE=num_ue;
ltePar.nBS=1;

LtePar.ChanModPar.ChanType='AWGN';   % 'AWGN', 'PDP', 'WINNER2'
LtePar.Scheduler='Round Robin';    % Round Robin , PF , max CI
%% SUSISO_1_4M
elseif strcmp(LtePar.SimType,'SUSISO_1_4M')
LtePar.Bandwidth=1.4e6;  % system bandwidth (Hz)  10MHz   
LtePar.UEpar.mode=1;   %  trasmission model: 1, single antenna; 2, transmit diversity; 3, OLSM; 4, CLSM; 5 MUMIMO;
LtePar.UEpar.nTX=1;   % number of antennas
LtePar.UEpar.speed=1/3.6;  % user speed
LtePar.UEpar.nCodewords=1;
LtePar.UEpar.nLayers=1;

LtePar.BSpar.channel_estimation_method='PERFECT';   %  'PERFECT', 'LS', 'MMSE'
LtePar.BSpar.receiver='MMSE';  % 'MMSE' 'ZF'
LtePar.BSpar.nRX=1;

 LtePar.nUE=1;
LtePar.nBS=1;

LtePar.ChanModPar.ChanType='AWGN';   % 'AWGN', 'PDP', 'WINNER2'
LtePar.Scheduler='Round Robin';    % Round Robin , PF , max CI
%% Fading Test
elseif strcmp(LtePar.SimType,'Fading_Test1x1')
LtePar.Bandwidth=1.4e6;  % system bandwidth (Hz)  1.4MHz   
LtePar.UEpar.mode=1;   %  trasmission model: 1, single antenna; 2, transmit diversity; 3, OLSM; 4, CLSM; 5 MUMIMO;
LtePar.UEpar.nTX=1;   % number of antennas
LtePar.UEpar.speed=1/3.6;  % user speed

LtePar.UEpar.nCodewords=1;
LtePar.UEpar.nLayers=1;



LtePar.BSpar.channel_estimation_method='PERFECT';   %  'PERFECT', 'LS', 'MMSE'
LtePar.BSpar.receiver='ZF';  % 'MMSE' 'ZF'
LtePar.BSpar.nRX=1;

 LtePar.nUE=1;
LtePar.nBS=1;

LtePar.ChanModPar.ChanType='TU';   % 'AWGN', 'TU','PedB', 'WINNER2'
LtePar.Scheduler='Round Robin';    % Round Robin , PF , max CI
LtePar.ChanModPar.FadeType='Block Fading';   % 'Fast Fading', 'Block Fading'
%LtePar.ChanModPar.CorrType='correlated';
LtePar.ChanModPar.CorrType='independent';

%% SU_MIMO
elseif strcmp(LtePar.SimType,'SU_MIMO2x2')

LtePar.Bandwidth=1.4e6;
LtePar.UEpar.mode=3;   %  trasmission model: 1, single antenna; 2, transmit diversity; 3, OLSM; 4, CLSM; 5 MUMIMO;
LtePar.UEpar.nTX=2;   % number of antennas
LtePar.UEpar.speed=2/3.6;  % user speed

LtePar.UEpar.nCodewords=1;
LtePar.UEpar.nLayers=2;

LtePar.BSpar.channel_estimation_method='PERFECT';   %  'PERFECT', 'LS', 'MMSE'
LtePar.BSpar.receiver='ZF';  % 'MMSE' 'ZF'
LtePar.BSpar.nRX=2;

 LtePar.nUE=1;
LtePar.nBS=1;

LtePar.ChanModPar.ChanType='TU';   %  AWGN TU  PedB PedA 
LtePar.Scheduler='Round Robin';    % Round Robin , PF , max CI

LtePar.ChanModPar.FadeType='Block Fading';   % 'Fast Fading', 'Block Fading'
%LtePar.ChanModPar.CorrType='correlated';
LtePar.ChanModPar.CorrType='independent';
elseif strcmp(LtePar.SimType,'SU_SIMO1x2')
    LtePar.Bandwidth=1.4e6;
    LtePar.UEpar.mode=2;   %  trasmission model: 1, single antenna; 2, receiver diversity; 3, OLSM; 4, CLSM; 5 MUMIMO;
    LtePar.UEpar.nTX=1;   % number of antennas
    LtePar.UEpar.speed=2/3.6;  % user speed

    LtePar.UEpar.nCodewords=1;
    LtePar.UEpar.nLayers=1;

    LtePar.BSpar.channel_estimation_method='PERFECT';   %  'PERFECT', 'LS', 'MMSE'
    LtePar.BSpar.receiver='ZF';  % 'MMSE' 'ZF'
    LtePar.BSpar.nRX=2;
    
   LtePar.nUE=1;
    LtePar.nBS=1;

    LtePar.ChanModPar.ChanType='TU';   %  AWGN TU  PedB PedA 
    LtePar.Scheduler='Round Robin';    % Round Robin , PF , max CI
    LtePar.ChanModPar.FadeType='Block Fading';   % 'Fast Fading', 'Block Fading'
    LtePar.ChanModPar.CorrType='independent';
    
elseif strcmp(LtePar.SimType,'SU_MIMO4x4')
LtePar.Bandwidth=1.4e6;
LtePar.UEpar.mode=3;   %  trasmission model: 1, single antenna; 2, transmit diversity; 3, OLSM; 4, CLSM; 5 MUMIMO;
LtePar.UEpar.nTX=4;   % number of antennas
LtePar.UEpar.speed=2/3.6;  % user speed

LtePar.UEpar.nCodewords=2;
LtePar.UEpar.nLayers=4;

LtePar.BSpar.channel_estimation_method='PERFECT';   %  'PERFECT', 'LS', 'MMSE'
LtePar.BSpar.receiver='ZF';  % 'MMSE' 'ZF'
LtePar.BSpar.nRX=4;

 LtePar.nUE=1;
LtePar.nBS=1;

LtePar.ChanModPar.ChanType='TU';   %  AWGN TU  PedB PedA 
LtePar.Scheduler='Round Robin';    % Round Robin , PF , max CI

LtePar.ChanModPar.FadeType='Block Fading';   % 'Fast Fading', 'Block Fading'
%LtePar.ChanModPar.CorrType='correlated';
LtePar.ChanModPar.CorrType='independent';
elseif strcmp(LtePar.SimType,'MU_SISO1x1')
    LtePar.Bandwidth=5e6;%1.4e6;
LtePar.UEpar.mode=1;   %  trasmission model: 1, single antenna; 2, transmit diversity; 3, OLSM; 4, CLSM; 5 MUMIMO;
LtePar.UEpar.nTX=1;   % number of antennas
LtePar.UEpar.speed=1/3.6;  % user speed

LtePar.UEpar.nCodewords=1;
LtePar.UEpar.nLayers=1;

LtePar.BSpar.channel_estimation_method='PERFECT';   %  'PERFECT', 'LS', 'MMSE'
LtePar.BSpar.receiver='ZF';  % 'MMSE' 'ZF'
LtePar.BSpar.nRX=1;

 LtePar.nUE=num_ue;
LtePar.nBS=1;

LtePar.ChanModPar.ChanType='TU';   %  AWGN TU  PedB PedA 
LtePar.Scheduler='max CI';    % Round Robin , PF , max CI

LtePar.ChanModPar.FadeType='Block Fading';   % 'Fast Fading', 'Block Fading'
%LtePar.ChanModPar.CorrType='correlated';
LtePar.ChanModPar.CorrType='independent';
else    
    error('not implemented yet');
end

if(LtePar.Bandwidth==1.4e6)
    LtePar.Nrb=6;
else
    LtePar.Nrb=(LtePar.Bandwidth*0.9)/LtePar.PRB;
end

LtePar.Ntot=LtePar.Nsc*LtePar.Nrb; % Number of subcarriers of all PRBs
LtePar.HARQ=false;    % false : No HARQ;

LtePar.OFDM_mode=0;   % 0: SC-FDMA  1: OFDM

%% FFT Lengths ,CP lengths 
LtePar.FrameDur=10e-3;   % Frame duration (s)

if(LtePar.Bandwidth==15e6)
    LtePar.Nfft=1536;   % Number of FFT points
else
    LtePar.Nfft=2^ceil(log2(LtePar.Ntot));
end

LtePar.Fs=LtePar.SubcarrierSpacing*LtePar.Nfft;  % Sampling frequency (Hz)
LtePar.CyclicPrefix='normal';
LtePar.Tg=zeros(1,2);    
LtePar.Tg(1)=160/(15000*2048);   % normal CP time - first symbol in subframe
LtePar.Tg(2)=144/(15000*2048);   % normal CP time - remaining symbols in subframe
LtePar.NsymbolSlot=7;
LtePar.NsymbolSubframe=2*7;

LtePar.Ng=zeros(1,2);
LtePar.Ng(1)=LtePar.Tg(1)*LtePar.Fs;
LtePar.Ng(2)=round(LtePar.Tg(2)*LtePar.Fs);
LtePar.Tb=1/LtePar.SubcarrierSpacing;
LtePar.Ts(1)=LtePar.Tb+LtePar.Tg(1);
LtePar.Ts(2)=LtePar.Tb+LtePar.Tg(2);
LtePar.Tslot=LtePar.Tg(1)+LtePar.Tg(2)*(LtePar.NsymbolSlot-1)+LtePar.NsymbolSlot*LtePar.Tb;
LtePar.Tsubframe=2*LtePar.Tslot;
LtePar.TxSymbols=2*(LtePar.Ng(1)+(LtePar.NsymbolSlot-1)*LtePar.Ng(2)+LtePar.Nfft*LtePar.NsymbolSlot);

%% modulation symbol alphabet
LtePar.symbol_alphabet=cell(1,6);
LtePar.symbol_alphabet{1}=[1+1i,-1-1i]/sqrt(2);     % BPSK
LtePar.symbol_alphabet{2}=[1+1i,1-1i,-1+1i,-1-1i]/sqrt(2);        % QPSK
LtePar.symbol_alphabet{4}=[1+1i,1+3i,3+1i,3+3i,1-1i,1-3i,3-1i,3-3i,...
                           -1+1i,-1+3i,-3+1i,-3+3i,-1-1i,-1-3i,-3-1i,-3-3i]/sqrt(10);      % 16QAM
LtePar.symbol_alphabet{6}=[3+3i,3+1i,1+3i,1+1i,3+5i,3+7i,1+5i,1+7i,...
                           5+3i,5+1i,7+3i,7+1i,5+5i,5+7i,7+5i,7+7i,...
                           3-3i,3-1i,1-3i,1-1i,3-5i,3-7i,1-5i,1-7i,...
                           5-3i,5-1i,7-3i,7-1i,5-5i,5-7i,7-5i,7-7i,...
                           -3+3i,-3+1i,-1+3i,-1+1i,-3+5i,-3+7i,-1+5i,-1+7i,...
                           -5+3i,-5+1i,-7+3i,-7+1i,-5+5i,-5+7i,-7+5i,-7+7i,...
                           -3-3i,-3-1i,-1-3i,-1-1i,-3-5i,-3-7i,-1-5i,-1-7i,...
                           -5-3i,-5-1i,-7-3i,-7-1i,-5-5i,-5-7i,-7-5i,-7-7i]/sqrt(42);    % 64QAM
                       
% symbol map to bit table

LtePar.bit_table=cell(1,6);
LtePar.bit_table{1}=[0,1];        % BPSK
LtePar.bit_table{2}=[0,0;0,1;1,0;1,1];       % QPSK
LtePar.bit_table{4}=[0,0,0,0; 0,0,0,1; 0,0,1,0; 0,0,1,1;...          % 16QAM
                     0,1,0,0; 0,1,0,1; 0,1,1,0; 0,1,1,1;...
                     1,0,0,0; 1,0,0,1; 1,0,1,0; 1,0,1,1;...
                     1,1,0,0; 1,1,0,1; 1,1,1,0; 1,1,1,1];
       
LtePar.bit_table{6}=[0,0,0,0,0,0; 0,0,0,0,0,1; 0,0,0,0,1,0; 0,0,0,0,1,1; 0,0,0,1,0,0; 0,0,0,1,0,1; 0,0,0,1,1,0; 0,0,0,1,1,1;...          %64QAM
                     0,0,1,0,0,0; 0,0,1,0,0,1; 0,0,1,0,1,0; 0,0,1,0,1,1; 0,0,1,1,0,0; 0,0,1,1,0,1; 0,0,1,1,1,0; 0,0,1,1,1,1;...
                     0,1,0,0,0,0; 0,1,0,0,0,1; 0,1,0,0,1,0; 0,1,0,0,1,1; 0,1,0,1,0,0; 0,1,0,1,0,1; 0,1,0,1,1,0; 0,1,0,1,1,1;...
                     0,1,1,0,0,0; 0,1,1,0,0,1; 0,1,1,0,1,0; 0,1,1,0,1,1; 0,1,1,1,0,0; 0,1,1,1,0,1; 0,1,1,1,1,0; 0,1,1,1,1,1;...
                     1,0,0,0,0,0; 1,0,0,0,0,1; 1,0,0,0,1,0; 1,0,0,0,1,1; 1,0,0,1,0,0; 1,0,0,1,0,1; 1,0,0,1,1,0; 1,0,0,1,1,1;...
                     1,0,1,0,0,0; 1,0,1,0,0,1; 1,0,1,0,1,0; 1,0,1,0,1,1; 1,0,1,1,0,0; 1,0,1,1,0,1; 1,0,1,1,1,0; 1,0,1,1,1,1;...
                     1,1,0,0,0,0; 1,1,0,0,0,1; 1,1,0,0,1,0; 1,1,0,0,1,1; 1,1,0,1,0,0; 1,1,0,1,0,1; 1,1,0,1,1,0; 1,1,0,1,1,1;...
                     1,1,1,0,0,0; 1,1,1,0,0,1; 1,1,1,0,1,0; 1,1,1,0,1,1; 1,1,1,1,0,0; 1,1,1,1,0,1; 1,1,1,1,1,0; 1,1,1,1,1,1];
                 
                 
                 
 %% 空分复用码本参数  {i,j}为码本索引为 i-1 ，层数为j的码本
% 2天线码本  
LtePar.codebook_pars.two=cell(5,2);
LtePar.codebook_pars.two{1,1}=1/sqrt(2)*[1;1];
LtePar.codebook_pars.two{1,2}=1/sqrt(2)*[1,0;0,1];
LtePar.codebook_pars.two{2,1}=1/sqrt(2)*[1;-1];
LtePar.codebook_pars.two{3,1}=1/sqrt(2)*[1;1i];
LtePar.codebook_pars.two{4,1}=1/sqrt(2)*[1;-1i];
LtePar.codebook_pars.two{5,1}=1/sqrt(2)*[1;0];
LtePar.codebook_pars.two{6,1}=1/sqrt(2)*[0;1];

% 4天线码本   
LtePar.codebook_pars.four=cell(24,4);

% 1 Layers
LtePar.codebook_pars.four{1,1}=1/2*[1;1;1;-1];
LtePar.codebook_pars.four{2,1}=1/2*[1;1;1i;1i];
LtePar.codebook_pars.four{3,1}=1/2*[1;1;-1;1];
LtePar.codebook_pars.four{4,1}=1/2*[1;1;-1i;-1i];
LtePar.codebook_pars.four{5,1}=1/2*[1;1i;1;1i];
LtePar.codebook_pars.four{6,1}=1/2*[1;1i;1i;1];
LtePar.codebook_pars.four{7,1}=1/2*[1;1i;-1;-1i];
LtePar.codebook_pars.four{8,1}=1/2*[1;1i;-1i;-1];
LtePar.codebook_pars.four{9,1}=1/2*[1;-1;1;1];
LtePar.codebook_pars.four{10,1}=1/2*[1;-1;1i;-1i];
LtePar.codebook_pars.four{11,1}=1/2*[1;-1;-1;-1];
LtePar.codebook_pars.four{12,1}=1/2*[1;-1;-1i;1i];
LtePar.codebook_pars.four{13,1}=1/2*[1;-1i;1;-1i];
LtePar.codebook_pars.four{14,1}=1/2*[1;-1i;1i;-1];
LtePar.codebook_pars.four{15,1}=1/2*[1;-1i;-1;1i];
LtePar.codebook_pars.four{16,1}=1/2*[1;-1i;-1i;1];
LtePar.codebook_pars.four{17,1}=1/2*[1;0;1;0];
LtePar.codebook_pars.four{18,1}=1/2*[1;0;-1;0];
LtePar.codebook_pars.four{19,1}=1/2*[1;0;1i;0];
LtePar.codebook_pars.four{20,1}=1/2*[1;0;-1i;0];
LtePar.codebook_pars.four{21,1}=1/2*[0;1;0;1];
LtePar.codebook_pars.four{22,1}=1/2*[0;1;0;-1];
LtePar.codebook_pars.four{23,1}=1/2*[0;1;0;1i];
LtePar.codebook_pars.four{24,1}=1/2*[0;1;0;-1i];

% 2 Layers
LtePar.codebook_pars.four{1,2}=1/2*[1,0;1,0;0,1;0,-1i];
LtePar.codebook_pars.four{2,2}=1/2*[1,0;1,0;0,1;0,1i];
LtePar.codebook_pars.four{3,2}=1/2*[1,0;-1i,0;0,1;0,1];
LtePar.codebook_pars.four{4,2}=1/2*[1,0;-1i,0;0,1;0,-1];
LtePar.codebook_pars.four{5,2}=1/2*[1,0;-1,0;0,1;0,-1i];
LtePar.codebook_pars.four{6,2}=1/2*[1,0;-1,0;0,1;0,1i];
LtePar.codebook_pars.four{7,2}=1/2*[1,0;1i,0;0,1;0,1];
LtePar.codebook_pars.four{8,2}=1/2*[1,0;1i,0;0,1;0,-1];
LtePar.codebook_pars.four{9,2}=1/2*[1,0;0,1;1,0;0,1];
LtePar.codebook_pars.four{10,2}=1/2*[1,0;0,1;1,0;0,-1];
LtePar.codebook_pars.four{11,2}=1/2*[1,0;0,1;-1,0;0,1];
LtePar.codebook_pars.four{12,2}=1/2*[1,0;0,1;-1,0;0,-1];
LtePar.codebook_pars.four{13,2}=1/2*[1,0;0,1;0,1;1,0];
LtePar.codebook_pars.four{14,2}=1/2*[1,0;0,1;0,-1;1,0];
LtePar.codebook_pars.four{15,2}=1/2*[1,0;0,1;0,1;-1,0];
LtePar.codebook_pars.four{16,2}=1/2*[1,0;0,1;0,-1;-1,0];

% 3 Layers
LtePar.codebook_pars.four{1,3}=1/2*[1,0,0;1,0,0;0,1,0;0,0,1];
LtePar.codebook_pars.four{2,3}=1/2*[1,0,0;-1,0,0;0,1,0;0,0,1];
LtePar.codebook_pars.four{3,3}=1/2*[1,0,0;0,1,0;1,0,0;0,0,1];
LtePar.codebook_pars.four{4,3}=1/2*[1,0,0;0,1,0;-1,0,0;0,0,1];
LtePar.codebook_pars.four{5,3}=1/2*[1,0,0;0,1,0;0,0,1;1,0,0];
LtePar.codebook_pars.four{6,3}=1/2*[1,0,0;0,1,0;0,0,1;-1,0,0];
LtePar.codebook_pars.four{7,3}=1/2*[0,1,0;1,0,0;1,0,0;0,0,1];
LtePar.codebook_pars.four{8,3}=1/2*[0,1,0;1,0,0;-1,0,0;0,0,1];
LtePar.codebook_pars.four{9,3}=1/2*[0,1,0;1,0,0;0,0,1;1,0,0];
LtePar.codebook_pars.four{10,3}=1/2*[0,1,0;1,0,0;0,0,1;-1,0,0];
LtePar.codebook_pars.four{11,3}=1/2*[0,1,0;0,0,1;1,0,0;1,0,0];
LtePar.codebook_pars.four{12,3}=1/2*[0,1,0;0,0,1;1,0,0;-1,0,0];

% 4 Layers
LtePar.codebook_pars.four{1,4}=1/2*[1,0,0,0;0,1,0,0;0,0,1,0;0,0,0,1];


%% Channel Model Parameters

LtePar.ChanModPar.winner_settings.PropaCondition='NLOS';
LtePar.ChanModPar.winner_settings.Scenario=11;
LtePar.ChanModPar.winner_settings.SampleDensity=2;
LtePar.ChanModPar.winner_settings.UniformTimeSampling='yes';
LtePar.ChanModPar.winner_settings.FixedPdpUsed='no';
LtePar.ChanModPar.winner_settings.FixedAnglesUsed='no';
LtePar.ChanModPar.winner_settings.PolariseArrays='yes';
LtePar.ChanModPar.winner_settings.ShadowingModelUsed='no';
LtePar.ChanModPar.winner_settings.PathLossModel='pathloss';
LtePar.ChanModPar.winner_settings.PathLossOption='CR_light';
LtePar.ChanModPar.winner_settings.RandomSeed=[];
LtePar.ChanModPar.winner_settings.UseManualPropCondition='yes';

if strcmp(LtePar.ChanModPar.ChanType,'TU')
    LtePar.ChanModPar.PDP_dB=[-5.7000 -7.6000 -10.1000 -10.2000 -10.2000 -11.5000 -13.4000 -16.3000 -16.9000 -17.1000 -17.4000,...
                          -19.0000 -19.0000 -19.8000 -21.5000 -21.6000 -22.1000 -22.6000 -23.5000 -24.3000; % Average power [dB]           % Typical Urban
                          0 0.2170 0.5120 0.5140 0.5170 0.6740 0.8820 1.2300 1.2870 1.3110 1.3490 1.5330 1.5350,...
                          1.6220 1.8180 1.8360 1.8840 1.9430 2.0480 2.1400];
    LtePar.ChanModPar.PDP_dB(2,:)=LtePar.ChanModPar.PDP_dB(2,:)*10^-6;
elseif strcmp(LtePar.ChanModPar.ChanType,'PedA')
    LtePar.ChanModPar.PDP_dB=[0 -9.7 -19.2 -22.8;
                            0 110*10^-9 190*10^-9 410*10^-9];
elseif strcmp(LtePar.ChanModPar.ChanType,'PedB')
    LtePar.ChanModPar.PDP_dB=[0 -0.9 -4.9 -8 -7.8 -23.9;
                            0 200*10^-9 800*10^-9 1200*10^-9 2300*10^-9 3700*10^-9];
elseif strcmp(LtePar.ChanModPar.ChanType,'VehA')
    LtePar.ChanModPar.PDP_dB=[0 -1 -9 -10 -15 -20;
                            0 310*10^-9 710*10^-9 1090*10^-9 1730*10^-9 2510*10^-9];
elseif strcmp(LtePar.ChanModPar.ChanType,'AWGN')
    ;
else
    error('this channel is not implemented yet!');
end


LtePar.ChanModPar.M=8;
if ~strcmp(LtePar.ChanModPar.ChanType,'AWGN')
    LtePar.psi=-pi+2*pi*rand(LtePar.BSpar.nRX,LtePar.UEpar.nTX,size(LtePar.ChanModPar.PDP_dB,2),LtePar.ChanModPar.M,LtePar.nUE);
    LtePar.theta=-pi+2*pi*rand(LtePar.BSpar.nRX,LtePar.UEpar.nTX,size(LtePar.ChanModPar.PDP_dB,2),LtePar.ChanModPar.M,LtePar.nUE);
    LtePar.phi=-pi+2*pi*rand(LtePar.BSpar.nRX,LtePar.UEpar.nTX,size(LtePar.ChanModPar.PDP_dB,2),LtePar.ChanModPar.M,LtePar.nUE);
end

if LtePar.UEpar.nTX>1
    LtePar.ChanModPar.TxCorrMatrix=eye(LtePar.UEpar.nTX);
else
    LtePar.ChanModPar.TxCorrMatrix=1;
end

if LtePar.BSpar.nRX>1
    LtePar.ChanModPar.RxCorrMatrix=eye(LtePar.BSpar.nRX);
else
    LtePar.ChanModPar.RxCorrMatrix=1;
end


%% CQI parameters
%LtePar.CQIpar(1).CQI=1;
%LtePar.CQIpar(1).modulation='QPSK';
LtePar.CQIpar(1).modulation_order=2;
LtePar.CQIpar(1).code_rate_x_1024=78;
LtePar.CQIpar(1).efficiency=0.1523;

LtePar.CQIpar(2).modulation_order=2;
LtePar.CQIpar(2).code_rate_x_1024=120;
LtePar.CQIpar(2).efficiency=0.2344;

LtePar.CQIpar(3).modulation_order=2;
LtePar.CQIpar(3).code_rate_x_1024=193;
LtePar.CQIpar(3).efficiency=0.3770;

LtePar.CQIpar(4).modulation_order=2;
LtePar.CQIpar(4).code_rate_x_1024=308;
LtePar.CQIpar(4).efficiency=0.6016;

LtePar.CQIpar(5).modulation_order=2;
LtePar.CQIpar(5).code_rate_x_1024=449;
LtePar.CQIpar(5).efficiency=0.8770;

LtePar.CQIpar(6).modulation_order=2;
LtePar.CQIpar(6).code_rate_x_1024=602;
LtePar.CQIpar(6).efficiency=1.1758;

LtePar.CQIpar(7).modulation_order=4;
LtePar.CQIpar(7).code_rate_x_1024=378;
LtePar.CQIpar(7).efficiency=1.4766;

LtePar.CQIpar(8).modulation_order=4;
LtePar.CQIpar(8).code_rate_x_1024=490;
LtePar.CQIpar(8).efficiency=1.9141;

LtePar.CQIpar(9).modulation_order=4;
LtePar.CQIpar(9).code_rate_x_1024=616;
LtePar.CQIpar(9).efficiency=2.4063;

LtePar.CQIpar(10).modulation_order=6;
LtePar.CQIpar(10).code_rate_x_1024=466;
LtePar.CQIpar(10).efficiency=2.7305;

LtePar.CQIpar(11).modulation_order=6;
LtePar.CQIpar(11).code_rate_x_1024=567;
LtePar.CQIpar(11).efficiency=3.3223;

LtePar.CQIpar(12).modulation_order=6;
LtePar.CQIpar(12).code_rate_x_1024=666;
LtePar.CQIpar(12).efficiency=3.9023;

LtePar.CQIpar(13).modulation_order=6;
LtePar.CQIpar(13).code_rate_x_1024=772;
LtePar.CQIpar(13).efficiency=4.5234;

LtePar.CQIpar(14).modulation_order=6;
LtePar.CQIpar(14).code_rate_x_1024=873;
LtePar.CQIpar(14).efficiency=5.1152;

LtePar.CQIpar(15).modulation_order=6;
LtePar.CQIpar(15).code_rate_x_1024=948;
LtePar.CQIpar(15).efficiency=5.5547;

%% beta for EESM
LtePar.beta=[0.18 0.5];

%% channel autocorrelation matrix
% LtePar.ChanModPar.normH=sqrt(sum(10.^(LtePar.ChanModPar.PDP_dB(1,:)/10)));
% NTap=size(LtePar.ChanModPar.PDP_dB,2);
% tap_delays=round(LtePar.ChanModPar.PDP_dB(2,:)*15000*1024);
% LtePar.ChanModPar.H=zeros(tap_delays(end)+1,1);
% for tap_i=1:NTap
%     if(tap_i==1)
%         prev_delay=-1;
%     else
%         prev_delay=tap_delays(tap_i-1);
%     end
%     if(tap_delays(tap_i)~=prev_delay)
%         H(tap_delays(tap_i)+1)=sqrt(10.^(LtePar.ChanModPar.PDP_dB(1,tap_i)./10));
%     end
% end


% LtePar.ChanModPar.H=H./LtePar.ChanModPar.normH;
% DFT_matrix=dftmtx(LtePar.Nfft);
% matrix_tmp=zeros(LtePar.Nfft);
% matrix_tmp(1:length(LtePar.ChanModPar.H),1:length(LtePar.ChanModPar.H))=diag(LtePar.ChanModPar.H.^2);
end







