%This script calculate: PSD, cross-PSD and Coherence in two LFP 
%during events of object exploration.
%It use coherencyc from chronux

close all
clear
clc

%% Parameters
%paramestros para definir segmentos a analizar
z = 2550; %extention of the segment to analyse
y = 0; %how much before the start I want to include
criterio = 400; %criteria, minimal length of an events to be included (samples)
point = 0; %1  if I want to use the median, 0 if I want the previous criteria

%parameters for multi-tapper fft
W = 1; %window size
T = 1; %temporal window fft multi-taper
p = 1; %Tapers=2*TW-p
fmx = 40; %maximal frequency
fm = 1250; %sampling frequency
dt = 1/fm;
fnyq = fm/5; %Nyquist

%current directory
cd ('/home/facu/Escritorio/Rata3/Analisis') %current folder

%% Open Channel from the HPC
[FileName,PathName,FilterIndex] = uigetfile('*.txt','Seleccionar canal del HIPOCAMPO');
lfpHP = load([PathName,FileName]);
lfpHP = detrend(lfpHP,0); %eliminate drift

%% Open Channel from the PFC
[FileName1,PathName1,FilterIndex1] = uigetfile('*.txt','Seleccionar canal de la PFC',PathName);
lfpPFC=load([PathName1,FileName1]);
lfpPFC=detrend(lfpPFC,0); %eliminate drift

%% Notch in 50 HZ
Wo = 50/fnyq; Bw=Wo/20;
[b,a] = iirnotch(Wo,Bw,1);  

lfpHP = filtfilt(b,a,lfpHP);
lfpPFC = filtfilt(b,a,lfpPFC);

%% Filtering the signals
%parametros del filtro y SF

LowCut=2/(fm/2);  % low-pass in 2Hz
HighCut=(fnyq-1)/(fm/2);  %high-pass in Nyquist (624 Hz)

%butterworth filter
[B,A]=butter(3,[LowCut, HighCut]);

%Filter in the HPC signal
lfpHP = filtfilt(B,A,lfpHP(:,1));

%Filter in the PFC signal
lfpPFC = filtfilt(B,A,lfpPFC(:,1));

%% Load the Events.mat file
fm_original = 32556; %frequency before downsampling
[FileName2,PathName2,FilterIndex2] = uigetfile('*.mat','Seleccionar Events.mat',PathName);
load ([PathName2,FileName2],'TimeStamps','interval_rojo','interval_verde','ID_separada');
size_T = length(lfpHP);
time=((TimeStamps(1):((1/1250)*1000000):(((size_T/1250)*1000000)+TimeStamps(1)-((1/1250)*1000000)))); %time vector starting from the start recording time
time=time'; %transpose to column

%% Find the position in time vector of the ON(Start) and OFF(end) of events
%time_rojo = closest time value from time vector in comparison to the events
%index_rojo = position of time_rojo in the time vector
%interval_rojo = ON and OFF of all the events
%Timestamps from Neuralynx Digital Inputs
timepoint_rojo = interval_rojo(:,1); %for ON
timepoint_rojo_off = interval_rojo(:,2); %for OFF
timepoint_verde = interval_verde(:,1);
timepoint_verde_off = interval_verde(:,2);

%This loops looks for the position in the time vector of the events
index_rojo=zeros(length(timepoint_rojo),1); 
time_rojo=zeros(length(timepoint_rojo),1);
index_rojo_off=zeros(length(timepoint_rojo_off),1);
time_rojo_off=zeros(length(timepoint_rojo_off),1);
%for incongruent object events
for i=1:length(timepoint_rojo)
    [time_rojo(i),index_rojo(i)] = min(abs(time-timepoint_rojo(i)));
    [time_rojo_off(i),index_rojo_off(i)] = min(abs(time-timepoint_rojo_off(i)));
end


index_verde=zeros(length(timepoint_verde),1);
time_verde=zeros(length(timepoint_verde),1);
index_verde_off=zeros(length(timepoint_verde_off),1);
time_verde_off=zeros(length(timepoint_verde_off),1);
%for congruent object evnts
for i=1:length(timepoint_verde)
    [time_verde(i),index_verde(i)] = min(abs(time-timepoint_verde(i)));
    [time_verde_off(i),index_verde_off(i)] = min(abs(time-timepoint_verde_off(i)));
end

idx_verde=[index_verde,index_verde_off]; %ON and OFF Congruent
idx_rojo=[index_rojo,index_rojo_off]; %ON and OFF Incongruent

clear time_rojo time_rojo_off time_verde time_verde_off
clear timepoint_rojo timepoint_rojo_off timepoint_verde timepoint_verde_off
clear index_verde index_verde_off index_rojo index_rojo_off interval_rojo interval_verde
clear LowCut HighCut FilterIndex  FilterIndex1 FilterIndex2 fmd

%% Generates the Non-Object exploration events
% I take +/-500ms from the ON and OFF of the events, to be sure.
V(:,1) = idx_verde(:,1)-625;
V(:,2) = idx_verde(:,2)+625;
R(:,1) = idx_rojo(:,1)-625;
R(:,2) = idx_rojo(:,2)+625;
idx_no = [R;V];
clear V R

HP_NO = lfpHP; %Generates LFP for Non-object exploration
PFC_NO = lfpPFC;

for t = 1:length(idx_no) %replace segments for NaN
HP_NO(idx_no(t,1):idx_no(t,2),1) = NaN;
PFC_NO(idx_no(t,1):idx_no(t,2),1) = NaN;
end

clear t idx_no

HP_NO=HP_NO(~isnan(HP_NO)); %now, eliminates NaN
PFC_NO=PFC_NO(~isnan(PFC_NO)); 

% Generates index to restrict NO events
[time_inicio,idx_inicio] = min(abs(time-TimeStamps(1,2))); %busca el indice del primer TTL
clear time_inicio
L=fix((length(HP_NO)-idx_inicio)/1250/2);
idx_no=[];
for i=1:L %this loop take 0.5 seeach 1 sec
    idx_no(i,1) = idx_inicio+((i-1)*1250);
    idx_no(i,2) = idx_no(i,1)+625;
end
clear i L

%shuffle the index to take segments from different part of the session
p = randperm (length(idx_no));
idx_no_tmp = [];
for c = 1:length(idx_no)
idx_no_tmp(c,:) = idx_no(p(1,c),:);
end

idx_no = idx_no_tmp;
clear c idx_no_tmp p

%concatenates Non-object segments
no_HP=[];
no_PFC=[];
for ii = 1:length(idx_no)
    %Hippocañpus
    ev_no_HP_tmp=HP_NO(idx_no(ii,1):idx_no(ii,2),1);
    H = hamming(length(ev_no_HP_tmp),'symmetric');
    ev_no_HP_tmp=ev_no_HP_tmp.*H;
    no_HP=[no_HP;ev_no_HP_tmp];
    clear ev_no_HP_tmp H

    %Prefrontal Cortex
    ev_no_PFC_tmp=PFC_NO(idx_no(ii,1):idx_no(ii,2),1);
    H = hamming(length(ev_no_PFC_tmp),'symmetric');
    ev_no_PFC_tmp=ev_no_PFC_tmp.*H;
    no_PFC=[no_PFC;ev_no_PFC_tmp];
    clear ev_no_PFC_tmp H    
end

%% Baseline
baseline_HP = lfpHP_filt(1:idx_inicio,1);
baseline_PFC = lfpPFC_filt(1:idx_inicio,1);

L = fix(idx_inicio/1250/2);

idx_baseline = [];

for i=1:L %this loop take 0.5 seeach 1 sec
    idx_baseline(i,1) = i*1250;
    idx_baseline(i,2) = idx_baseline(i,1)+625;
end
clear i L

%Suffle Baseline Events
p = randperm (length(idx_baseline));
idx_baseline_tmp = [];
for c = 1:length(idx_baseline)
idx_baseline_tmp(c,:) = idx_baseline(p(1,c),:);
end

idx_baseline=idx_baseline_tmp;
clear c idx_baseline_tmp p

%Concatenate Baseline segments
base_HP=[];
base_PFC=[];
for ii = 1:length(idx_baseline)
    %Hippocampus
    ev_base_HP_tmp=baseline_HP(idx_baseline(ii,1):idx_baseline(ii,2),1);
    H = hamming(length(ev_base_HP_tmp),'symmetric');
    ev_base_HP_tmp=ev_base_HP_tmp.*H;
    base_HP=[base_HP;ev_base_HP_tmp];
    clear ev_base_HP_tmp H

    %Prefrontal
    ev_base_PFC_tmp=baseline_PFC(idx_baseline(ii,1):idx_baseline(ii,2),1);
    H = hamming(length(ev_base_PFC_tmp),'symmetric');
    ev_base_PFC_tmp=ev_base_PFC_tmp.*H;
    base_PFC=[base_PFC;ev_base_PFC_tmp];
    clear ev_base_PFC_tmp H    
end

%% Discard of segments according criteria
%Median calculation
median_V=median((idx_verde(:,2)-idx_verde(:,1)),1);
median_R=median((idx_rojo(:,2)-idx_rojo(:,1)),1);

if median_V < median_R
    median = median_V;
else
    median = median_R;
end

if point == 1
    criterio = median;
else
    criterio = criterio;
end

clear median_V median_R

%loop congruent events
idx_verde_tmp = [ ];
for c = 1 : length(idx_verde)
    if(idx_verde(c,2)-idx_verde(c,1)) > criterio
    idx_verde_tmp=[idx_verde_tmp; idx_verde(c,:)]; 
    end
end
idx_verde = idx_verde_tmp;
clear c idx_verde_tmp
%suffle of congruent events
p = randperm (length(idx_verde));
idx_verde_tmp = [];
for c = 1:length(idx_verde)
idx_verde_tmp(c,:) = idx_verde(p(1,c),:);
end
idx_verde=idx_verde_tmp;
clear c idx_verde_tmp p

%loop incongruent events
idx_rojo_tmp=[];
for c=1:length(idx_rojo)
    if (idx_rojo(c,2)-idx_rojo(c,1))>criterio
        idx_rojo_tmp=[idx_rojo_tmp;idx_rojo(c,:)];
    end
end
idx_rojo=idx_rojo_tmp;
clear c idx_rojo_tmp
%suffle of incongruent events
p = randperm (length(idx_rojo));
idx_rojo_tmp = [];
for c = 1:length(idx_rojo)
idx_rojo_tmp(c,:) = idx_rojo(p(1,c),:);
end

idx_rojo = idx_rojo_tmp;
clear c idx_rojo_tmp p

%% Concatenation of Object exploration events
%Congruent events
verde_HP=[];
verde_PFC=[]; 
for ii=1:(length(idx_verde));
    %Hippocamppus
    ev_verde_HP_tmp = lfpHP(idx_verde(ii,1):idx_verde(ii,2),1); %temporal
    H = hamming(length(ev_verde_HP_tmp),'symmetric'); %Hamming window
    ev_verde_HP_tmp = ev_verde_HP_tmp.*H; %applying of Hamming
    verde_HP = [verde_HP;ev_verde_HP_tmp]; %concatenation
    clear ev_verde_HP_tmp H

    %Prefrontal
    ev_verde_PFC_tmp = lfpPFC(idx_verde(ii,1):idx_verde(ii,2),1);
    H = hamming(length(ev_verde_PFC_tmp),'symmetric');
    ev_verde_PFC_tmp = ev_verde_PFC_tmp.*H;
    verde_PFC = [verde_PFC;ev_verde_PFC_tmp];
clear ev_verde_PFC_tmp H
end

%Incongruent events
rojo_HP=[];
rojo_PFC=[];
for ii=1:(length(idx_rojo));
    %Hippocampus
    ev_rojo_HP_tmp=lfpHP_filt(((idx_rojo(ii,1)):(idx_rojo(ii,2))),1);
    H = hamming(length(ev_rojo_HP_tmp),'symmetric');
    ev_rojo_HP_tmp=ev_rojo_HP_tmp.*H;
    rojo_HP=[rojo_HP;ev_rojo_HP_tmp];
    clear ev_rojo_HP_tmp H

    %Prefrontal
    ev_rojo_PFC_tmp=lfpPFC_filt(((idx_rojo(ii,1)):(idx_rojo(ii,2))),1);
    H = hamming(length(ev_rojo_PFC_tmp),'symmetric');
    ev_rojo_PFC_tmp=ev_rojo_PFC_tmp.*H;
    rojo_PFC=[rojo_PFC;ev_rojo_PFC_tmp];
    clear ev_rojo_PFC_tmp H
end


%% Made LFP segments the same size

T = [];
if size(rojo_HP,1) > size(verde_HP,1)
    T = size(verde_HP,1);
else
    T = size(rojo_HP,1);
end

rojo_HP = rojo_HP(1:T,1);
rojo_PFC = rojo_PFC(1:T,1);
verde_HP = verde_HP(1:T,1);
verde_PFC = verde_PFC(1:T,1);
no_PFC = no_PFC(1:T,1);
no_HP = no_HP(1:T,1);
base_HP = base_HP(1:T,1);
base_PFC = base_PFC(1:T,1);
clear PFC_NO HP_NO baseline_HP baseline_PFC lfpPFC lfpHP

%% Normalization of signal to max or min
%Congruent HP
max_ev_verde_HP = max(verde_HP);
min_ev_verde_HP = abs(min(verde_HP));

if max_ev_verde_HP>min_ev_verde_HP 
    verde_HP=verde_HP/max_ev_verde_HP; 
else
    verde_HP=verde_HP/min_ev_verde_HP;
end

%Congruent PFC
max_ev_verde_PFC=max(verde_PFC);
min_ev_verde_PFC=abs(min(verde_PFC));

if max_ev_verde_PFC>min_ev_verde_PFC
    verde_PFC=verde_PFC/max_ev_verde_PFC; 
else
    verde_PFC=verde_PFC/min_ev_verde_PFC;
end

clear max_ev_verde_HP min_ev_verde_HP max_ev_verde_PFC min_ev_verde_PFC ii

%Incongruent HPC
max_ev_rojo_HP=max(rojo_HP);
min_ev_rojo_HP=abs(min(rojo_HP));

if max_ev_rojo_HP>min_ev_rojo_HP
    rojo_HP=rojo_HP/max_ev_rojo_HP; 
else
    rojo_HP=rojo_HP/min_ev_rojo_HP;
end

%Inongruent PFC
max_ev_rojo_PFC=max(rojo_PFC);
min_ev_rojo_PFC=abs(min(rojo_PFC));

if max_ev_rojo_PFC>min_ev_rojo_PFC
    rojo_PFC=rojo_PFC/max_ev_rojo_PFC; 
else
    rojo_PFC=rojo_PFC/min_ev_rojo_PFC;
end

clear max_ev_rojo_PFC min_ev_rojo_PFC max_ev_rojo_HP min_ev_rojo_HP ii

%No-Objects PFC
max_ev_no_PFC=max(no_PFC);
min_ev_no_PFC=abs(min(no_PFC));

if max_ev_no_PFC>min_ev_no_PFC
    no_PFC=no_PFC/max_ev_no_PFC; 
else
  no_PFC=no_PFC/min_ev_no_PFC;
end

%No-Objects HPC
max_ev_no_HP=max(no_HP);
min_ev_no_HP=abs(min(no_HP));

if max_ev_no_HP>min_ev_no_HP
    no_HP=no_HP/max_ev_no_HP; 
else
  no_HP=no_HP/min_ev_no_HP;
end

clear max_no_rojo_PFC min_no_rojo_PFC max_no_rojo_HP min_ev_no_HP i PFC_NO HP_NO

% Baseline HPC
max_ev_base_HP=max(base_HP);
min_ev_base_HP=abs(min(base_HP)); 

if max_ev_base_HP>min_ev_base_HP 
    base_HP=base_HP/max_ev_base_HP; 
else
    base_HP=base_HP/min_ev_base_HP;
end

% Baseline PFC
max_ev_base_PFC=max(base_PFC);
min_ev_base_PFC=abs(min(base_PFC));

if max_ev_base_PFC>min_ev_base_PFC
    base_PFC=base_PFC/max_ev_base_PFC; 
else
    base_PFC=base_PFC/min_ev_base_PFC;
end

clear max_ev_base_HP min_ev_base_HP max_ev_base_PFC min_ev_base_PFC ii

%% Coherence (C), Phase Coherence (phi), Cross-PSD (S12), PSD (S1 y S2)
W = 1;
p = fix((T/1250)/4)-(2*W*(T/1250));

params = struct('tapers',[W T/1250 p],'pad',-1,'Fs',1250,'fpass',[0 40],'err',[2 2],'trialave',1);

[C_red,phi_red,S12_red,S1_red,S2_red,f_red,confC_red,phistd_red,Cerr_red] = coherencyc(rojo_HP,rojo_PFC,params);
[C_green,phi_green,S12_green,S1_green,S2_green,f_green,confC_green,phistd_green,Cerr_green] = coherencyc(verde_HP,verde_PFC,params);
[C_no,phi_no,S12_no,S1_no,S2_no,f_no,confC_no,phistd_no,Cerr_no] = coherencyc(no_HP,no_PFC,params);
[C_base,phi_base,S12_base,S1_base,S2_base,f_base,confC_base,phistd_base,Cerr_base] = coherencyc(base_HP,base_PFC,params);

%% Save
 save(([PathName,'\data.mat']),'C_base','phi_base','S12_base','S1_base','S2_base','f_base','confC_base','phi_green','phi_red','C_green','C_green_S','C_red','C_red_S','S1_green','S1_red','S2_green','S2_red','S1_green_S','S1_red_S','S2_green_S','S2_red_S','S12_green_S','S12_red_S','S12_green','S12_red','f_red','f_green','C_no','phi_no','S12_no','S1_no','S2_no','f_no')
