%---------------
% This script take two LFP signals, apply a filter in theta band (6-10 Hz),
% extract the amplitude, using Hilbert transform and cross-correlate.
%---------------
clear
clc
close all

%% Parameters
cd ('myPath') %current folder
dt = 1/1250;
fm = 1250; %sampling frequency
fnyq = fm/2; %Nyquist

x = 1250;
y = 1250;

%Inclusion criteria
criterio = 350; %in samples
point = 1; %1 if I want events' median as a criteria, 0 if I want the previosly defined

%Cross-correlation Lags
lag = round(fm/2); %change the denominator to define el lag (fm/X)

%High and Lowpass frequencies
low = 6;
high = 10;

%% Upload LFPs
%Hippocampus
[FileName,PathName,FilterIndex] = uigetfile('*.txt','Seleccionar canal del HIPOCAMPO');
lfpHP1=load([PathName,FileName]);
lfpHP=detrend(lfpHP1,0); %saca el drift
clear lfpHP1

%Prefrontal Cortex
[FileName1,PathName1,FilterIndex1] = uigetfile('*.txt','Seleccionar canal de la PFC',PathName);
lfpPFC1=load([PathName1,FileName1]);
lfpPFC=detrend(lfpPFC1,0);
clear lfpPFC2

%% Filters application
LowCut=low/(fm/2);  %highpass
HighCut=high/(fm/2);  %lowpass
%butterworth
[B,A]=butter(3,[LowCut, HighCut]);

%Hippocampus
lfpHP_filt= filtfilt(B,A,lfpHP(:,1));
%Prefrontal Cortex
lfpPFC_filt= filtfilt(B,A,lfpPFC(:,1));

%Rename
HP = lfpHP_filt;
PFC = lfpPFC_filt;
clear lfpHP_filt lfpPFC_filt low high LowCut HighCut B A

%% Upload Events.mat
[FileName2,PathName2,FilterIndex2] = uigetfile('*.mat','Seleccionar Events.mat',PathName);
load ([PathName2,FileName2],'TimeStamps','interval_rojo','interval_verde','ID_separada');
size_T=length(lfpHP);
time =((TimeStamps(1):((1/1250)*1000000):(((size_T/1250)*1000000)+TimeStamps(1)-((1/1250)*1000000)))); %time vector
time = time'; 
%% Find Position of events in time vectors
timepoint_rojo=interval_rojo(:,1);
timepoint_rojo_off=interval_rojo(:,2);
timepoint_verde=interval_verde(:,1);
timepoint_verde_off=interval_verde(:,2);
%Incongruent
index_rojo=zeros(length(timepoint_rojo),1); %store of position in time
time_rojo=zeros(length(timepoint_rojo),1);
index_rojo_off=zeros(length(timepoint_rojo_off),1);
time_rojo_off=zeros(length(timepoint_rojo_off),1);
for i=1:length(timepoint_rojo)
    [time_rojo(i),index_rojo(i)] = min(abs(time-timepoint_rojo(i)));
    [time_rojo_off(i),index_rojo_off(i)] = min(abs(time-timepoint_rojo_off(i)));
end
%Congruent
index_verde=zeros(length(timepoint_verde),1);
time_verde=zeros(length(timepoint_verde),1);
index_verde_off=zeros(length(timepoint_verde_off),1);
time_verde_off=zeros(length(timepoint_verde_off),1);
for i=1:length(timepoint_verde)
    [time_verde(i),index_verde(i)] = min(abs(time-timepoint_verde(i)));
    [time_verde_off(i),index_verde_off(i)] = min(abs(time-timepoint_verde_off(i)));
end
clear time_rojo time_rojo_off time_verde time_verde_off
clear timepoint_rojo timepoint_rojo_off timepoint_verde timepoint_verde_off

idx_verde=[index_verde,index_verde_off]; %joint Start(ON) and End(OFF) of Congurnet Events
idx_rojo=[index_rojo,index_rojo_off]; %joint Start(ON) and End(OFF) of Incongurnet Events
clear index_verde index_verde_off index_rojo index_rojo_off interval_rojo interval_verde
clear LowCut HighCut FilterIndex  FilterIndex1 FilterIndex2 fmd A B

%% Not take into account events smaller than criteria
%Median
median_V=median((idx_verde(:,2)-idx_verde(:,1)),1);
median_R=median((idx_rojo(:,2)-idx_rojo(:,1)),1);

if median_V<median_R
    median=median_V;
else
    median=median_R;
end

if point==1
    criterio=median;
else
    criterio=criterio;
end

clear median_V median_R

%Congruents
idx_verde_tmp=[ ];
for c=1:length(idx_verde)
    if(idx_verde(c,2)-idx_verde(c,1))>criterio
    idx_verde_tmp=[idx_verde_tmp; idx_verde(c,:)]; 
    end
end
idx_verde=idx_verde_tmp;
clear c idx_verde_tmp

%Incongruent
idx_rojo_tmp=[];
for c=1:length(idx_rojo)
    if (idx_rojo(c,2)-idx_rojo(c,1))>criterio
        idx_rojo_tmp=[idx_rojo_tmp;idx_rojo(c,:)];
    end
end
idx_rojo=idx_rojo_tmp;
clear c idx_rojo_tmp

%% Take LFP segments and calculates Hilbert
HP_R = [];
HP_V = [];
PFC_R = [];
PFC_V = [];
%Congruent Events
for i = 1:length(idx_verde)
    %Indexes
    ii = idx_verde(i,1);
    iii = idx_verde(i,2);
    %Restriction of Signal
    HP_tmp = HP(ii : iii,1);
    PFC_tmp = PFC(ii:iii,1);
    %Hippocampus
    %hamming y and conncatenation
    H = hamming(length(HP_tmp),'symmetric'); %Hamming window generation
    HP_tmp=HP_tmp.*H; %apply of hamming
    HP_V=[HP_V;HP_tmp]; %concatenation
    %Prefrontal
    H = hamming(length(PFC_tmp),'symmetric');
    PFC_tmp=PFC_tmp.*H;
    PFC_V=[PFC_V;PFC_tmp];
     clear min max ii iii H HP_tmp PFC_tmp
end

%Incongruent Events
for i = 1:length(idx_rojo)
    %Index
    ii = idx_rojo(i,1);
    iii = idx_rojo(i,2);

    HP_tmp = HP(ii : iii,1);
    PFC_tmp = PFC(ii:iii,1);
    
    H = hamming(length(HP_tmp),'symmetric');
    HP_tmp=HP_tmp.*H;
    HP_R=[HP_R;HP_tmp];
    
    H = hamming(length(PFC_tmp),'symmetric');
    PFC_tmp=PFC_tmp.*H;
    PFC_R=[PFC_R;PFC_tmp];
    clear ii iii H HP_tmp PFC_tmp
end
clear i

%% Equal length of signals
if length(HP_R)>length(HP_V)
    T = length(HP_V);
else
    T = length(HP_R);
end

HP_R = HP_R(1:T,1);
HP_V = HP_V(1:T,1);
PFC_R = PFC_R(1:T,1);
PFC_V = PFC_V(1:T,1);

%% Normalization to Max or Min
%Hippocampus
%Congruent
max = max(HP_V); %Max
min = abs(min(HP_V)); %absolute minnimal value
if max > min 
   HP_V = HP_V/max; 
else
    HP_V = HP_V/min;
end
clear max min
%Incongruent
max = max(HP_R);
min = abs(min(HP_R));
if max > min
   HP_R = HP_R/max; 
else
    HP_R = HP_R/min;
end
clear max min

%Prefrontal
%Congruent
max = max(PFC_V); 
min = abs(min(PFC_V)); 
if max > min
   PFC_V = PFC_V/max; 
else
    PFC_V = PFC_V/min; 
end
clear max min
%Incongruent
max = max(PFC_R);
min = abs(min(PFC_R));
if max > min 
   PFC_R = PFC_R/max; 
else
    PFC_R = PFC_R/min;
end
clear max min

%% Hilbert transform
hilbert_HP_R = hilbert(HP_R);
hilbert_HP_V = hilbert(HP_V);
hilbert_PFC_R = hilbert(PFC_R);
hilbert_PFC_V = hilbert(PFC_V);

%% Amplitude
V1 = abs(hilbert_HP_V); %Hippocampus Congruent
R1 = abs(hilbert_HP_R); %Hippocampus Incongruent
V2 = abs(hilbert_PFC_V); %Prefrontal Congruent
R2 = abs(hilbert_PFC_R); %Prefrontal Incongruent

%substraccion of DC component
V1 = V1-mean(V1,1);
R1 = R1-mean(R1,1);
V2 = V2-mean(V2,1);
R2 = R2-mean(R2,1);

%% Ccross-correlacion
[V,lagV] = xcorr(V1,V2,lag,'coeff'); %congruent
[R,lagR] = xcorr(R1,R2,lag,'coeff'); %incongruent

lagV=(lagV./fm)*1000; %converts lags to miliseconds
lagR=(lagR./fm)*1000; %converts lags to miliseconds

%% Save
save(([PathName,'cross_amp6.mat']),'V','R','lagR','lagV')