% This script thake the ON or the OFF of the exploration event, it split
% the events dthree segments and cut surrounding LFP to calculate the 
% individual and average Spectrogram
%------
% Inputs : LFP and Events(timestamps in sec)
% Outputs: Individual  and Average Coherograms across the session
%------
clear
clc
close all

%% Parameters
currentdir=('D:\myPath');
listing=dir(currentdir);

%temporal limits to look for the maximal value in the spectrogram
% ON
desde = 2625;
hasta = 3250;
% OFF
% desde = 1875;
% hasta = 2500;

% Position of the frequency range to find the maximal value
abajo = 32;
arriba = 34;

% numero = 5;
scale = 18; %smooth of the coherogram
fs = 1250; %sampling frequecy
dt = 1/fs;

%Channels to upload
chHP = [24;20;20;20;20;20;20;20;20;20;20;2;2;2;18;18;20;20;20;20;18];
chPFC = [1;1;5;1;13;5;5;5;1;5;1;17;17;25;13;5;13;5;5;11;12];

vec1 = [];
vec2 = [];
vec3 = [];
vec4 = [];
vec5 = [];
vec6 = [];


promC = [];
promI = [];

intervalosC = [];
intervalosI = [];

durations = cell(length(listing)-2,2);

for y = 1:length(listing)-2
    %% Load and detrend of LFP
    %HIPOCAMPO
    
    file = dir(fullfile([currentdir,'\',listing(y+2).name,],'*.lfp'));
    
    %load .lfp
    lfpHP = LoadBinary([file.folder,'\',file.name],'channels',chHP(y),'frequency',fs,'nChannels',32);
    lfpHP = detrend(lfpHP,0); %saca el drift
    
    lfpPFC = LoadBinary([file.folder,'\',file.name],'channels',chPFC(y),'frequency',fs,'nChannels',32);
    lfpPFC = detrend(lfpPFC,0); %saca el drift
    
    %% Load the Events.mat file
    %     fm_original=32556; %frecuencia de muestreo original
    events =dir(fullfile([currentdir,'\',listing(y+2).name],'*.mat'));
    load([currentdir,'\',listing(y+2).name,'\Events_1.mat'],'TimeStamps','interval_rojo','interval_verde','ID_separada');
    size_T=length(lfpHP);
    time=((TimeStamps(1):((1/1250)*1000000):(((size_T/1250)*1000000)+TimeStamps(1)-((1/1250)*1000000)))); %time vector generation
    time=time'; 
    
    %% Find the Start(ON) and End(OFF) of the exploration events
    timepoint_rojo=interval_rojo(:,1);
    timepoint_rojo_off=interval_rojo(:,2);
    timepoint_verde=interval_verde(:,1);
    timepoint_verde_off=interval_verde(:,2);
    
    %Incongruent
    index_rojo=zeros(length(timepoint_rojo),1); %Esta variable almacenara el indice de cada timestamp dentro del vector time.
    time_rojo=zeros(length(timepoint_rojo),1); %esta variable no estoy seguro para que sirve. No lo uso mas adelante
    index_rojo_off=zeros(length(timepoint_rojo_off),1);
    time_rojo_off=zeros(length(timepoint_rojo_off),1);
    
    
    for i=1:length(timepoint_rojo)
        [time_rojo(i),index_rojo(i)] = min(abs(time-timepoint_rojo(i)));
        [time_rojo_off(i),index_rojo_off(i)] = min(abs(time-timepoint_rojo_off(i)));
    end
    intervalosI = [intervalosI;(index_rojo_off-index_rojo)/1250];
    
    %Congruent
    index_verde=zeros(length(timepoint_verde),1);
    time_verde=zeros(length(timepoint_verde),1);
    index_verde_off=zeros(length(timepoint_verde_off),1);
    time_verde_off=zeros(length(timepoint_verde_off),1);
    
    for i=1:length(timepoint_verde)
        [time_verde(i),index_verde(i)] = min(abs(time-timepoint_verde(i)));
        [time_verde_off(i),index_verde_off(i)] = min(abs(time-timepoint_verde_off(i)));
        
    end
    clear i
    intervalosC = [intervalosC;(index_verde_off-index_verde/1250)];
    
    verdes = [index_verde,index_verde_off];
    rojos = [index_rojo,index_rojo_off];
    clear time_verde  time_verde_off %index_verde_off index_verde
    clear time_rojo time_rojo_off %index_rojo_off index_rojo
    clear interval_rojo interval_verde ID_separada
    
    
    index_verde = index_verde(3:end,1);
    index_verde_off = index_verde_off(3:end,1);
    
    index_rojo = index_rojo(3:end,1);
    index_rojo_off = index_rojo_off(3:end,1);
    
    % Selection of events longer than 0.5 sec
    m = [];
    n = [];
    for h = 1:length(index_rojo)
        if index_rojo_off(h)-index_rojo(h) > 625
            m = [m,index_rojo(h)];
            n = [n;(index_rojo_off(h)-index_rojo(h))/fs];
        end
    end
    index_rojo = m';
    durations{y,1} = n;
    clear m n    
    
    n = [];
    m = [];
    for h = 1:length(index_verde)
        if index_verde_off(h)-index_verde(h) > 625
            m = [m,index_verde(h)];
            n = [n;(index_verde_off(h)-index_verde(h))/fs];
            
        end
    end
    index_verde = m';
    durations{y,2} = n;
    clear m n
    
    % definition of the number of events per segment
    if length(index_rojo) > length(index_verde)
        xx = length(index_verde);
        index_rojo = index_rojo(1:xx);
    else
        xx = length(index_rojo);
        index_verde = index_verde(1:xx);
    end
    
    if mod(xx/3,1) == 0
        numero = xx/3;
    else
        numero = round(xx/3)-1;
    end
    clear xx
    
    % parameters of filtering y SF
    LowCut=2/(fs/2);  %high-pass
    HighCut=624/(fs/2);  %low-pass
    % butterworth filter
    [B,A]=butter(3,[LowCut, HighCut]);
    
    %HP
    lfpHP= filtfilt(B,A,lfpHP(:,1));
    %PFC
    lfpPFC= filtfilt(B,A,lfpPFC(:,1));
    clear A B
    
    %% Perievents Histograms generation
    %Congruent
    tmp1 = [];
    tmp2 = [];
    tmp3 = [];
    tmp4 = [];
    
    W = 625; %window size
    window = hanning(W); %window of elimination in sec
    noverlap = round(W*0.90);
    p = 2;
    m = 50;
    freq = [0:0.1:m];
    
    
    count1 = 0;
    count2 = 0;
    count3 = 0;
    count4 = 0;
    
    
    for i = 1 : length(index_verde)
        if i <= numero
            %coherence
            lfp1 = lfpHP(index_verde(i)-fs*p:index_verde(i)+fs*p);
            lfp2 = lfpPFC(index_verde(i)-fs*p:index_verde(i)+fs*p);
            [c,wcs,fc] = wcoherence(lfp1,lfp2,fs,'FrequencyLimits',[0 50],'NumScalesToSmooth',scale);
            if isempty(tmp1)
                tmp1 = c;
                count1 = count1 + 1;
            else
                tmp1 = tmp1 + c;
                count1 = count1 + 1;
            end
            clear lfp1 lfp2 c
            
        elseif i <= numero*2  && i > numero
            %coherence
            lfp1 = lfpHP(index_verde(i)-fs*p:index_verde(i)+fs*p);
            lfp2 = lfpPFC(index_verde(i)-fs*p:index_verde(i)+fs*p);
            [c,wcs,fc] = wcoherence(lfp1,lfp2,fs,'FrequencyLimits',[0 50],'NumScalesToSmooth',scale);
            if isempty(tmp2)
                tmp2 = c;
                count2 = count2 + 1;
            else
                tmp2 = tmp2 + c;
                count2 = count2 + 1;
            end
            clear lfp1 lfp2 c
            
        elseif i <= numero*3 && i > numero*2
            %coherence
            lfp1 = lfpHP(index_verde(i)-fs*p:index_verde(i)+fs*p);
            lfp2 = lfpPFC(index_verde(i)-fs*p:index_verde(i)+fs*p);
            [c,wcs,fc] = wcoherence(lfp1,lfp2,fs,'FrequencyLimits',[0 50],'NumScalesToSmooth',scale);
            if isempty(tmp3)
                tmp3 = c;
                count3 = count3 + 1;
            else
                tmp3 = tmp3 + c;
                count3 = count3 + 1;
            end
            clear lfp1 lfp2 c
        end
    end
    CxyC1 = tmp1./count1;%COHERENCE CONGRUENT PART 1
    CxyC2 = tmp2./count2;%COHERENCE CONGRUENT PART 2
    CxyC3 = tmp3./count3;%COHERENCE CONGRUENT PART 3
    
    
    clear tmp1 tmp2 tmp3 tmp4 count1 count2 count3 count4 i
    
    %iNCONGRUENT
    tmp1 = [];
    tmp2 = [];
    tmp3 = [];
    tmp4 = [];
    count1 = 0;
    count2 = 0;
    count3 = 0;
    count4 = 0;

    
    for i = 1 : length(index_rojo)
        if i <= numero
            %coherence
            lfp1 = lfpHP(index_rojo(i)-fs*p:index_rojo(i)+fs*p);
            lfp2 = lfpPFC(index_rojo(i)-fs*p:index_rojo(i)+fs*p);
            [c,wcs,fc] = wcoherence(lfp1,lfp2,fs,'FrequencyLimits',[0 50],'NumScalesToSmooth',scale);
            if isempty(tmp1)
                tmp1 = c;
                count1 = count1 + 1;
            else
                tmp1 = tmp1 + c;
                count1 = count1 + 1;
            end
            
        elseif i <= numero*2  && i > numero
            %coherence
            lfp1 = lfpHP(index_rojo(i)-fs*p:index_rojo(i)+fs*p);
            lfp2 = lfpPFC(index_rojo(i)-fs*p:index_rojo(i)+fs*p);
            [c,wcs,fc] = wcoherence(lfp1,lfp2,fs,'FrequencyLimits',[0 50],'NumScalesToSmooth',scale);
            if isempty(tmp2)
                tmp2 = c;
                count2 = count2 + 1;
            else
                tmp2 = tmp2 + c;
                count2 = count2 + 1;
            end
        elseif i <= numero*3 && i > numero*2
            %coherence
            lfp1 = lfpHP(index_rojo(i)-fs*p:index_rojo(i)+fs*p);
            lfp2 = lfpPFC(index_rojo(i)-fs*p:index_rojo(i)+fs*p);
            [c,wcs,fc] = wcoherence(lfp1,lfp2,fs,'FrequencyLimits',[0 50],'NumScalesToSmooth',scale);
            if isempty(tmp3)
                tmp3 = c;
                count3 = count3 + 1;
            else
                tmp3 = tmp3 + c;
                count3 = count3 + 1;
            end

        end
        clear lfp1 lfp2 c wcs
    end
    CxyI1 = tmp1./count1;%COHERENCE INCONGRUENT PART 1
    CxyI2 = tmp2./count2;%COHERENCE INCONGRUENT PART 2
    CxyI3 = tmp3./count3;%COHERENCE INCONGRUENT PART 3
    
    
    clear tmp1 tmp2 tmp3 tmp4 count1 count2 count3 count4 i
    

    promC = [promC;max(CxyC1(abajo:arriba,desde:hasta),[],'all'),max(CxyC2(abajo:arriba,desde:hasta),[],'all'),max(CxyC3(abajo:arriba,desde:hasta),[],'all')],%,max(max(CxyC4(abajo:arriba,desde:hasta)))];
    promI = [promI;max(CxyI1(abajo:arriba,desde:hasta),[],'all'),max(CxyI2(abajo:arriba,desde:hasta),[],'all'),max(CxyI3(abajo:arriba,desde:hasta),[],'all')],%,max(max(CxyI4(abajo:arriba,desde:hasta)))];
    
    
    if y == 1 %In this variables matrices are summed
        vec1 = CxyC1;
        vec2 = CxyC2;
        vec3 = CxyC3;
        vec5 = CxyI1;
        vec6 = CxyI2;
        vec7 = CxyI3;
    else
        vec1 = vec1 + CxyC1;
        vec2 = vec2 + CxyC2;
        vec3 = vec3 + CxyC3;
        vec5 = vec5 + CxyI1;
        vec6 = vec6 + CxyI2;
        vec7 = vec7 + CxyI3;
    end
    
    clear FileName FileName2 FilterIndex2 fm_original freq HighCut ini inicio lfpHP lfpPFC LowCut
    clear max min noverlap path PathName2 size_T time timepoint_rojo timepoint_rojo_off timepoint_verde timepoint_verde_off
    clear TimeStamps W window Filename intervalR intervalV r v
    
    % Save of individual matrices
    save (['D:\para cami\',file.name(1:end-4),'.mat'],'CxyC1','CxyC2','CxyC3','CxyI1','CxyI2','CxyI3','fc');
    
    clear CxyC1 CxyC2 CxyC3 CxyI1 CxyI2 CxyI3 t1 t2 t3 t4 wcs
    clear index_rojo index_rojo_off index_verde index_verde_off xx numero
    
end
clear chHP chPFC fs y currentdir
[H,p] = ttest(promI(:,1),promC(:,1))

% Average of Martrices
CxyC1 = vec1./(length(listing)-2);
CxyC2 = vec2./(length(listing)-2);
CxyC3 = vec3./(length(listing)-2);
CxyI1 = vec5./(length(listing)-2);
CxyI2 = vec6./(length(listing)-2);
CxyI3 = vec7./(length(listing)-2);

clear vec1 vec2 vec3 vec4 vec5 vec6 vec7 vec8

%% Save Final Matrix
save('C:\myPath')
