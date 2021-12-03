% This script thake the ON or the OFF of the exploration event, cut
% surrounding LFP and calculate the individual and average Spectrogram
%------
% Inputs : LFP and Events(timestamps in sec)
% Outputs: Individual and Average Power, Cross-Power and Coherograms
%------

clear
clc
close all

%% Parameters
currentdir=('D:myPath');
listing=dir(currentdir);

%temporal limits to look for the maximal value in the spectrogram
% ON
desde = 2500;
hasta = 3250;
% %desde OFF
% desde = 1875;
% hasta = 2500;

fs = 1250; %sampling frequecy
dt = 1/fs;

%Channels to upload
chHP = [24;20;20;20;20;20;20;20;20;20;20;2;2;2;18;18;20;20;20;20;18];
chPFC = [1;1;5;1;13;5;5;5;1;5;1;17;17;25;13;5;13;5;5;11;12];

% variables to fill with the data
vec1 = [];
vec2 = [];
vec3 = [];
vec4 = [];
vec5 = [];
vec6 = [];
vec7 = [];
vec8 = [];
vec9 = [];
vec10 = [];

vecx = [];

% maximal values
promC = [];
promI = [];

intervalosC = [];
intervalosI = [];

for y = 1:length(listing)-2
    %% Load and detrend of LFP
    %HIPOCAMPUS
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
    size_T=length(lfpHP); %samples
    time=((TimeStamps(1):((1/1250)*1000000):(((size_T/1250)*1000000)+TimeStamps(1)-((1/1250)*1000000)))); %time vector
    time=time';
    
    %% Find the Start(ON) and End(OFF) of the exploration events
    timepoint_rojo = interval_rojo(:,1);
    timepoint_rojo_off = interval_rojo(:,2);
    timepoint_verde = interval_verde(:,1);
    timepoint_verde_off = interval_verde(:,2);

    %Incongruent
    index_rojo = zeros(length(timepoint_rojo),1); 
    time_rojo = zeros(length(timepoint_rojo),1);
    index_rojo_off = zeros(length(timepoint_rojo_off),1);
    time_rojo_off = zeros(length(timepoint_rojo_off),1);
    
    
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
    

    %parameters of filtering y SF
    LowCut=2/(fs/2);  %high-pass
    HighCut=624/(fs/2);  %low-pass
    %butterworth filter
    [B,A]=butter(3,[LowCut, HighCut]);
    
    %HP
    lfpHP= filtfilt(B,A,lfpHP(:,1));
    %PFC
    lfpPFC= filtfilt(B,A,lfpPFC(:,1));
    clear A B

    
    %% Selection of events longer than 0.5 sec
    m = [];
    for h = 1:length(index_rojo)
        if index_rojo_off(h)-index_rojo(h) > 625
            m = [m,index_rojo_off(h)];
%             m = [m,index_rojo(h)];

        end
    end
    index_rojo = m';
    clear m

    m = [];
    for h = 1:length(index_verde)
        if index_verde_off(h)-index_verde(h) > 625
            m = [m,index_verde_off(h)];
%             m = [m,index_verde(h)];
        end
    end
    index_verde = m';
    clear m    
    
    r = length(index_rojo);
    v = length(index_verde);
    red = index_rojo;
    green = index_verde;
    
    %Suffle the events 20 times
for yy = 1:20
         red = shuffle(red);
         green = shuffle(green);
end
    
    if r>v
        index_rojo = red(1:v,1);
    else
        index_verde = green(1:r,1);
    end
    clear green red r v
    
    %% Perievents Histograms generation
    %Congruent
    tmp = [];
    tmp1 = [];
    tmp2 = [];
    tmp3 = [];
    tmp4 = [];
    W = 625;
    window = hanning(W); %window of elimination in sec
    noverlap = round(W*0.90);
    p = 2;
    m = 50;
    freq = [0:0.5:m];

    
    for i = 1 : length(index_verde)
        %Hippocampus
        lfp = lfpHP(index_verde(i)-fs*p:index_verde(i)+fs*p);
        [s,f,t] = spectrogram(lfp,window,noverlap,freq,fs,'power');
        if i == 1
            tmp = zeros(size(s));
            tmp = tmp + log(abs(s));
        else
            tmp = tmp + log(abs(s));
        end
        clear lfp s
        
        
        %Prefrontal
        lfp = lfpPFC(index_verde(i)-fs*p:index_verde(i)+fs*p);
        [s,f,t] = spectrogram(lfp,window,noverlap,freq,fs,'power');
        if i == 1
            tmp1 = zeros(size(s));
            tmp1 = tmp1 + log(abs(s));
        else
            tmp1 = tmp1 + log(abs(s));
        end
        clear lfp s
        
        %cross
        lfp1 = lfpHP(index_verde(i)-fs*p:index_verde(i)+fs*p);
        lfp2 = lfpPFC(index_verde(i)-fs*p:index_verde(i)+fs*p);
        
        [s,f,t] = xspectrogram(lfp1,lfp2,window,noverlap,freq,fs,'power');
        if i == 1
            tmp2 = zeros(size(s));
            tmp2 = tmp2 + log(abs(s));
        else
            tmp2 = tmp2 + log(abs(s));
        end
        clear lfp1 lfp2 s
        
        
        %coherence
        lfp1 = lfpHP(index_verde(i)-fs*p:index_verde(i)+fs*p);
        lfp2 = lfpPFC(index_verde(i)-fs*p:index_verde(i)+fs*p);
        [c,wcs,fc] = wcoherence(lfp1,lfp2,fs,'FrequencyLimits',[0 50],'NumScalesToSmooth',18);
        %         c = c/cb;
        [Pxy,F] = cpsd(lfp1,lfp2,window,noverlap,freq,fs);
        if i == 1
            tmp3 = c;
            tmp4 = Pxy;
        else
            tmp3 = tmp3 + c;
            tmp4 = [tmp4,Pxy];
        end
        clear lfp1 lfp2 Pxy c
    end
    PxxC = tmp./length(index_verde);%POWER HIPOCAMPÖ CONGRUENTE
    PyyC = tmp1./length(index_verde);%POWER HIPOCAMPÖ CONGRUENTE
    PxyC = tmp2./length(index_verde);%CROSS HIP-PFC CONGRUENTE
    CxyC = tmp3./length(index_verde);%COHERENCIA CONGRUENTE
    
    PhaseC = mean((180/pi)*angle(tmp4),2);
    
    
    clear tmp tmp1 tmp2 tmp3 tmp4 i
    
    %Incongruent
    tmp = [];
    tmp1 = [];
    tmp2 = [];
    tmp3 = [];
    tmp4 = [];
    
    tmpx = [];
    for i = 1 : length(index_rojo)
        
        %Hippocampus
        lfp = lfpHP(index_rojo(i)-fs*p:index_rojo(i)+fs*p);
        [s,f,t] = spectrogram(lfp,window,noverlap,freq,fs,'power');
        if i == 1
            tmp = zeros(size(s));
            tmp = tmp + log(abs(s));
            
            tmpx = real(log(s));
        else
            tmp = tmp + log(abs(s));
            tmpx = tmpx + real(log(s));
        end
        clear lfp s
        
        
        %Prefrontal
        lfp = lfpPFC(index_rojo(i)-fs*p:index_rojo(i)+fs*p);
        [s,f,t] = spectrogram(lfp,window,noverlap,freq,fs,'power');
        if i == 1
            tmp1 = zeros(size(s));
            tmp1 = tmp1 + log(abs(s));
        else
            tmp1 = tmp1 + log(abs(s));
        end
        clear lfp s
        
        %cross
        lfp1 = lfpHP(index_rojo(i)-fs*p:index_rojo(i)+fs*p);
        lfp2 = lfpPFC(index_rojo(i)-fs*p:index_rojo(i)+fs*p);
        
        [s,f,t] = xspectrogram(lfp1,lfp2,window,noverlap,freq,fs,'power');
        if i == 1
            tmp2 = zeros(size(s));
            tmp2 = tmp2 + log(abs(s));
        else
            tmp2 = tmp2 + log(abs(s));
        end
        clear lfp1 lfp2 s
        
        %coherence
        lfp1 = lfpHP(index_rojo(i)-fs*p:index_rojo(i)+fs*p);
        lfp2 = lfpPFC(index_rojo(i)-fs*p:index_rojo(i)+fs*p);
        [c,wcs,fc] = wcoherence(lfp1,lfp2,fs,'FrequencyLimits',[0 50],'NumScalesToSmooth',18);
        [Pxy,F] = cpsd(lfp1,lfp2,window,noverlap,freq,fs);
        %         c = c/cb;
        if i == 1
            tmp3 = c;
            
            tmp4 = Pxy;
            
        else
            tmp3 = tmp3 + c;
            tmp4 = [tmp4,Pxy];
        end
        clear lfp1 lfp2 c Pxy
    end
    PxxI = tmp./length(index_rojo);%POWER HIPOCAMPÖ CONGRUENTE
    
    PxxIx = tmpx./length(index_rojo);
    
    PyyI = tmp1./length(index_rojo);%POWER HIPOCAMPÖ CONGRUENTE
    PxyI = tmp2./length(index_rojo);%CROSS HIP-PFC CONGRUENTE
    CxyI = tmp3./length(index_rojo);%COHERENCIA CONGRUENTE
    
    PhaseI = mean((180/pi)*angle(tmp4),2);

    
    clear tmp tmp1 tmp2 tmp3 tmp4 i cb wcsb fcb
    
    promC = [promC;max(max(CxyC(31:35,desde:hasta)))];
    promI = [promI;max(max(CxyI(31:35,desde:hasta)))];

    if y == 1
        vec1 = PxxC;
        vec2 = PyyC;
        vec3 = PxyC;
        vec4 = CxyC;
        vec5 = PhaseC;
        
        vec6 = PxxI;
        vec7 = PyyI;
        vec8 = PxyI;
        vec9 = CxyI;
        vec10 = PhaseI;
        
        vecx = tmpx;
        
    else
        vec1 = vec1 + PxxC;
        vec2 = vec2 + PyyC;
        vec3 = vec3 + PxyC;
        vec4 = vec4 + CxyC;
        vec5 = [vec5,PhaseC];
        
        vec6 = vec6 + PxxI;
        vec7 = vec7 + PyyI;
        vec8 = vec8 + PxyI;
        vec9 = vec9 + CxyI;
        vec10 = [vec10,PhaseI];
        
        vecx = vecx + PxxIx;
        
    end
    
    clear FileName FileName2 FilterIndex2 fm_original freq HighCut ini inicio lfpHP lfpPFC LowCut
    clear max min noverlap path PathName2 size_T time timepoint_rojo timepoint_rojo_off timepoint_verde timepoint_verde_off
    clear TimeStamps W window Filename intervalR intervalV r v
    
    % Save the individual matrix
    save ([file.folder,'\espectrograma_individual_OFF.mat'],'PxxC','PyyC','PxyC','CxyC','PhaseC','PxxI','PyyI','PxyI','CxyI','PhaseI','PxxIx','index_rojo','index_rojo_off','index_verde','index_verde_off');

    clear PxxC PyyC PxyC CxyC PhaseC PxxI PyyI PxyI CxyI PhaseI PxxIx
    clear index_rojo index_rojo_off index_verde index_verde_off
    
end
clear chHP chPFC fs y currentdir


% Average of all the outputs
PxxC = vec1/(length(listing)-2); %Power HPC Congruent
PyyC = vec2/(length(listing)-2); %Power PFC Congruent
PxyC = vec3/(length(listing)-2); %Cross-Power HPC-PFC Congruent
CxyC = vec4/(length(listing)-2); % Coherence HPC-PFC Congruent
PhaseC = vec5;

PxxIx = vecx/(length(listing)-2);

PxxI = vec6/(length(listing)-2); %Power HPC Incongruent
PyyI = vec7/(length(listing)-2); %Power PFC Incongruent
PxyI = vec8/(length(listing)-2); %Cross-Power HPC-PFC Incongruent
CxyI = vec9/(length(listing)-2); % Coherence HPC-PFC Incongruent
PhaseI = vec10;
clear vec1 vec2 vec3 vec4 vec5 vec6 vec7 vec8 vec9 vec10

%% Save
currentdir=('D:\myPath');
save ([currentdir,'\espectrograma_OFF.mat'])
