%-----------------------
% This script analyse the phase-locking of PFC Units to HPC theta rhythm
% Upload .clu and .res files to use spkies' timestamps
% Upload .mat file to use events' timestamps
% Upload Hippocampal LFP, apply hilbert transform in the filtered signal (6-10 Hz)
% Rayleigh to calculate probability of circular distribution with SU phases.
% Output:
%       1) p = probability of uniform distribution
%       2) theta = Vectors' angle
%       3) rbar = Module of the Vector
%       4) delta = dispersion
%-----------------------

clear
clc
close all
%% Parameters

tipo_trial = 2;         % Number of events
freqmLFP = 1250;        % LFP sampling frequency
criterio = 25;          % Minimal number of Spiks to perform Rayleight

% Frequency band to filter the signal
bajo = 6;
alto = 10;

nro_tetrodos = 4; %	Number of tetrodes
ch = [18]; % Hippocampal channel

%% Upload of .mat and .lfp files

folder=('myPath'); %change deppending on the session
[nombre,path] = uigetfile('/*.mat','Seleccionar archivo con los eventos a analizar', folder);

if isequal(nombre,0) || isequal(path,0)
    disp('Apretaste cancelar')
    return
else
    disp(['Sesion seleccionada:   ', nombre(end-9:end-4)])
end

% Load of tje events
nomarch = [path,nombre];
load(nomarch,'EventStrings','EventStrings','TimeStamps','timestamps','ttl2bin','interval_rojo','interval_verde'); %carga las variables al workspace

nombre_LFP = [path(end-16:end-1),'.lfp']; %define el nombre del archivo LFP a cargar
path_LFP = [path];

nomarch = [path_LFP, nombre_LFP];
fid = fopen(nomarch); %load of LFP
data = fread(fid, [1,inf], 'int16');
fclose(fid);

% Load of LFP
canales = reshape(data, [32,(length(data)/32)]); %
ch_LFP = [canales(ch(1),:)'];
clear canales

%% Single-Units (SU)
% One cell per SU
clusters = cell(4,1); %per tetrode
nro_neuronas = zeros(4,1);

for nr = 1 : 4
    file_t = fopen([path,path(end-16:end-1),'.res.',num2str(nr)]); %Open .res (32550 HZ)
    file_c = fopen([path,path(end-16:end-1),'.clu.',num2str(nr)]); %Open .clu

    timestamps_n = fscanf(file_t,'%i'); %Timestamps of each SU
    clust = fscanf(file_c,'%i'); %Cluster of each Spike

    fm = (1/32556)*1000000; %Trqnsformqtion of Sampling Frequwncy to microsec
    timestamps_n = timestamps_n*fm; %from time to position

    clusters{nr} = [clust(2:end),timestamps_n+TimeStamps(1)]; %indexing of clusters and Spks
    nro_neuronas(nr) = max(clust(2:end))-1;
end

for u = 1 : length(nro_neuronas)
    if nro_neuronas(u) < 0
        nro_neuronas(u) = 0;
    end
end


%% Split SU according to their cluster id
nr_t1 = cell(nro_neuronas(1),1);
nr_t2 = cell(nro_neuronas(2),1);
nr_t3 = cell(nro_neuronas(3),1);
nr_t4 = cell(nro_neuronas(4),1);

%Tetrode 1
for ind = 2 : nro_neuronas(1) + 1 
    aux = find(clusters{1}(:,1) == ind);
    neurona = [clusters{1}(aux,1),clusters{1}(aux,2)];
    nr_t1{ind-1} = neurona;
end
%Tetrode 2
for ind = 2 : nro_neuronas(2) + 1
    aux = find(clusters{2}(:,1) == ind);
    neurona = [clusters{2}(aux,1),clusters{2}(aux,2)];
    nr_t2{ind-1} = neurona;
end
%Tetrode 3
for ind = 2 : nro_neuronas(3) + 1
    aux = find(clusters{3}(:,1) == ind);
    neurona = [clusters{3}(aux,1),clusters{3}(aux,2)];
    nr_t3{ind-1} = neurona;
end
%Tetrode 4
for ind = 2 : nro_neuronas(4) + 1
    aux = find(clusters{4}(:,1) == ind);
    neurona = [clusters{4}(aux,1),clusters{4}(aux,2)];
    nr_t4{ind-1} = neurona;
end
%Store
neuronas = cell(nro_tetrodos,1);
neuronas{1} = nr_t1;
neuronas{2} = nr_t2;
neuronas{3} = nr_t3;
neuronas{4} = nr_t4;

%% Label for the data
label1 = 'Congruent';
label2 = 'Icongruent';
label_t = cell(tipo_trial,1);
label_t{1} = label1;
label_t{2} = label2;


label1 = 'Tetrode 1';
label2 = 'Tetrode 2';
label3 = 'Tetrode 3';
label4 = 'Tetrode 4';
label_c = cell(nro_tetrodos,1);
label_c{1} = label1;
label_c{2} = label2;
label_c{3} = label3;
label_c{4} = label4;

clear label1 label2 label3 label4

%% Find Start (ON) and End (OFF) of events
timepoint_rojo=interval_rojo(:,1);
timepoint_rojo_off=interval_rojo(:,2);
timepoint_verde=interval_verde(:,1);
timepoint_verde_off=interval_verde(:,2);

%time vector for the LFP
size_T=length(ch_LFP); 
tiempo=((TimeStamps(1):((1/1250)*1000000):(((size_T/1250)*1000000)+TimeStamps(1)-((1/1250)*1000000)))); %genera vector tiempo utilizando los TimeStamps
tiempo=tiempo';


% Find closest position of the events the LFP time vector
index_rojo=zeros(length(timepoint_rojo),1); 
time_rojo=zeros(length(timepoint_rojo),1); 
% incongruent ON
for i=1:length(timepoint_rojo)
    [time_rojo(i),index_rojo(i)] = min(abs(tiempo-timepoint_rojo(i)));
end
% incongruent OFF
index_rojo_off=zeros(length(timepoint_rojo_off),1);
time_rojo_off=zeros(length(timepoint_rojo_off),1);
for i=1:length(timepoint_rojo_off)
    [time_rojo_off(i),index_rojo_off(i)] = min(abs(tiempo-timepoint_rojo_off(i)));
end
% congruent ON
index_verde=zeros(length(timepoint_verde),1);
time_verde=zeros(length(timepoint_verde),1);
for i=1:length(timepoint_verde)
    [time_verde(i),index_verde(i)] = min(abs(tiempo-timepoint_verde(i)));
end
% congruent OFF
index_verde_off=zeros(length(timepoint_verde_off),1);
time_verde_off=zeros(length(timepoint_verde_off),1);
for i=1:length(timepoint_verde_off)
    [time_verde_off(i),index_verde_off(i)] = min(abs(tiempo-timepoint_verde_off(i)));
end

% Only segments larger that 1 seconds is taken into account

idx_verde=[index_verde,index_verde_off];
idx_rojo=[index_rojo,index_rojo_off];
clear index_verde index_verde_off index_rojo index_rojo_off

%Congruent
idx_verde_tmp=[];
for c=1:length(idx_verde)
    if (idx_verde(c,2)-idx_verde(c,1))>1250
        idx_verde_tmp=[idx_verde_tmp;idx_verde(c,:)];
    end
end
idx_verde_tmp=idx_verde;
clear c idx_verde_tmp

%Incongruent
idx_rojo_tmp=[];
for c=1:length(idx_rojo)
    if (idx_rojo(c,2)-idx_rojo(c,1))>1250
        idx_rojo_tmp=[idx_rojo_tmp;idx_rojo(c,:)];
    end
end
idx_rojo_tmp=idx_rojo;
clear c idx_rojo_tmp

%% Generation of LFP for each condition
signal = cell(nro_tetrodos,tipo_trial); 

% store of time
 t_signalV = cell(1,length(idx_verde)); %congruent
 t_signalR = cell(1,length(idx_rojo)); %incongruent
 
% LOOP: It iterates channel by channel.
% It generates a matrix where the LFP is store
% To do so, it takes each segment and relativize the time to the start

for idx = 1 : length(ch)
    
    %congruent
    senial_V=cell(1,length(idx_verde));
    for i = 1 : length(idx_verde)
        senial_V{1,i}=ch_LFP(idx_verde(i,1):idx_verde(i,2),idx); %LFP
        t_signalV{1,i} = tiempo(idx_verde(i,1):idx_verde(i,2),1); %time
    end
    signal{idx,1} = senial_V; %store of LFP
    
    
    %incongruent
    senial_R = cell(1, length(idx_rojo));
    for i = 1 : length(idx_rojo)
          senial_R{1,i}=ch_LFP(idx_rojo(i,1):idx_rojo(i,2),idx); %SEÑAL
          t_signalR{1,i} = tiempo(idx_rojo(i,1):idx_rojo(i,2),1); %TIEMPO DE LA SEÑAL
    end
    signal{idx,2} = senial_R; %store LFP
      
end


%% Rayleight

LowCut = bajo / (freqmLFP/2);  % filter
HighCut = alto / (freqmLFP/2);
          
% Phase-locking to each segment
fase_nr_V = cell(nro_tetrodos,max(nro_neuronas));%to store phase Congruent
fase_nr_R = cell(nro_tetrodos,max(nro_neuronas));%to store phase Incongruent

for ind = 1 : nro_tetrodos
    
    nr = neuronas{ind};

 %Congruent   
    for ind2 = 1 : length(idx_verde)
        
        [B,A] = butter(4,[LowCut, HighCut]);
        lfp_filt = filtfilt(B,A,signal{1,1}{1,ind2}); %LFP filtering
        
        % Hilbert and phase
        hilb = hilbert(lfp_filt); %hilbert
        fase = angle(hilb);%phase
        clear hilb
        
        
        int = t_signalV{1,ind2};
        
        aux = cell(length(nr),1);
        for idx = 1 : length(nr)
            N = find(nr{idx}(:,2) > int(1) & nr{idx}(:,2) < int(end));
            nr_int = nr{idx}(N,2);
            aux{idx} = nr_int;
        end
        
        for k = 1 : length(nr)
            for h = 1 : length(aux{k})
                idx_nr = find(aux{k}(h) <= int, 1 ,'first');
                fase_nr_V{ind,k} = [fase_nr_V{ind,k},fase(idx_nr)];
            end
        end
    end
    
 %Incongruent
    for ind2 = 1 : length(idx_rojo)
        
        if isempty(ind2)
            break
        end
        
        
        [B,A] = butter(4,[LowCut, HighCut]);
        lfp_filt1 = filtfilt(B,A,signal{1,2}{1,ind2}); %LFP filtering
        
        % Hilbet y phase
        hilb = hilbert(lfp_filt);
        fase = angle(hilb);
        clear hilb
        
        
        int = t_signalR{1,ind2};
        
        %Find Spikes withing the LFP time interval
        aux = cell(length(nr),1);
        for idx = 1 : length(nr)
            N = find(nr{idx}(:,2) > int(1) & nr{idx}(:,2) < int(end));
            nr_int = nr{idx}(N,2);
            aux{idx} = nr_int;
        end
        
        %find each spike position in the time vector
        for k = 1 : length(nr)
            for h = 1 : length(aux{k})
                idx_nr = find(aux{k}(h) <= int, 1 ,'first');
                fase_nr_R{ind,k} = [fase_nr_R{ind,k},fase(idx_nr)];
            end
        end
    end
end

nro_spk = []; %phases where each spike occurs
fase_V = cell(nro_tetrodos,max(nro_neuronas)); %store phases congruent
fase_R = cell(nro_tetrodos,max(nro_neuronas)); %store phases incongruent
phase_V = []; 
phase_R = [];
for idx = 1 : nro_tetrodos
    
    for idx2 = 1 : nro_neuronas(idx)
        nro_spikesV = length(fase_nr_V{idx,idx2}); %congruent
        nro_spikesR = length(fase_nr_R{idx,idx2}); %incongruent
        aux = [nro_spikesV,nro_spikesR];
        nro_spk = [nro_spk;aux];
        
        if nro_spikesV < nro_spikesR
            fase_V {idx,idx2} = fase_nr_V{idx,idx2}(1:nro_spikesV);
            fase_R {idx,idx2} = fase_nr_R{idx,idx2}(1:nro_spikesV);
            phase_V = [phase_V,fase_nr_V{idx,idx2}(1:nro_spikesV)]; 
            phase_R = [phase_R,fase_nr_R{idx,idx2}(1:nro_spikesV)]; 
        else
            fase_V {idx,idx2} = fase_nr_V{idx,idx2}(1:nro_spikesR);
            fase_R {idx,idx2} = fase_nr_R{idx,idx2}(1:nro_spikesR);
            phase_V = [phase_V,fase_nr_V{idx,idx2}(1:nro_spikesR)]; 
            phase_R = [phase_R,fase_nr_R{idx,idx2}(1:nro_spikesR)];             
        end

    end
end

[p(1),theta(1),rbar(1),delta(1)] = rayleigh(phase_V'); %p=probability; theta=non-weighted circular mean ; rbar=mean resultant length (rbar);  delta = dispersion
[p(2),theta(2),rbar(2),delta(2)] = rayleigh(phase_R'); %p=probability; theta=non-weighted circular mean ; rbar=mean resultant length (rbar);  delta = dispersion
figure(1), rose(phase_V,30),xlabel({['Prob = ',num2str( p(1)),' mod =',num2str( theta(1))],['n = ',num2str(length(phase_V'))]},'FontSize',8,'Fontweight','demi'),title('Roseta de Verdes')
figure(2), rose(phase_R,30),xlabel({['Prob = ',num2str( p(2)),' mod =',num2str( theta(2))],['n = ',num2str(length(phase_R'))]},'FontSize',8,'Fontweight','demi'),title('Roseta de Rojos')

%% Save
filename = strcat(path,nombre_LFP(2:end-4));
save(filename,'senial_R','senial_V','signal', 'fase_nr_V', 'fase_nr_R','t_signalV', 't_signalR', 'fase_V','fase_R','phase_V','phase_R') %,'Rayleigh_todo','fase_acum');