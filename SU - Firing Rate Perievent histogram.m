%-----------------
% This code analyse the activity of several Single Units.
% To do so, upload .mat files where the timesamps of the events are store.
% Also, upload the .clu and .res files to use the spikes' timestamps.
% The .clu and .res files are constructed with Neurosuite package.
% Output:
%       1)Index of Excited and inhibited SU in both conditions
%       1)Sorted Perievent histogram (in both conditions)
%-----------------

clear
clc
close all

%% parameters (in seconds)
tpo_preEvento = 20; %time before the event
tpo_postEvento = 20;  %time after the event
bin = 0.1; %size of the bin
orden = 6; %smooth
Neuronas_totales = 0; %counter of SU
tpo_media_pre = -0.2; %time before the event to calculate mean for sorting
tpo_media_post = 0.2; %time after the event to calculate mean for sorting
tpo_preHisto = 2; %limits of the histogram, before
tpo_postHisto = 2; %limits of the histogram, after
tiempo_histo = (-tpo_preHisto:bin:tpo_postHisto-bin); %time vector
tiempo_histo_grande = (-tpo_preEvento:bin:tpo_postEvento-bin); %mean time vector for zscore calculation
z = 1;
%Index of Multiunits to eliminate them
MU = [7;8;16;28;31;34;35;38;45;47;64;69;79;90;92;96;102;104;138;140;149;166;169;171;174;179;190;192;193;194;195;196;215;219;230];

%Selection of Path where all the sessions are store
%For windows
[myPath] = uigetdir('myPath','Seleccione carpeta');

carpeta1 = dir(fullfile(myPath));
carpeta = carpeta1(3:end); %Path for all my sessions


for carp=1:length(carpeta)
     %.mat files with events
    load([myPath,'\', carpeta(carp).name])   

    
    %I keep only the events larger than the median
    largoR = interval_rojo(:,2)-interval_rojo(:,1);
    largoV = interval_verde(:,2)-interval_verde(:,1);
    MedianR = median(largoR,1);
    MedianV = median(largoV,1);
    % for incongruent events
    interval_Rojo = [];
    for i = 1:length(interval_rojo)
        tmp = [];
        if largoR(i) >= MedianR
            tmp = interval_rojo(i,:);
        end
        interval_Rojo = [interval_Rojo;tmp];
    end
    clear largoR MedianR
    % for congruent events
    interval_Verde = [];
    for i = 1:length(interval_verde)
        tmp = [];
        if largoV(i) >= MedianV
            tmp = interval_verde(i,:);
        end
        interval_Verde = [interval_Verde;tmp];
    end
    clear largoV MedianV i 
    
    %%Shuffle of the events
    %Incongruent
	r = randperm(length(interval_Rojo));
    for i = 1: length(interval_Rojo)
        tmp (i,:) = interval_Rojo(r(i),:);
    end
    interval_Rojo = tmp;    
    clear i r tmp
    %Congruent
	r = randperm(length(interval_Verde));
    for i = 1: length(interval_Verde)
        tmp (i,:) = interval_Verde(r(i),:);
    end
    interval_Verde = tmp;
    clear i r tmp
    
    %keep the same number of events across conditions
    x = length(interval_Verde);
    y = length(interval_Rojo);
 
    if x<y
        interval_Rojo = interval_Rojo(1:x,:);
        interval_Verde = interval_Verde(1:x,:);
    else
        interval_Verde(1:y,:);
        interval_Rojo = interval_Rojo(1:y,:);
    end
    clear x y
    
    %Start(ON)of the events
    rojo_on = interval_rojo(:,1)/1000000; %transformation to seconds
    verde_on = interval_verde(:,1)/1000000;
   
    %loop to upload the 
    for TT=1:length(neuronas)
        for neur=1:length(neuronas{TT,1})
            timestamps_neu=neuronas{TT,1}{neur,1}(:,2)/1000000; %transformation to seconds
            Neuronas_totales=Neuronas_totales+1; %counter of SU
                     
            %Construction of histogram locked to ON of the event
            % for incongruent events
            spikes_all_trials_rojo=[];
            spikes_all_trials_verde=[];
            for rojo=1:length(rojo_on)
                pre_local = rojo_on(rojo) - tpo_preEvento;
                post_local = rojo_on(rojo) + tpo_postEvento;
                idx_rojo = find(timestamps_neu>pre_local & timestamps_neu<post_local);
                spikes_trial = timestamps_neu(idx_rojo)-rojo_on(rojo);
                spikes_all_trials_rojo = [spikes_all_trials_rojo;spikes_trial];
            end
            clear spikes_trial
            
            %for the other object event
            for verde=1:length(verde_on)
                pre_local = verde_on(verde) - tpo_preEvento;
                post_local = verde_on(verde) + tpo_postEvento;
                idx_verde = find(timestamps_neu>pre_local & timestamps_neu<post_local);
                spikes_trial= timestamps_neu(idx_verde)-verde_on(verde);
                spikes_all_trials_verde = [spikes_all_trials_verde;spikes_trial];
            end
            clear spikes_trial
            nro_bines=(tpo_preEvento+tpo_postEvento)/bin;
            histo_rojo=(hist(spikes_all_trials_rojo, nro_bines)/bin/length(rojo_on)); 
            histo_verde=(hist(spikes_all_trials_verde, nro_bines)/bin/length(verde_on));
            %calculation of the zscore (-/+ 20 seconds from the ON)
            histo_merge=[histo_rojo histo_verde];
            frec_media=mean(histo_merge);
            desvio=std(histo_merge);
            
            % Restrict between -/+ 2 seconds for graph
            idx_pre_histo=find(tiempo_histo_grande<= -tpo_preHisto, 1,'last');
            idx_post_histo=find(tiempo_histo_grande<= tpo_postHisto, 1,'last');
            histo_rojo_corto=histo_rojo(:,idx_pre_histo:idx_post_histo-1);
            histo_verde_corto=histo_verde(:,idx_pre_histo:idx_post_histo-1);
            
            histo_rojo_zscore = (histo_rojo_corto- frec_media)/desvio; %transformation to zscore        
            histo_verde_zscore = (histo_verde_corto- frec_media)/desvio;
            
            %smooth
            histo_rojo_smooth = smooth(histo_rojo_corto, orden);
            histo_rojo_zscore_smooth = smooth(histo_rojo_zscore, orden);
            
            histo_verde_smooth=smooth(histo_verde_corto, orden);
            histo_verde_zscore_smooth=smooth(histo_verde_zscore, orden);
            
            %% Store
            %  Incongruent
            histo_Rtodos_S(Neuronas_totales,:) = histo_rojo_smooth;
            histo_Rtodos_zscore_S(Neuronas_totales,:) = histo_rojo_zscore_smooth;
            %  Congruent
            histo_Vtodos_S(Neuronas_totales,:) = histo_verde_smooth;
            histo_Vtodos_zscore_S(Neuronas_totales,:) = histo_verde_zscore_smooth;
            
            clear histo_rojo histo_rojo_zscore histo_verde histo_verde_zscore histo_rojo_smooth histo_rojo_zscore_smooth
        end
    end
end

%Elimination of MU Incongruent
histo_Rtodos_S(MU,:)=[];
histo_Rtodos_zscore_S(MU,:)=[];
%Elimination of MU Congruent
histo_Vtodos_S(MU,:)=[];
histo_Vtodos_zscore_S(MU,:)=[];

neuronas_prueba = (1:1:(Neuronas_totales-length(MU)));
Neuronas_totales = Neuronas_totales-length(MU);

%% Detection of excited or inhibited SU
idx_pre_media = find(tiempo_histo >= tpo_media_pre,1, 'first');
idx_post_media = find(tiempo_histo >= tpo_media_post,1, 'first');

% Incongruent
tpo_promediar = histo_Rtodos_zscore_S(:,idx_pre_media:idx_post_media);
promedio = mean(tpo_promediar')';
R_S = promedio; %average of SU in the Incongruent condition
[fr,xr] = ksdensity(R_S,'Function','cdf');

[Exc_rojos,Y] = find(promedio>z); %index of excited SU
[Inh_rojos,Y] = find(promedio<-z); %index of inhibited SU
clear tmp

tmp = histo_Vtodos_zscore_S(Exc_rojos,idx_pre_media:idx_post_media); 
tmp = mean(tmp')';
Exc_rojos = [Exc_rojos,promedio(Exc_rojos),tmp]; %store Excited SU in the other condition
clear tmp

tmp = histo_Vtodos_zscore_S(Inh_rojos,idx_pre_media:idx_post_media); 
tmp = mean(tmp')';
Inh_rojos = [Inh_rojos,promedio(Inh_rojos),tmp]; %store Inhibited SU in the other contion
clear tmp

ROJOS_EXC = Exc_rojos;
ROJOS_INH = Inh_rojos;

%Sorting the histogram according the average in the Incongruent condition
[prom_ordenRZS, idx_mediaRZS] = sort(promedio, 'descend')

histo_Rtodos_zscore_S_Ord = zeros(size(histo_Rtodos_zscore_S));
histo_Vtodos_zscore_S_OrdR = zeros(size(histo_Rtodos_zscore_S));
for orden=1:Neuronas_totales
    histo_Rtodos_zscore_S_Ord(orden,:) = histo_Rtodos_zscore_S(idx_mediaRZS(orden),:);
    histo_Vtodos_zscore_S_OrdR(orden,:) = histo_Vtodos_zscore_S(idx_mediaRZS(orden),:);
end
clear tpo_promediar promedio 


%Congruent condtion
tpo_promediar = histo_Vtodos_zscore_S(:,idx_pre_media:idx_post_media);
promedio = mean(tpo_promediar')';
V_S = promedio;
[fv,xv] = ksdensity(V_S,'Function','cdf');

[Exc_verdes,Y] = find(promedio>z); %index of SU excited in Congruent
[Inh_verdes,Y] = find(promedio<-z); %index of SU inhibited in Congruent
clear tmp

tmp = histo_Rtodos_zscore_S(Exc_verdes,idx_pre_media:idx_post_media); %mean of excited SU in the other condition
tmp = mean(tmp')';
Exc_verdes = [Exc_verdes,promedio(Exc_verdes), tmp]; %store excited SU in congruent

clear tmp
tmp = histo_Rtodos_zscore_S(Inh_verdes,idx_pre_media:idx_post_media);%mean of inhibited SU in the other condition
tmp = mean(tmp')';
Inh_verdes = [Inh_verdes,promedio(Inh_verdes),tmp]; %store inhibited SU in congruent
clear tmp

VERDES_EXC = Exc_verdes;
VERDES_INH = Inh_verdes;

%Sorting the histogram according the average in the Congruent condition
[prom_ordenVZS, idx_mediaVZS]=sort(promedio, 'descend');
histo_Vtodos_zscore_S_Ord=zeros(size(histo_Vtodos_zscore_S));
histo_Rtodos_zscore_S_OrdV=zeros(size(histo_Vtodos_zscore_S));
for orden=1:Neuronas_totales
    histo_Vtodos_zscore_S_Ord(orden,:)=histo_Vtodos_zscore_S(idx_mediaVZS(orden),:);
    histo_Rtodos_zscore_S_OrdV(orden,:)=histo_Rtodos_zscore_S(idx_mediaVZS(orden),:);
end
clear tpo_promediar promedio 
