function [ PRI_array, threshold_min_residue ] = sdif_constant_and_stagger( TOA_array,sdif_threshold_multiplying_factor_array);
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

 display('in main function');
TOA_array_residual=TOA_array;


scanTime=max(TOA_array)-min(TOA_array);
width=scanTime/(length(TOA_array));
maxLevel=floor(length(TOA_array)/5);            %%% minimum 3 pulses are needed to form sequent of PRI --- this is the highest --- minimum 5 pulses .. changed on 1st April
PRI_array=[];
hist_level=1;

%% for testing %%

%maxLevel=2  ;

%% remove above line afterwards
PRI_rejected=[];
maxpeakRejected=[];


for threshold_iterator=1:length(sdif_threshold_multiplying_factor_array)
    sdif_threshold_multiplying_factor=sdif_threshold_multiplying_factor_array(threshold_iterator);
    clear TOA_array_residual;
    clear TOA_array_residual_final_from_sdif;
    TOA_array_residual=TOA_array;
    hist_level=1;
    PRI_array=[];
    PRI_rejected=[];
    maxpeakRejected=[];
    
while hist_level<=maxLevel
    dTOA=[];
%     display('in while loop');
    if length(TOA_array_residual)==0                  %% all PRIs extracted
        display('all PRIs extracted');
        display('PRI_array');
        display(PRI_array);
    break;
    end
    TOA_length=length(TOA_array_residual);
%     display('hist_level');
%     display(hist_level);
    
    for i=1:length(TOA_array_residual)-hist_level
        dTOA(i)=TOA_array_residual(i+hist_level)-TOA_array_residual(i);
    end

length(TOA_array_residual));
% [M C]=hist(dTOA,length(TOA_array_residual));

%% Modified on 30 March
[M C]=hist(dTOA,length(TOA_array));
        maxpeak=max(M);
        display('maxpeak');
        display(maxpeak);
        maxpeakBin=find(M==maxpeak);
        display('maxpeakBin');
        display(maxpeakBin);
        %maxpeakDTOA_index=round(maxpeakBin);                  %%% changed 0n 28 feb
        PRI_indicated=C(maxpeakBin);
        display('PRI_indicated');
        display(PRI_indicated);
        histBinWidth=(max(dTOA)-min(dTOA))/length(TOA_array_residual);
        %errorRange=histBinWidth*2;
        errorRange=0.1*PRI_indicated;

% display('errorRange');
% display(errorRange);
                if(length(maxpeakBin)>1)
                    hist_level=hist_level+1;
                    continue;
                end
                %peakDetectionThreshold=(scanTime*0.02)/(hist_level*PRI_indicated);
                %peakDetectionThreshold=0.9*(length(TOA_array_residual)-hist_level)*exp(-maxpeakBin/(0.05*length(TOA_array_residual)));
%                 peakDetectionThreshold=0.5*(length(TOA_array_residual)-hist_level)*exp(-(0.1*PRI_indicated)/length(TOA_array_residual));

%%%%%% -------
        %% need to change threshold multiplying factor little bit %%
        
       % ---- %%%%%%%

                peakDetectionThreshold=sdif_threshold_multiplying_factor*(length(TOA_array_residual)-hist_level)*exp(-(maxpeakBin)/length(M));
                display('peakDetectionThreshold');
                display(peakDetectionThreshold);
               for i=1:length(PRI_rejected)
                   if PRI_indicated-PRI_rejected(i) < 0.2*PRI_indicated  && length(TOA_array)>=length(TOA_array_residual)
                        maxpeak=maxpeak+maxpeakRejected(i);
                   end
               end
               if maxpeak< peakDetectionThreshold
                   display('peak is small');
                   PRI_rejected=[PRI_rejected PRI_indicated];
                   maxpeakRejected=[maxpeakRejected maxpeak];
               end
            if maxpeak>=peakDetectionThreshold
                display('peak detected for PRI');
                display(PRI_indicated);
            %[res TOA_array_residual]=sequenceSearch(PRI_indicated,errorRange,TOA_array_residual);
            [res TOA_array_residual]=sort_search(PRI_indicated,TOA_array_residual);
            display('in sdif--- res=');
            display(res);
            if strcmp(res,'SUCCESS') && length(dTOA)>0.05*length(TOA_array) && maxpeak> 0.01*length(TOA_array)
                PRI_array=[PRI_array PRI_indicated];
%                 display('TOA_array_residual');
%                 display(TOA_array_residual);
%                 display('dTOA');
%                 display(dTOA);
               hist_level=1;
               display('hist_level reset to 1');
               continue;
            end
            end
            if length(TOA_array_residual)<0.02*length(TOA_array)
                break;
          
            else
                display('hist_level');
                display(hist_level);
                display('length of residual array');
                display(length(TOA_array_residual));
                display('going for next');
            end
                
                
        hist_level=hist_level+1;
end
PRI_array_sdif=PRI_array;
    TOA_array_residual_final_from_sdif=TOA_array_residual;
    mylength=length(TOA_array_residual);
    len(threshold_iterator)=mylength;   
    threshold_criterion(threshold_iterator)=len(threshold_iterator);
    %threshold_criterion(threshold_iterator)=len(threshold_iterator)+length(PRI_array);  %% larger length of PRI_array is desirable
    clear M;
    clear C;
end
min_residual_threshold=min(threshold_criterion);
threshold_array_index=min(find(threshold_criterion==min_residual_threshold));
threshold_min_residue=sdif_threshold_multiplying_factor_array(threshold_array_index);


