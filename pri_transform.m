function [absD, Tao_K_t,peaks,peakBins_clustered,threshold]=pri_transform(tao_max, TOA_array)

display('in pri_transform');
N=length(TOA_array);

tao_min=0;


scanTime= max(TOA_array)-min(TOA_array);

K=20*length(TOA_array)+1;
b=(tao_max-tao_min)/K;

D=zeros(K,1); %%% doubt 

for k=1:K                           %%% Loop 1
Tao_K(k)=(k-1/2)*b+tao_min;
Tao_K_low(k)=Tao_K(k)-b/2;
Tao_K_up(k)=Tao_K(k)+b/2;
end

   
    n=2;
   
    
    while(n<N)                      %%% Loop 2   --- confirm loop exit criterion
%         display('in Loop 1');
         m=n-1;
        while(m>=1)                 %%% Loop 3
           %display('in Loop 2');

            tao=TOA_array(n)-TOA_array(m);      

            if(tao<=tao_min)
                m=m-1;
                continue;
            end
                if(tao>tao_max)
%                     n=n+1;
                    break;
                end
  
                
                    for k=1:K

                        if(tao>Tao_K_low(k) && tao<=Tao_K_up(k))
                        %%% Processing to be done
                            %% update Dk
                            D(k)=D(k)+exp((1i*2*pi*TOA_array(n))/(tao)); %% 1i is same as i or j --- t(n)-t(m)= tao


                        end
                    end
%                     display('exit loop 3')
%                 else
%                 m=m-1;
                %end
                
%             else
%             n=n+1;    
%             end

         m=m-1;
        end
        n=n+1;
    end
    absD=abs(D);
    Tao_K_t=Tao_K';
    plot(Tao_K_t, absD);
    xlabel('PRI');
    ylabel(' absolute value of PRI Transform');
    
    %% thresholding
    scanTime=max(TOA_array)-min(TOA_array);
    for k=1:length(D)

	threshold(k)=0.3*scanTime/Tao_K(k);         %%% as mix of all 3 PRI types is 								%present
    end
    

    peaks=[];
    peakBins=[];
    for k=1: length(D)
            if absD(k)>threshold(k) && absD(k)>0.05*length(TOA_array)   %% additional percentage condition  

                peaks=[peaks absD(k)];
                peakBins=[peakBins Tao_K(k)];

            end
    end
    
    plot(Tao_K_t, absD);
    hold
    plot(Tao_K,threshold,'-.');
    xlabel('PRI');
    ylabel('absolute value of PRI Transform');
    legend('PRI Transform','threshold function');
    %% added 21 March
    
    peakBins_clustered=[];
    for i=1: length(peakBins)
        if i==1
            peakBins_clustered=peakBins(i);
        elseif abs(peakBins(i)-peakBins(i-1))<=0.01*(peakBins(i)+peakBins(i-1))
            %display('clustering');
            peakBins_clustered=peakBins_clustered(peakBins_clustered~=peakBins(i-1));
            peakBins_clustered=[peakBins_clustered peakBins(find(peaks==max(peaks(i),peaks(i-1))))];
        else
            %display('in else');
            peakBins_clustered=[peakBins_clustered peakBins(i)];
        end
    end
    