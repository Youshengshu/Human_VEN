% Note:AP threshold不是固定的值，而是找slope的拐点
cd ("G:\VEN project_KW\Data analysis\Patch-seq mutimodal\Etype\matfile_AP_Para\matfile_FI")

clear all;close all;clc
dbstop if error
WantState=1;%1/2：分析with/without holding状态下的数据
Num=6;%选择大部分细胞都会达到的一个AP Num
paraRecord=[];

%% cluster infomation
totaldata=readmatrix('clusteringOutput_using.csv', 'OutputType', 'string');
total_filename=totaldata(:,1);
etype=str2double(totaldata(:,2));% 4 cluster or 3cluster
for kk=2%:5
    idx=[];
    idx=find(ismember(etype,kk)==1);
    Fidx{kk}=total_filename(idx);
        % Fidx{kk}=dir('DATA_FIcurve*.mat');
    
    for k=1:length(Fidx{kk})
                    % filename=Fidx{kk}(k).name;
        filename=Fidx{kk}(k);
        filename=convertStringsToChars(filename);
        load(filename)
        disp([num2str(k),'/',num2str(length(Fidx{kk})),'--',filename])
        dt=tspanV(2)-tspanV(1);
        
        
        %% extract each AP waveform to calculate HW AHP slope change within a AP train
        Wantidx=find(ismember(spkFreq,[Num*2:2:Num*2+20]));
        [~,iidx]=min(StimAmp(Wantidx));


      if ~isempty(Wantidx(iidx))
        for i=Wantidx(iidx)
            pretime=0.1;
            AP_Delay=spkTime{i}(1)-pretime;
            ISI=diff(spkTime{i});
            Fre_median=median(1./ISI);
            ISI=ISI(find(isoutlier(ISI,"grubbs")==0));
            ISIvari=std(ISI);
            AP_Freq1st=1./(spkTime{i}(2)-spkTime{i}(1));
            
            
            data=waveformV(:,i);
            [spkwaveform,ttspanspk,HW, AHP, Threshold, Amp,max_slope,min_slope,...
                HWidx,AHPidx, Thresholdidx,max_slopeidx,min_slopeidx] = slowAP_Statistic(data,spkTime{i},tspanV,1:length(spkTime{i}));
            spkwaveform_temp{i}=spkwaveform;
            AHP=AHP(~isnan(AHP));
            AHP(find(AHP==0))=[];
            HW=HW(~isnan(HW));
            Threshold=Threshold(~isnan(Threshold));
            max_slope=max_slope(~isnan(max_slope));
            min_slope=min_slope(~isnan(min_slope));
            Amp=Amp(~isnan(Amp));
            
            HW_Adap=HW(end-1)./HW(1);
            Threshold_Adap=Threshold(end-2)./Threshold(1);
            %AHP_Adap(i)=(AHP(end-1)-AHP(1))/AHP(1);
            AHP_Adap=mean(AHP(1:end-1));
            maxSlope_Adap=max_slope(end)./max_slope(1);
            minSlope_Adap=min_slope(end)./min_slope(1);
            Amp_Adap=Amp(end)./Amp(1);
            Ampvari=std(Amp);

%             HW_Adap=HW(2)./HW(1);
%             Threshold_Adap=Threshold(2)./Threshold(1);
%             %AHP_Adap(i)=(AHP(end-1)-AHP(1))/AHP(1);
%             AHP_Adap=mean(AHP(1:2));
%             maxSlope_Adap=max_slope(2)./max_slope(1);
%             minSlope_Adap=min_slope(2)./min_slope(1);
%             Amp_Adap=Amp(end)./Amp(1);
%             Ampvari=std(Amp);
%             HW_vari=std(HW);

          
            
            figure(1),clf
            for jj=1:size(spkwaveform,2)
                subplot(2,1,1)
                plot(tspanV,data),hold on
                if ~isnan(HWidx(jj,:))
                    plot(ttspanspk(HWidx(jj,:))+spkTime{i}(jj)-pretime,spkwaveform(HWidx(jj,:),jj),'ro')
                    hold on
                end
                if ~isnan(AHPidx(jj));
                    plot(ttspanspk(AHPidx(jj))+spkTime{i}(jj)-pretime,spkwaveform(AHPidx(jj),jj),'ko')
                    hold on
                    plot(ttspanspk(Thresholdidx(jj))+spkTime{i}(jj)-pretime,spkwaveform(Thresholdidx(jj),jj),'co')
                    
                end
                subplot(2,1,2)
                slope=[nan;diff(waveformV(:,i))/(dt*1000)];
                plot(tspanV,slope','DisplayName','waveformV')                
                box off
            end
            
            figure(2)
            subplot(6,8,k)
            plot(tspanV,data),hold on
            title(filename(14:end), 'Interpreter', 'none')
            ylim([-80,60])
            box off
        end
        
        paraRecord(k,1)=AP_Delay;
        paraRecord(k,2)=AP_Freq1st;
        paraRecord(k,3)=ISIvari;       
        paraRecord(k,4)=Amp_Adap;
        paraRecord(k,5)=Ampvari;
        paraRecord(k,6)=Fre_median;
        paraRecord(k,7)=HW_Adap;
        paraRecord(k,8)=Threshold_Adap;
        paraRecord(k,9)=AHP_Adap;
        paraRecord(k,10)=maxSlope_Adap;
        paraRecord(k,11)=minSlope_Adap;
      end
                % pause
        
    end
end
paramName={'AP_Delay','AP_Freq1st','ISIvari','Amp_Adap','Ampvari','Fre_median','HW_Adap',...
    'Threshold_Adap','AHP_Adap','maxSlope_Adap','minSlope_Adap'};

save('Results_FiringPattern.mat','Fidx',...
    'paraRecord','ttspanspk','paramName')






