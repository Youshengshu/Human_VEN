clear all;clc
close all;
dbstop if error
% dbstop at 45 in AP_Statistic

cd ("G:\VEN project_KW\Data analysis\Patch-seq mutimodal\Etype\matfile_AP_Para\matfile_FI")

WantState=1;%1/2：分析with/without holding状态下的数据
Num=6;%选择大部分细胞都会达到的一个AP Num
spkparamRecord=[];
spkwaveformRecord=[];

for kk=1
        Fidx{kk}=dir('DATA_FIcurve*.mat');
    for k=1:length(Fidx{kk})
%         filename=Fidx(k).name;
        filename=Fidx{kk}(k).name;
       
        load(filename)
        disp([num2str(k),'/',num2str(length(Fidx{kk})),'--',filename])
        dt=tspanV(2)-tspanV(1);
        
        stateIdx{1}=find(round(Vrest)<=-60);%筛选Vrest<-60mV的trial
        stateIdx{2}=find(abs(round(Ihold))<=2);%筛选without holding的trial
        %% ------------- spike train para---------------------------
        uniqueStim=unique(StimAmp);
        if ~isempty(uniqueStim)
            StimAmp_act=zeros(size(StimAmp))';
            spkFreq_act=zeros(size(StimAmp))';
            for i=1:length(uniqueStim)
                idx=find(StimAmp==uniqueStim(i));
                StimAmp_act(i)=mean(StimAmp(idx));
                spkFreq_act(i)=round(mean(spkFreq(idx)));
            end
            
            Fun=@(a,b,x) max([zeros(length(x),1),a.*x(:)+b],[],2);
            fo= fitoptions('Method','NonlinearLeastSquares',...
                'Lower',[0.001,-1000],...
                'Upper',[10,1000], ...
                'Startpoint',[0.1,-10]);
            
            [tempx,idx]=sort(StimAmp_act);
            tempy=spkFreq_act(idx);
            tempy(end-2:end)=[];tempx(end-2:end)=[];
            idx=find(tempy<=25);
            
            x=[0;tempx(idx)];
            y=[0;tempy(idx)];
            
            [Model,~]=fit(x(:),y(:),fittype(Fun),fo);
            xx=linspace(0,x(end),100);
            yy=feval(Model,xx);
            
            Fre_threshold=max(spkFreq(find(StimAmp==min(rheobase))));
            maxFre=max(spkFreq);
            
%             FIcurveParamRecord(k,:)=[-(Model.b/Model.a),Model.a,Fre_threshold,maxFre];
        end
        figure(1),clf
        plot(StimAmp_act,spkFreq_act,'o','markersize',6),hold on
        plot(xx,yy,'m')
        xlim([0,400])
        ylim([0,30])
        box off
        pause(0.1)

        %% ------------- Bin spkshape,1~5 AP first spike ---------------------------
        Wantidx=find(ismember(spkFreq,[Num*2:2:Num*2+15]));
        [~,iidx]=min(StimAmp(Wantidx));
        
        if ~isempty(Wantidx)
            Wantidx=Wantidx(iidx);
            spkwaveform_temp=[];
            HWtemp=[];
            AHPtemp=[];
            Thresholdtemp=[];
            Amptemp=[];
            max_slopetemp=[];
            min_slopetemp=[];
            slope=[];
            
            for jj=1%:length(spkTime{Wantidx})
                data=waveformV(:,Wantidx);
                [spkwaveform,ttspanspk,HW, AHP, Threshold, Amp,AIS_maxslope,Soma_maxslope,max_slope,min_slope...
                    HWidx,AHPidx, Thresholdidx,min_slope2idx] = AP_Statistic(data,spkTime{Wantidx}(jj),tspanV,1);
                
                
                spkwaveform_temp(:,jj)=spkwaveform;
                HWtemp(jj,:)=HW;
                AHPtemp(jj)=AHP;
                Thresholdtemp(jj)=Threshold;
                Amptemp(jj)=Amp;
                max_slopetemp(jj)=max_slope;
                min_slopetemp(jj)=min_slope;
                AIS_slope(jj)=AIS_maxslope;
                Soma_slope(jj)=Soma_maxslope;
                
                figure(100),clf
                plot(data)
                
                figure(200)
                H=[];
                %                 H(1)=subplot(2,length(spkTime{Wantidx}),jj);
                H(1)=subplot(3,1,1);
                plot(ttspanspk,spkwaveform_temp(:,jj)),hold on
                plot(ttspanspk(HWidx),spkwaveform_temp(HWidx,jj),'ro')
                plot(ttspanspk(AHPidx),spkwaveform_temp(AHPidx,jj),'ko')
                plot(ttspanspk(Thresholdidx),spkwaveform_temp(Thresholdidx,jj),'co')
                ylim(H(1),[-70,50])
                box off
                
                %                 H(2)=subplot(2,length(spkTime{Wantidx}),jj+length(spkTime{Wantidx}));
                H(2)=subplot(3,1,2);
                slope(:,jj)=[nan;diff(spkwaveform_temp(:,jj))/(dt*1000)];
                slope2=[nan;diff(slope(:,jj))/(dt*1000)]; slope2=smooth(slope2,5);
                slope3=[nan;diff(slope2)/(dt*1000*100)];slope3=smooth(slope3,5);
                box off
                                plot(ttspanspk,slope,'k'),hold on
                plot(ttspanspk,slope2','b'),hold on
                plot(ttspanspk(min_slope2idx),slope2(min_slope2idx),'ro')
%                 ylim(H(2),[-800,800])
                box off
                %                 plot(ttspanspk,spkwaveform),hold on
                %                 plot(ttspanspk(min_slope2idx),spkwaveform(min_slope2idx),'ro')
                H(3)=subplot(3,1,3);
                plot(spkwaveform_temp(:,jj),slope(:,jj)),hold on
%                 ylim(H(3),[-150,250])
                xlim(H(3),[-60,60])
                drawnow
             
            end
           
            
            spkparamRecord(k,1)=mean(HWtemp);
            spkparamRecord(k,2)=mean(AHPtemp);
            spkparamRecord(k,3)=mean(Thresholdtemp);
            spkparamRecord(k,4)=mean(Amptemp);
            spkparamRecord(k,5)=Vrest(Wantidx);
            spkparamRecord(k,6)=min(rheobase);
            spkparamRecord(k,7)=Fre_threshold;
            spkparamRecord(k,8)=Soma_slope(1);
            spkparamRecord(k,9)=max_slopetemp(1);
            spkparamRecord(k,10)=min_slopetemp(1);
            
            spkwaveformRecord(k,:)=spkwaveform_temp(:,1);%第一个spike waveform
            
            ParamName={'HW', 'AHP', 'Threshold', 'Amp','Vm','Rheobase','Fre_threshold','Soma_slope','max slope','min slope'};
            
            save('Results_spkWaveform.mat','Fidx',...
                'tspanV','spkwaveformRecord','spkparamRecord','ttspanspk','ParamName')
        else
            ttspanspk=[];
            spkparamRecord(k,1)=nan;
            spkparamRecord(k,2)=nan;
            spkparamRecord(k,3)=nan;
            spkparamRecord(k,4)=nan;
            spkparamRecord(k,5)=nan;
            spkparamRecord(k,6)=nan;
            spkparamRecord(k,7)=nan;
            spkparamRecord(k,8)=nan;
            spkparamRecord(k,9)=nan;
            spkparamRecord(k,10)=nan;
            spkwaveformRecord(k,:)=nan(1376,1);
            
        end
%         pause
    end
end


