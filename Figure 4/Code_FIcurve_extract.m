clear all;clc
close all;

Axon_data=readtable('Axon_origion_FI.xlsx');
sample=Axon_data.filename;
for k=13:length(sample)
    
        filename=[cell2mat(Axon_data.filename(k)),'.mat'];        
        cellType=Axon_data.cell_type(k);
        Dist=Axon_data.Axon_distance(k);        

        %% loading data
        if ~exist(filename)
            continue
        end
        load(filename)
        disp([num2str(k),'--',num2str(length(sample)),'--',filename])
        %% ------------- reading data ---------------------------
        clear(('*Ch31*'));
        varia=who('*Ch*');
        
        for i=1:length(varia)
            if isfield(eval(varia{i}),'units')
                str=eval([varia{i},'.units']);
                
                if ismember('pA',str)
                    Im=eval([varia{i} '.values']);
                    file=eval([varia{i} '.interval']);
                    dtI=eval([varia{i} '.interval']);                    
                end
                
                if ismember('mV',str)
                    Vm=eval([varia{i} '.values']);
                    dtV=eval([varia{i} '.interval']);
                end
            end
        end
        % 插值至50K Hz采样率
        if round(1/(dtI*1000))==25
            tspan=linspace(1,size(Im,1)*dtI,size(Im,1));
            ttspan=linspace(1,size(Im,1)*dtI,2*size(Im,1));
            Im=interp1(tspan,Im,ttspan,'spline');
            Im=Im';
            dtI=dtI/2;
        end
        if round(1/(dtV*1000))==25
            tspan=linspace(1,size(Vm,1)*dtV,size(Vm,1));
            ttspan=linspace(1,size(Vm,1)*dtV,2*size(Vm,1));
            Vm= interp1(tspan,Vm,ttspan,'spline');
            Vm=Vm';
            dtV=dtV/2;
        end
%         pause
        clear (varia{:})
        %%  finding key time point
        pretime=0.1;
        postime=0.8;
        TimeStim = eventedgefinder(Im,5,1/dtI,0.3,0.5,1,0).*dtI;%holding 情况下找的，只能找到Im大于0的stim        
        TimeStim=TimeStim(:,1);
        disp(length(TimeStim))
        Ihold=[];
        for i=1:length(TimeStim)
            baseidxI=round(TimeStim(i)/dtI-pretime/dtI:TimeStim(i)/dtI);
            Ihold=[Ihold,mean(Im(baseidxI))];
        end
        %% waveform extract
        waveformV=[];
        waveformI=[];
        Vrest=[];
        

        TimeStim_hold=eventedgefinder(Im-min(Ihold),5,1/dtI,0.3,0.5,1,0).*dtI;%without holding current找stim
        Dura=round(range(TimeStim_hold,2)*100)./100;
        TimeStim=unique(roundn(union(TimeStim,TimeStim_hold(:,1)),-3));  %with/without holding current找到的stim求并集
        
        TimeStim(find(isnan(TimeStim)))=[];
        TimeStim(find((TimeStim==0)))=[];
        
        Ihold=[];
        
        for i=1:length(TimeStim)
            dataidxV=round(TimeStim(i)/dtV-pretime/dtV:TimeStim(i)/dtV+postime/dtV);
            baseidxV=round(TimeStim(i)/dtV-pretime/dtV:TimeStim(i)/dtV);
            dataidxI=round(TimeStim(i)/dtI-pretime/dtI:TimeStim(i)/dtI+postime/dtI);
            baseidxI=round(TimeStim(i)/dtI-pretime/dtI:TimeStim(i)/dtI);
            
            Vrest=[Vrest,mean(Vm(baseidxV))];
            Ihold=[Ihold,mean(Im(baseidxI))];
            waveformV=[waveformV,Vm(dataidxV)];
            %             waveformI=[waveformI,Im(dataidxI)-mean(Im(baseidxI))];
            waveformI=[waveformI,Im(dataidxI)];
        end
        tspanV=linspace(-pretime,postime,size(waveformV,1));
        tspanI=linspace(-pretime,postime,size(waveformI,1));
        
    
        
        %% make same necessary Tag
        % spk shape
        spkFreq=[];
        spkTime=[];
        rheobase=[];
        
        for i=1:size(waveformV,2)
            peaktime=peakfinder(waveformV(:,i)-Vrest(i),10,-10, 1, 0).*dtV;
            wrongidx=[];
            for j=1:length(peaktime)
                idx=round(peaktime(j)/dtV-0.01/dtV:peaktime(j)/dtV);
                temp=min(waveformV(idx,i));
                if temp>waveformV(idx(end),i)-18;
                    wrongidx=[wrongidx;j];
                end
                if waveformV(idx(end),i)<-20&peaktime(j)>0.5+pretime;
                    wrongidx=[wrongidx;j];
                end
            end
            peaktime(wrongidx)=[];
            
            if ~isempty(peaktime);peaktime(peaktime>0.51+pretime|peaktime<pretime+0.0005)=[];end
            %             idx=find(diff(peaktime)>0.001);
            %             idx=[1;idx(:)+1];
            %             peaktime=peaktime(idx);
            
            %             StimAmp=mean(waveformI(tspanI>0.1&tspanI<0.4,:))-mean(waveformI(tspanI<0,:));
%             StimAmp=round(mean(waveformI(tspanI>0.1&tspanI<0.4,:))./10)*10;%%不考虑holding current
            StimAmp=round(mean(waveformI(tspanI>0.1&tspanI<0.4,:)-mean(waveformI(tspanI>-0.1&tspanI<0,:)))./10)*10;%%考虑holding current
            spkFreq=[spkFreq;length(peaktime)/0.5];
            spkTime{i}=peaktime;
            
            figure(1),clf
            subplot(2,1,1)
            plot(tspanI,waveformI(:,i)),hold on
            title(filename, 'Interpreter', 'none')
            ylim([-250,100])
            xlim([-0.1,0.9])
            box off
            subplot(2,1,2)
            plot(tspanV,waveformV(:,i)),hold on
            plot(tspanV(round(peaktime/dtV)),waveformV(round(peaktime/dtV),i),'ro')
            ylim([-75,60])
            xlim([-0.1,0.9])
            box off
            pause
        end
        
        temp=waveformV-repmat(Vrest,size(waveformV,1),1);
        VAHP=trapz(temp(tspanV>0.51&tspanV<0.8,:),1);
        % %
        %         Vposi=mean(waveformV(tspanV>0.4&tspanV<Dura,:))-mean(waveformV(tspanV<0,:));
        %         Vposi(spkFreq~=0)=nan;
        
        % rheobase  notice that this method not suitable for cVENs
%         Idx=find(spkFreq==min(spkFreq(spkFreq>0)));
%         rheobase=min(StimAmp(Idx));

        Idx=find(spkFreq>0);
        rheobase=min(StimAmp(Idx));
        
        save(['DATA_FIcurve_',filename],'cellType','Dist','tspanV','tspanI','waveformV','waveformI','TimeStim','Dura',...
            'Vrest','Ihold','VAHP','spkFreq','spkTime','StimAmp','rheobase')
%                 pause
        
        %         figure(1),clf
        %         plot(StimAmp,VAHP,'ro')
        %         drawnow
    end
   

