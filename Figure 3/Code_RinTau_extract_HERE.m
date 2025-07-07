clear all;clc
cd ("G:\VEN project_KW\Data analysis\Patch-seq mutimodal\Etype\matfile_AP_Para\matfile_FI")
Fidx{1}=dir('KW_20230412_S3C4.mat');

%Fidx{3}=dir('three*.mat');
%%
for kk=1:length(Fidx)
    
    for k=1:length(Fidx{kk})
        %% loading data
        filename=Fidx{kk}(k).name;
        load(filename)
        disp([num2str(k),'/',num2str(length(Fidx{kk})),'--',num2str(kk),'/',num2str(length(Fidx)),'--',filename])
        %% ------------- reading data ---------------------------
        clear(('*Ch31*'));
        varia=who('*Ch*');
        
        for i=1:length(varia)
            if isfield(eval(varia{i}),'units')
                str=eval([varia{i},'.units']);
                
                if ismember('pA',str)
                    Im=eval([varia{i} '.values']);
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
        
        clear (varia{:})
        
        %%  finding key time point
        TimeStim=[];
        Im=Im-mean(Im(1:100));
        pretime=0.1;
        postime=0.8;
        TimeStim = eventedgefinder(Im,-5,1/dtI,0.4,0.5,-1,0).*dtI;%holding 情况下找的，只能找到Im大于0的stim
        TimeStim=TimeStim(:,1);
%         if ~isnan(TimeStim)
%             TimeStim=TimeStim(:,1);
%             TimeStim(find(round(TimeStim)==0))=[];
%             else
%             TimeStim=nan;
%         end
%             % 找到holding current(利用向上的刺激trial)
%             pretime=0.1;
%             postime=0.8;
%             TimeStimup = eventedgefinder(Im,5,1/dtI,0.3,0.5,1,0).*dtI;%holding 情况下找的，只能找到Im大于0的stim
%             Ihold=[];
%             for i=1:length(TimeStimup)
%                 baseidxI=round(TimeStimup(i)/dtI-pretime/dtI:TimeStimup(i)/dtI);
%                 Ihold=[Ihold,mean(Im(baseidxI))];
%             end
%             %% waveform extract
%             TimeStim_hold= eventedgefinder(Im-min(Ihold),-5,1/dtI,0.3,0.5,-1,0).*dtI;%去掉holding current找stim
%             TimeStim_hold=TimeStim_hold(:,1);
%             TimeStim=unique(roundn(union(TimeStim,TimeStim_hold),-3));  %去掉holding current找到的stim和不去掉找到的stim求并集
%             TimeStim(find(TimeStim==0))=[];
            disp(length(TimeStim))
            if ~isnan(TimeStim)
                Ihold=[];
                waveformV=[];
                waveformI=[];
                Vrest=[];
                StimAmp=[];
                for i=1:length(TimeStim)
                    dataidxV=round(TimeStim(i)/dtV-pretime/dtV:TimeStim(i)/dtV+postime/dtV);
                    baseidxV=round(TimeStim(i)/dtV-pretime/dtV:TimeStim(i)/dtV);
                    dataidxI=round(TimeStim(i)/dtI-pretime/dtI:TimeStim(i)/dtI+postime/dtI);
                    baseidxI=round(TimeStim(i)/dtI-pretime/dtI:TimeStim(i)/dtI);
                    
                    Vrest(i)=median(Vm(baseidxV));
                    Ihold(i)=round(mean(Im(baseidxI))/10)*10;                   
                    waveformV=[waveformV,Vm(dataidxV)];
                    waveformI=[waveformI,Im(dataidxI)-median(Im(baseidxI))];
                end
                tspanV=linspace(-pretime,postime,size(waveformV,1));
                tspanI=linspace(-pretime,postime,size(waveformI,1));
                
                % exclude wrong trace(in case of epileptic network activity)
                wrongidx=[];
                temp=mean(waveformI(tspanV>0&tspanV<0.5,:),1);
                wrongidx=[wrongidx;find(temp>-50)'];
                
                temp2=mean(waveformI(tspanV>0.5,:),1);
                wrongidx=[wrongidx;find(round(temp2-temp)==0)'];
                TimeStim(wrongidx)=[];
                waveformI(:,wrongidx)=[];
                waveformV(:,wrongidx)=[];
                Vrest(wrongidx)=[];
                
                %找到reboud AP
                reboundidx=[];
                for i=1:length(TimeStim)
                    peaktime=peakfinder(waveformV(:,i),20,20, 1, 0).*dtV;
                    peaktimeBase=peakfinder(waveformV(tspanV<0,i),20,20, 1, 0).*dtV;
                    if ~isempty(peaktime)&isempty(peaktimeBase)
                        reboundAP=length(peaktime);
                    else
                        reboundAP=nan;
                    end
                    reboundidx=[reboundidx;reboundAP];
                end
                
                
                %% make same necessary Tag
                % spk shape
                figure(1),clf
                RinRecord=[];
                tauRecord=[];
                Cm=[];
                sagRatio=[];
                Max=[];
                
                Fun=@(tau,A,B,x) A.*exp(-x./tau)+B;
                fo= fitoptions('Method','NonlinearLeastSquares',...
                    'Lower',[0,1,-100],...
                    'Upper',[10,100,10], ...
                    'Startpoint',[0.1,15,-70]);
                
                figure(1),clf
                for i=1:size(waveformV,2)
                    x=tspanV(tspanV>0.001&tspanV<0.05);
                    y=waveformV(tspanV>0.001&tspanV<0.05,i);
                    [model,gof]=fit(x(:),y(:),Fun,fo);
                    yy=feval(model,x);
                    gof.rsquare
                    if gof.rsquare>0.8
                        tau=(model.tau).*1000;%单位为ms                        
                    else                        
                        tau=nan;
                    end
                    tauRecord=[tauRecord;tau];
                    Max(i)=Vrest(i)-min(waveformV(tspanV>0&tspanV<0.2,i));
                    SS(i)=Vrest(i)-median(waveformV(tspanV>0.4&tspanV<0.45,i));
                    sagRatio=[sagRatio;(Max(i)-SS(i))./Max(i)];
                    
                    SSI=-median(waveformI(tspanI>0.4&tspanI<0.45,i));
                    Rin=Max(i)./SSI*1000;%单位为兆欧
                    RinRecord=[RinRecord;Rin];
                    Cm=[Cm;tau./Rin*1000];%单位为pF
                    
                    subplot(211)
                    plot(tspanV,waveformI(:,i)),hold on
                    subplot(212)
                    plot(tspanV,waveformV(:,i)),hold on
                    plot(x,yy,'k')
                    drawnow
                    %             pause
                end
                
                % StimAmp
                StimAmp=mean(waveformI(tspanI>0.1&tspanI<0.4,:))-mean(waveformI(tspanI<0,:));
                StimAmp=round(StimAmp/10)*10;
                
                save(['DATA_RinTau_',filename],'tspanV','tspanI','waveformV','waveformI','TimeStim',...
                    'StimAmp','Vrest','Max','Ihold','RinRecord','Cm','tauRecord','sagRatio','reboundidx')
                %         pause
            end
        
        
    end
end