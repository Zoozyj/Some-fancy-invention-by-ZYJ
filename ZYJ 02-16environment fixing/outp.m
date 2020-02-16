classdef sink < handle
    
    properties
        tc
    end

    methods
        
        %% Constructor
        function this = sink(tc)
            if nargin = 1
                this.tc = tc;
            else
                error('Wrong number of input argument for class ue constructor')
            end
        end
        
        
        %% Plot signals output TRP and signal input of UE in time domain
        function plot_signal_per_ue(obj)
        
            % Time domain signal (real part) before precoding for each of 6 UE closest TRP
            figure;
            for i=1:6
                y = signalInTRPs(:,i); 
                subplot(2,3,i);%  hold on; 
                plot(real(y)); 
                grid on;
                axis tight;
            end

            % Time domain signal at UE location (after applied channel) from the 6
            % closest TRP
            figure; 
            for i=1:6
                y = signalOutTRPs(:,1,i); 
                subplot(2,3,i);%  hold on; 
                plot(real(y)); 
                grid on;
                axis tight;
            end
            
        end
        
        %% plot correlations -  simple receiver
        function plot_ue_corr(obj)
            figure; 
            for i=1:6
                x = [0:length(signalInTRPs(:,1,1))-1].'/(30.720e6*8)*3e8;
            %    x = [0:length(signalInTRPs(:,1,1))-1].';
                y = xcorr(signalInTRPs(:,1),signalInTRPs(1:160+144*2+2048*3+2048,i)); 
            %    y = xcorr(signalInTRPs(:,1),signalInTRPs(:,i)); 
                y = abs(y(length(signalInTRPs(:,1)):end));
                subplot(2,3,i);%  hold on; 
            %    plot(x,y); 
                plot(x(1:500),y(1:500)); 
                grid on;
                axis tight;
            end

            figure;
            for i=1:6
                x = [0:length(signalInTRPs(:,1,1))-1].'/(30.720e6*8)*3e8;
            %    y = xcorr(sum(signalOutTRPs,3),signalInTRPs(:,i)); 
                y = xcorr(sum(signalOutTRPs,3),signalInTRPs(1:160+144*2+2048*3+2048,i)); 
                y = y(length(signalInTRPs(:,1)):end);
                subplot(2,3,i);%  hold on; 
                plot(x(1:100),abs(y(1:100)));
                grid on;
                axis tight;
            %    y = xcorr(signalOutTRPs(:,1,i),signalInTRPs(:,i)); 
                y = xcorr(signalOutTRPs(:,1,i),signalInTRPs(1:160+144*2+2048*3+2048,i)); 
                y = y(length(signalInTRPs(:,1)):end);
                hold on;
                plot(x(1:100),abs(y(1:100)),'r--');
            end
        end

        function [LocRate, Fail_Mat, Fail_Pos ]=PlotPosErr(obj, ChanNum,PosBuf,UEPositions,simConfig,fignum,kk,num)
            % Fail_Mat = simConfig.Fail;
            % Fail_Pos = simConfig.FailPos;
            % flag = simConfig.flag;
            % c = 299792458.0;
            Position_Error_All = [];
            Rx_title = [];
            Rx_title1 = [];
            Rx_title2 = [];
            
            for j=1:length(simConfig.flag)
                Position_Error_UEs=[];
                for n=1:simConfig.numUEs
                    
                    estPos             = PosBuf{j}{n};                                             %to read position based on channel number and configuration
                    PosErr             = [];
                    if ~isempty(estPos)
                        PosErr         = [PosErr;sqrt((UEPositions(n,1)-estPos(:,1)).^2+(UEPositions(n,2)-estPos(:,2)).^2)];
                        %                     PosErr         = [PosErr;(sqrt((UEPositions(n,1)-estPos(:,1)).^2+(UEPositions(n,2)-estPos(:,2)).^2))/c];
                    else
                        PosErr         = [PosErr;10000];
                    end
                    Position_Error_UEs = [Position_Error_UEs; PosErr];
                end
                Position_Error_All = [Position_Error_All ,Position_Error_UEs];
            end
            
            % only if the error is less than 500 meter it is a valid detection
            
            Position_Error_Config = [];
            linespec    = {'b-', 'g-', 'r-', 'k*-', 'r*-','mo-','co-','gd-','bh-','y*-','y*-','k-','m-','r-','y-','m-','k-','c-','bh-','y*-','r-'}; % different color for each configuration
            Legend     = {'Original','2 PRS Optimized', '2 PRS manual','3 PRS Manual','3PRS Optim Sequence','3 PRS 3-2-3' ,...
                '3 PRS 2-1-2','3 PRS 3-1-3','4PRS manual','4 PRS optimized','3PRS','2 PRS test','2 PRS ZC-Walsh','ZC-seq-997',...
                '2 PRS Test 2','extra 2 PRS test','Barker','Walsh Hadamard','ZC-seq-63','ZC-seq-839','Barker 8PSK'};
            
            
            for ii=1:length(simConfig.flag)
                % Ind = find(Position_Error_All(:,ii)<500);
                Ind              = find(Position_Error_All(:,ii)<250);
                fail               = find(Position_Error_All(:,ii)>250);
                % fail               = find(Position_Error_All(:,ii)>500);
                failedUEpos = UEPositions(fail,:);
                ii
                %% to be confirmed later where to use it since it needs to be updated for all configurations
                Fail_Mat{ii}= fail;
                Fail_Pos{ii}= failedUEpos;
                LocRate(ii)           = numel(Ind)/simConfig.numUEs;%/ChanNum;
                Position_Error_Config = Position_Error_All(:,ii);
                %% to be added when plotting for 3 PRS and 4 PRS
                % if flag(ii)==1
                %     figure(fignum+num*kk)
                %     hold on
                %     grid on
                %     [posCDF, posError]   = homemade_ecdf(Position_Error_Config(Ind)');     %homemade_ecdf works only for row vectors not column vectors thats why a transpose is needed at input
                %     stairs(posError, posCDF,linespec{ii},'LineWidth',1)
                %     legend(Legend(flag(ii)))
                %     Rx_title = Legend(flag(ii));
                %     figure(fignum+1+num*kk)
                %     hold on
                %     grid on
                %     [posCDF, posError]   = homemade_ecdf(Position_Error_Config(Ind)');
                %     stairs(posError, posCDF,linespec{ii},'LineWidth',1)
                %     legend(Legend(flag(ii)))
                %     Rx_title1 = Legend(flag(ii));
                %     figure(fignum+2+num*kk)
                %     hold on
                %     grid on
                %     [posCDF, posError]   = homemade_ecdf(Position_Error_Config(Ind)');
                %     stairs(posError, posCDF,linespec{ii},'LineWidth',1)
                %      legend(Legend(flag(ii)))
                %      Rx_title2 = Legend(flag(ii));
                % end
                
                if ((flag(ii)<4 && flag(ii)>=1) || flag(ii) == 12 || flag(ii) == 13 || flag(ii) == 14 || flag(ii) == 15 || flag(ii) == 16 || flag(ii) == 17 || flag(ii) == 18 || flag(ii) == 19 || flag(ii) == 20 || flag(ii) == 21)
                    figure(fignum+num*kk)
                    hold on
                    grid on
                    [posCDF, posError]   = homemade_ecdf(Position_Error_Config(Ind)');
                    stairs(posError, posCDF,linespec{flag(ii)},'LineWidth',1)
                    % axis   ([0, 500, 0, 1])
                    axis   ([0, 250, 0, 1])
                    
                    title(['SNR',num2str(kk)])
                    xlabel ('Horizontal error [m]');
                    ylabel ('CDF');
                    Rx_title = [Rx_title, Legend(flag(ii))];
                    legend(Rx_title)
                    savefig(['SNR',num2str(fignum+num*kk)])
                    
                    % elseif ((flag(ii)<9 && flag(ii)>3) ||flag(ii) == 11)
                    % figure(fignum+1+num*kk)
                    % hold on
                    % grid on
                    % [posCDF, posError]   = homemade_ecdf(Position_Error_Config(Ind)');         %homemade_ecdf works only for row vectors not column vectors thats why a transpose is needed at input
                    % stairs(posError, posCDF,linespec{flag(ii)},'LineWidth',1)
                    % title(['SNR',num2str(kk)])
                    % axis   ([0, 250, 0, 1])
                    % xlabel ('Horizontal error [m]');
                    % ylabel ('CDF');
                    % Rx_title1 = [Rx_title1, Legend(flag(ii))];
                    % legend(Rx_title1)
                    % savefig(['SNR',num2str(fignum+num*kk+1)])
                    % elseif flag(ii)>8 && flag(ii)<11
                    % figure(fignum+2+num*kk)
                    % hold on
                    % grid on
                    % [posCDF, posError]   = homemade_ecdf(Position_Error_Config(Ind)');
                    % stairs(posError, posCDF,linespec{flag(ii)},'LineWidth',1)
                    % title(['SNR',num2str(kk)])
                    % axis   ([0, 250, 0, 1])
                    % xlabel ('Horizontal error [m]');
                    % ylabel ('CDF');
                    % Rx_title2 = [Rx_title2, Legend(flag(ii))];
                    % legend(Rx_title2)
                    % savefig(['SNR',num2str(fignum+num*kk+2)])
                    % end
                    
                end
            end
        
end