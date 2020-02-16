classdef sink < handle
    
    properties
        tc
    end

    methods
        
        %% Constructor
        function this = sink(tc)
            if nargin == 1
                this.tc = tc;
            else
                error('Wrong number of input argument for class ue constructor')
            end
        end
        
        
        %% Plot signals output TRP and signal input of UE in time domain
%         function plot_signal_per_ue(obj)
%         
%             % Time domain signal (real part) before precoding for each of 6 UE closest TRP
%             figure;
%             for i=1:6
%                 y = signalInTRPs(:,i); 
%                 subplot(2,3,i);%  hold on; 
%                 plot(real(y)); 
%                 grid on;
%                 axis tight;
%             end
% 
%             % Time domain signal at UE location (after applied channel) from the 6
%             % closest TRP
%             figure; 
%             for i=1:6
%                 y = signalOutTRPs(:,1,i); 
%                 subplot(2,3,i);%  hold on; 
%                 plot(real(y)); 
%                 grid on;
%                 axis tight;
%             end
%             
%         end
%         
%         %% plot correlations -  simple receiver
%         function plot_ue_corr(obj)
%             figure; 
%             for i=1:6
%                 x = [0:length(signalInTRPs(:,1,1))-1].'/(30.720e6*8)*3e8;
%             %    x = [0:length(signalInTRPs(:,1,1))-1].';
%                 y = xcorr(signalInTRPs(:,1),signalInTRPs(1:160+144*2+2048*3+2048,i)); 
%             %    y = xcorr(signalInTRPs(:,1),signalInTRPs(:,i)); 
%                 y = abs(y(length(signalInTRPs(:,1)):end));
%                 subplot(2,3,i);%  hold on; 
%             %    plot(x,y); 
%                 plot(x(1:500),y(1:500)); 
%                 grid on;
%                 axis tight;
%             end
% 
%             figure;
%             for i=1:6
%                 x = [0:length(signalInTRPs(:,1,1))-1].'/(30.720e6*8)*3e8;
%             %    y = xcorr(sum(signalOutTRPs,3),signalInTRPs(:,i)); 
%                 y = xcorr(sum(signalOutTRPs,3),signalInTRPs(1:160+144*2+2048*3+2048,i)); 
%                 y = y(length(signalInTRPs(:,1)):end);
%                 subplot(2,3,i);%  hold on; 
%                 plot(x(1:100),abs(y(1:100)));
%                 grid on;
%                 axis tight;
%             %    y = xcorr(signalOutTRPs(:,1,i),signalInTRPs(:,i)); 
%                 y = xcorr(signalOutTRPs(:,1,i),signalInTRPs(1:160+144*2+2048*3+2048,i)); 
%                 y = y(length(signalInTRPs(:,1)):end);
%                 hold on;
%                 plot(x(1:100),abs(y(1:100)),'r--');
%             end
%         end

function [pos_err] = analyze_result(obj, ue_pos, est_coord) 
    pos_err=[];
    for n = obj.tc.env.ue_nbr
        pos_err = [pos_err sqrt((ue_pos(n,1) - est_coord(n,1)).^2 + ...
                                (ue_pos(n,2) - est_coord(n,2)).^2 + ...
                                (ue_pos(n,3) - est_coord(n,3)).^2)];
    end
end

function [v_f,v_x] = homemade_ecdf(v_data)
    
    nb_data = numel(v_data);

    v_sorted_data = sort(v_data);
    v_unique_data = unique(v_data);
    nb_unique_data = numel(v_unique_data);
    v_data_ecdf = zeros(1,nb_unique_data);
    for index = 1:nb_unique_data
        current_data = v_unique_data(index);
        v_data_ecdf(index) = sum(v_sorted_data <= current_data)/nb_data;
    end

    v_x = [v_unique_data(1) v_unique_data];
    v_f = [0 v_data_ecdf];

end

function [LocRate, Fail_Mat, Fail_Pos ] = plot_pos_err(obj, ChanNum, PosBuf, UEPositions, simConfig, fignum, kk, num)
            Position_Error_All = [];
            Rx_title = [];
            Rx_title1 = [];
            Rx_title2 = [];

            Position_Error_Config = [];
%             linespec    = {'b-', 'g-', 'r-', 'k*-', 'r*-','mo-','co-','gd-','bh-','y*-','y*-','k-','m-','r-','y-','m-','k-','c-','bh-','y*-','r-'}; % different color for each configuration
%             Legend     = {'Original','2 PRS Optimized', '2 PRS manual','3 PRS Manual','3PRS Optim Sequence','3 PRS 3-2-3' ,...
%                 '3 PRS 2-1-2','3 PRS 3-1-3','4PRS manual','4 PRS optimized','3PRS','2 PRS test','2 PRS ZC-Walsh','ZC-seq-997',...
%                 '2 PRS Test 2','extra 2 PRS test','Barker','Walsh Hadamard','ZC-seq-63','ZC-seq-839','Barker 8PSK'};
            
            
            for ii=1:length(simConfig.flag)
                Ind              = find(Position_Error_All(:,ii)<250);
                fail               = find(Position_Error_All(:,ii)>250);
                failedUEpos = UEPositions(fail,:);

                Fail_Mat{ii}= fail;
                Fail_Pos{ii}= failedUEpos;
                LocRate(ii)           = numel(Ind)/simConfig.numUEs;%/ChanNum;
                Position_Error_Config = Position_Error_All(:,ii);
                
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
                    
                    
                end
            end
        
end
    end
end