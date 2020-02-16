classdef ue < handle
    
    properties
        conf        % Struct with configurations for the UE
        nr
    end
    
    methods
        
        %% Constructor of ue
        function this = ue(conf, nr)
            if nargin == 2
                this.conf = conf;
                this.nr = nr;
            else
                error('Wrong number of input argument for class ue constructor')
            end
        end
        
        %% "Chest" (read channel from channel model) and combining 
        function [out] = comb(obj, ue_rx, H, snr, nb_tx)
            switch obj.conf.combiner 
                case 1  %Zero forcing 
                    W = pinv(H); 
                case 2 % MMSE 
                    N0 = 1./snr;
                    C  = cov(nb_tx);
                    W =C*H'*inv((H*C*H'+N0*eye(size(H*H',1))));   
                case 3 % Matched filter 
                    W = H';
            end
            out  = W*ue_rx;
        end

        %% Decode PRS - measure TOA sample per "best" beam
        function [OTOA,corr_test] = decode_prs(obj, rx_sig, ref_sig)
            N_ofdm = obj.nr.N_ofdm;
            N_s_ofdm = 2048; % 
            N_s_cp_ofdm = [160, 144, 144, 144, 144, 144, 144, 160, 144, 144, 144, 144, 144, 144];%obj.nr.N_s_cp_ofdm;
            prs_symbols = [3,4,5,6,7,8,9,10];
            %The search grid is depend on the range of the secnario.
            Len_corr = obj.conf.searchGrid*2+1;
            searchGrid=obj.conf.searchGrid;    % 200   periodicity of the sequence in time domain
            interp_fac=1;                      % interpolation factor
             
            % sf_t subframe in time indexed per symbol with removed CP
            num_cells = length(ref_sig(1,:));
            sf_t = nan(num_cells, N_ofdm, N_s_ofdm );
            ref_t = nan(num_cells, N_ofdm, N_s_ofdm );
            % Remove CP per symbol for signals from all cells
            ofdm_start = nan(1,N_ofdm);
            for s = 1:N_ofdm
                ofdm_start(s) = (s-1)*N_s_ofdm + 1 + sum(N_s_cp_ofdm(1:s));
            end
            for s = 1:N_ofdm
                for c = 1:num_cells
                    sf_t(c,s,:) = rx_sig(ofdm_start(s) : (ofdm_start(s) + N_s_ofdm - 1), c);
                    ref_t(c,s,:) = ref_sig(ofdm_start(s) : (ofdm_start(s) + N_s_ofdm - 1), c );
                end
            end
            
            % Correlate the reference signal with the UE Rx signal
            OTOA = nan(1, num_cells);
            switch obj.conf.corr_method               
                case 'time'   
                   corr_test = zeros(Len_corr*interp_fac,num_cells);                  %Buffe for the six corr 
                   for c = 1:num_cells
                        % Correlation per symbol per cell
                        corr_s = zeros(Len_corr,length(prs_symbols));  %Buffer for corr of each symbol
                        
                        for s = prs_symbols+1  
                           
                                corr_s(:,s-prs_symbols(1)) = xcorr(squeeze(sf_t(c, s, :)), squeeze(ref_t(c, s, :)),searchGrid);

                        end
                        corr   = mean(corr_s');    %Can be improved 
                        
                        corr_itp=interp(corr,interp_fac);
                        OTOA(c)    = obj.ToAEstimation(corr_itp,interp_fac,0); 
                        
                        %[~,OTOAIntpl(k)]= ToAEstimation(Corr,1); 
                        corr_test(:,c)=abs(corr_itp);                               %dddd   
                   end 

                   
                   
                  % didn't check this part   
                case 'Frequency Rx'                                        %That is not the case that we need to concider now 
                    sf_f = nan(num_cells,num_ofdm,num_fft);
                    ref_f = nan(num_cells,num_ofdm,num_fft);
                    for c = 1:num_cells 
                        corr_s = [];
                        for s = 1:length(N_ofdm)
                            sf_f(c,s,:) = fft(sf_t(c,s,:), N_fft);
                            ref_f(c,s,:) = conj(fft(ref_t(c, s, :), N_fft));
                            corr_s = [corr_s (ref_f(c,s,:).*sf_f(s,:))'];
                            % Correlation_time = ifft(Correlation);
                            % CorrSym(c, :) = Correlation_time;
                        end
                        corr_mean = [corr_mean mean(corr_s)];
                        OTOA(c)            = obj.ToAEstimation(corr_mean(c),0);
                    end  
            end            
        end
        
        %% Estimate TimeOfArrival from a cross correlation ouptput
        function TOA=ToAEstimation(obj,Corr,interp_fac,Ind)
                
            eta = obj.conf.eta;   %0.8
            dr = obj.conf.dr;     %1
            Fs = obj.nr.Fsamp;    % dddd Zhang: not correct

            [maxVal, maxCorrPos] = max(abs(Corr));
            
            PAR        = maxVal/mean(abs(Corr));
            Delay      = [];
            TOA=nan;
            if PAR > eta  % PRS presented
                
                if Ind == 0   
                    
                    % estimate the ToA of the first path
                    aa = 1:1:length(Corr); % ascending vector for indices comparison
                    Normalized_Corr = abs(Corr)./maxVal;
                    bb = find(Normalized_Corr>=eta); % find qualified value (better than threshold)
                    cc = setdiff(aa, bb); % select all values not qualified by threshold
                    Normalized_Corr(cc) = 0;
                    
                    [~,firstCorrPos] = findpeaks(Normalized_Corr);
                    CorrFlag = isempty(firstCorrPos);
                    
                    if CorrFlag == 0
                        TOA = (firstCorrPos(1)-(100*2)*interp_fac/2)/Fs/interp_fac; % dddddddd -1 here bcs in MATLAB all start from 1   Zhang: Fs here is wrong!
                        %Fs;
                        %TOA             %ddddddd
                    else
                        firstCorrPos         = bb;
                        TOA                    = (firstCorrPos(1)-(100*2)*interp_fac/2)/Fs/interp_fac;
                    end
 
                    
                    
                    
                else 
                    
                    % interpolation at found max positions, take care of the indices
                    startPos = max(maxCorrPos-SearchGrid,1); % return 1 in case of negative value
                    endpos = min(maxCorrPos+SearchGrid,length(Corr));
                    tempcorr = Corr(startPos:endpos);
                    CorrIntpl = interp(tempcorr,dr);                       %Zhang: dr=1 and it does not change anything
                    
                    if obj.conf.doMPS == true
                        [TOAOffset, ~] = MPD(CorrIntpl,obj.conf);
                        TOA = ((startPos-1)*dr+TOAOffset-1)/Fs/dr;
                    else
                        aa  = 1:1:length(CorrIntpl); % ascending vecor for indices comparison 
                        Normalized_Corr = abs(CorrIntpl)./maxVal;
                        bb = find(Normalized_Corr>=eta); % find qualified value (better than threshold)
                        cc = setdiff(aa, bb); % select values not qualified by threshold
                        Normalized_Corr(cc) = 0;
                        [~,firstCorrPos] = findpeaks(Normalized_Corr);
                        CorrFlag = isempty(firstCorrPos);
                        
                        if CorrFlag == 0
                            TOA = ((startPos-1)*dr+firstCorrPos(1)-1)/Fs/dr; % -1 here bcs in MATLAB all start from 1
                        else
                            firstCorrPos = bb;
                            TOA = ((startPos-1)*dr+firstCorrPos(1)-1)/Fs/dr;
                        end
                        %             firstCorrPos  = find(abs(CorrIntpl)./max(abs(CorrIntpl))>=eta2);   % First path detected
                        %             TOA           = ((startPos-1)*dr+firstCorrPos(1)-1)/Fs/dr;
                        
                    end
                    
                end              
            end
          
        end
        
        %% Correlation for all beams and transmission points in vectors
        function outp = xcorr_prs(obj, ref, inp)
            ll = size(inp);
            outp = [];
            for i = 1:ll(end)
                outp = [outp xcorr(ref(:,i), inp(:,i))];
            end
        end
        
    end
    
end