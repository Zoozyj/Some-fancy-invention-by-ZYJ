%% Class for # gnbs or # trps.
% Contains method to generate data of one subframe

classdef ls < handle
    
    properties % Public
        ls_conf  % Configuration of the location server, this
        trp_conf % Configuration of each trp
        nr
    end
    
    methods
        
        %% Constructor of gnb for the trp_conf configuration
        function this = ls(trp_conf, nr)
            
            if nargin == 3
                this.ls_conf = ls_conf;
                this.trp_conf = trp_conf;
                this.nr = nr;
            else
                error('Wrong number of arguments for gnb creation.')
            end
            
        end
        
        function coord = localize(obj, OTOA)
            
            switch obj.ls_conf.method
                
                case 'one'
                    estUEPos = estToA_v2(OTOA);
                    
                case 'two'
                    error('not implemented');
            
            end
            
        end
        
        function estUEPos = estToA_v2(obj, OTOA)
            c = obj.nr.c; % Speed of light in m/s
            
            % Find the eNodeBs corresponding to the OTOA measurements
            NodeBPosition= {};
            OTOADetected    = [];
            estToANum       = 0;
            %  a = find (OTOA>=8.8191e-06);
            a = find (OTOA>=1.2018e-05);
            OTOA(a) = NaN;
            
            for i = 1:length(obj.trp_conf.eNodeBs)
                if ~isnan(OTOA(i))
                    estToANum=estToANum+1;
                    NodeBPosition{estToANum} = obj.trp_conf.ArrayPosition{i};
                    OTOADetected(estToANum)    = OTOA(i);
                end
            end
            estUEPos        = nan(1,2);
            
            %OTDOA estiamtion when detected at least three e-NodeBs
            if estToANum>2
                StartStUEPos      = zeros(1,2);
                for k = 1:estToANum
                    StartStUEPos  = StartStUEPos+NodeBPosition{k};
                end
                StartStUEPos      = StartStUEPos./estToANum;
                sumToAerr         = @(UEPos)ToACostFunc(UEPos,estToANum, NodeBPosition,OTOADetected, c);
                [estUEPos, error] = minimize(sumToAerr, zeros(1,2), [], [], [],[], [-3000,-3000], [3000,3000]); % -3000 & 3000 bcs our map axis within these numbers
            end
            
        end
        
        % Estimated UE positions that minimizes the sum OTDOA error
        function sumToAerr=ToACostFunc(obj, UEPos, estToANum, eNodeBsDetected,OTOADetected, c)
            
            for k = 1:estToANum
                positions(k,:)    = NodeBPosition{k};
                tempOTOA(k)       = OTOADetected(k);
            end
            
            xt = UEPos(1);
            yt = UEPos(2);
            sumToAerr = 0;
            for i = 1:estToANum
                for j = i+1:estToANum
                    sumToAerr=sumToAerr+abs(sqrt((positions(j,1)-xt)^2+(positions(j,2)-yt)^2)...
                        -sqrt((positions(i,1)-xt)^2+(positions(i,2)-yt)^2)-c*(tempOTOA(j)-tempOTOA(i)));
                    
                    % sumToAerr=sumToAerr+(abs(sqrt((positions(j,1)-xt)^2+(positions(j,2)-yt)^2)...
                    %     -sqrt((positions(i,1)-xt)^2+(positions(i,2)-yt)^2)-c*(tempOTOA(j)-tempOTOA(i))))/c;
                end
            end
            
        end
        
        
    end
end