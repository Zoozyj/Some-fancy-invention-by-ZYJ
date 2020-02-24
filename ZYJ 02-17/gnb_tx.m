%% Class for # gnbs or # trps.
% Contains method to generate data of one subframe

classdef gnb_tx < handle
        
    properties % Public 
        trp_conf    % Configuration of each trp
        sig_trp     % Signal out from each trp per beam and sector
    end
    
    properties(Access = private)
        nr % NR specific parameters
    end
    
    methods 
        %% Constructor of gnb for the trp_conf configuration
        function this = gnb_tx(trp_conf, nr)
            
            if nargin == 2
                this.trp_conf = trp_conf;
                this.nr = nr;
            else
                error('Wrong number of arguments for gnb creation.')
            end
            
        end
        
        %% Method to generate one subframe of data
        function sig = gen_sf(obj, slotNumber)
            
            
            % Loop over TRPs
            for TRPIndex = 1:length(obj.trp_conf)
                
                % Generate the PRS to be transmitted. Symbols between the fourth and
                % the eleventh are allocated for PRS transmission.  index: 3-10
                [signalIn,signalInSymbols] = nrPRS(obj.trp_conf(TRPIndex).CellID,...
                    obj.nr.MaxDLNumberRB,...
                    obj.nr.DLNumberRB,...
                    obj.nr.PRSNumberRB,...
                    slotNumber);
                signalInTRPs(:,TRPIndex) = signalIn;
                
                % Apply a different precoding codeword per allocated PRS symbol.
                Ns = length(signalIn);
                Nt = length(obj.trp_conf(TRPIndex).ArrayCodebook(:,1));
                signalInPrecoded = zeros(Ns,Nt);
                for symbolIndex = 3:10
                    symbolStart = signalInSymbols(symbolIndex + 1);        %the index in the matrix is (symbolIndex+1), since symbolIndex are 0,1,2,3,4.....
                    if (symbolIndex + 1 < length(signalInSymbols))
                        symbolEnd = signalInSymbols(symbolIndex + 2) - 1;
                    else
                        symbolEnd = length(signalIn) - 1;                  % Zhang: In which situation that this case would happen? length(SignalInSymbols)>=0?
                    end
                    codewordIndex = mod(symbolIndex - 3,8);
                    
                    signalInPrecoded(symbolStart + 1:symbolEnd + 1,:) = signalIn(symbolStart + 1:symbolEnd + 1)...
                        * obj.trp_conf(TRPIndex).ArrayCodebook(:,codewordIndex + 1).';
                end
                sig(:,:,TRPIndex) = signalInPrecoded;
            end
        end
        
    end
end