%% Class for # location server or any estimator calculating the coordinates
% Contains methods to from time data calculate the coordinates of a UE.
% Deployment may be in location server / LMF, in UE modem or in any 3rd
% party application either on UE side or in cloud. However for different
% deployment different informations/parameters are availible.

classdef loc_est < handle
    
    properties % Public
        ls_conf     % Configuration of the location server, this
        trp_conf    % Configuration of each trp
        geo_max_dim % Maximum geographical dimenstions
        nr          % nr specific configurations
    end
    
    methods
        
        %% Constructor of gnb for the trp_conf configuration
        function this = loc_est(ls_conf, trp_conf, geo_max_dim, nr)
            
            if nargin == 4
                this.ls_conf = ls_conf;
                this.trp_conf = trp_conf;
                this.geo_max_dim = geo_max_dim;
                this.nr = nr;
            else
                error('Wrong number of arguments for gnb creation.')
            end
            
        end
        
        %% Selection of method for localization
        function coord = localize(obj, OTOA)
            
            switch obj.ls_conf.method
                
                case 'one'
                    coord = obj.estToA_v2(OTOA);
                case 'two'
                    error('not implemented');
            
            end
            
        end
        
        %% The estToA_v2 is copy from LTE pos simulator at earlier master thesis work
        function estUEPos = estToA_v2(obj, OTOA)
            c = physconst('LightSpeed'); % Speed of light in m/s
            ss = size(obj.trp_conf);
            N_trp=ss(2);
            
            trp_pos=zeros(3,N_trp); %extract gnb postion form the object to a matrix
            
            for i = 1:N_trp
                    trp_pos(:,i) = obj.trp_conf(i).ArrayPosition;
            end
             
            StartStUEPos    = mean(trp_pos,2) ;                %mean
   
            sum_toa_err = @(UEPos) obj.toa_cost_fnc(UEPos, N_trp , trp_pos', OTOA , c);
                
            [estUEPos, error] = minimize(sum_toa_err, StartStUEPos, [], [], [],[], [0;0;0], obj.geo_max_dim);      

            
        end
 
        
        % Sum errors used to minimizes the OTDOA error
        function sum_toa_err = toa_cost_fnc(obj, ue_pos, toa_nbr, trp_pos, toa, c)
            
            xt = ue_pos(1);
            yt = ue_pos(2);
            zt = ue_pos(3);
            
            sum_toa_err = 0;
            for i = 1:toa_nbr
                for j = i+1:toa_nbr
                    
                    sum_toa_err = sum_toa_err + abs(...
                        sqrt( (trp_pos(j,1)-xt)^2 + (trp_pos(j,2)-yt)^2 + (trp_pos(j,3)-zt)^2)...
                       -sqrt( (trp_pos(i,1)-xt)^2 + (trp_pos(i,2)-yt)^2 + (trp_pos(i,3)-zt)^2)...
                       -c*(toa(j)-toa(i)));
                    
                    
                end
            end
            
        end
        
        
        
    end
end