%% Class for # location server or any estimator calculating the coordinates
% Contains methods to from time data calculate the coordinates of a UE.
% Deployment may be in location server / LMF, in UE modem or in any 3rd
% party application either on UE side or in cloud. However for different
% deployment different informations/parameters are availible.

classdef ls < handle
    
    properties % Public
        ls_conf     % Configuration of the location server, this
        trp_conf    % Configuration of each trp
        geo_max_dim % Maximum geographical dimenstions
        nr          % nr specific configurations
    end
    
    methods
        
        %% Constructor of gnb for the trp_conf configuration
        function this = ls(ls_conf, trp_conf, geo_max_dim, nr)
            
            if nargin == 4
                this.ls_conf = ls_conf;
                this.trp_conf = trp_conf;
                this.geo_max_dim = geo_max_dim;
                this.nr = nr;
            else
                error('Wrong number of arguments for gnb creation.')
            end
            
        end
        
        %% Select the smallest OTOA among three measurement
        function OTOA=OTOASelect(obj, OTOA)
            Num_trp=length(obj.trp_conf);
            if sum(Num_trp == [9,12,15,18])
                switch obj.ls_conf.OTOAselect
                    case 'one'             %Select the smallest one from 3 OTOAs coresponding to 3 sector from one gnb
                        tt=nan(1,6);
                        I=nan(1,6);
                        for ii=1:Num_trp/3
                            ran=OTOA(ii*3-2 : ii*3);
                            [tt(ii),I(ii)]=min(ran);
                            I(ii)=I(ii)+(ii-1)*3;
                        end
                        OTOA=tt;
                        obj.trp_conf=obj.trp_conf(I);
                        
                        
                    case 'two'             %Select the smallest two from 3 OTOAs coresponding to 3 sector from one gnb
                        tt1=nan(1,6);I1=nan(1,6); % the smallest one
                        tt2=nan(1,6);I2=nan(1,6); % the second smallest one
                        for ii=1:Num_trp/3
                            ran=OTOA(ii*3-2 : ii*3);
                            [tt1(ii),I1(ii)]=min(ran);
                            ran(I1(ii))=Inf;
                            [tt2(ii),I2(ii)]=min(ran);
                            I1(ii)=I1(ii)+(ii-1)*3;
                            I2(ii)=I2(ii)+(ii-1)*3;
                        end
                        I=[I1 I2];
                        OTOA=[tt1 tt2];
                        obj.trp_conf=obj.trp_conf(I);
                        
                    case 'all'
                    otherwise
                        warning('ls_conf.OTOAselect >>>> Unexpected input arguement.')
                end
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