 classdef environment < handle
    
    properties (Constant)
        c = physconst('LightSpeed'); % ~ Speed of electromagnetic signals in vacuum - m/s
    end
   
    properties 
        trp_pos     % 3D array with transmission point coordinates {x, y, z}
        trp_bearing % Bearing of the transmission points
        max_dim     % Absolute upper bound of [x, y, z]
        channel     % Settings for nrCDLChannel channel model
        ue_pos      % 3D array with UE coordinates {x, y, z} or 'random'
        trpselect
        num_trp
        nr
        Fsamp
        len_channel
    end
    
    methods 
        
        %% Environment constructor
        function obj = environment(trp_pos, trp_bearing, max_dim, channel, ue_pos , trpselect , num_trp , Fsamp ,len_channel)
            
            if nargin == 9
                obj.trp_pos = trp_pos;
                obj.trp_bearing = trp_bearing;
                obj.max_dim = max_dim;
                obj.channel = channel;
                obj.ue_pos = ue_pos;
                obj.trpselect = trpselect;
                obj.num_trp = num_trp;
                obj.Fsamp = Fsamp;
                obj.len_channel=len_channel;
            else
                error('Wrong number of arguments for environment creation.');
            end
        end
        
        %% set environment 
        function Pos=setenvironment(obj)
            Pos=struct('Num',[],'UTs',[],'trp',[]);
            
            % Drop UE and store the info into Pos.UTs
            Pos.UTs=obj.dropUE;                
             %figRoom = obj.plotRoom(obj.trp,obj.UTs,1.0);
            
            % configure all the 12 gnbs (36 trps)
            Pos.trp = obj.dropTRPs( Pos.UTs.UTDirectionOfTravel, Pos.UTs.Mobility, Pos.UTs.ArrayPosition); 
 
            % Select 6 closest TRPs...
            Pos.trp = obj.selectClosestTRPs(Pos.trp , Pos.UTs);
            
            % ...and reassign CellIDs.
            for TRPIndex = 1:length(Pos.trp)
                Pos.trp(TRPIndex).CellID = mod(TRPIndex-1,6);
            end
           
        end
        
        %% Method to place new UE "2"
        function UTs=dropUE(obj)            
                  
             switch obj.ue_pos 
                 case 'random'
                 UTs = struct('ArrayPosition',[rand(1)*obj.max_dim(1); rand(1)*obj.max_dim(2); 1.5],...
                             'UTDirectionOfTravel',[rand(1)*360 - 180;0],...
                             'Mobility',3000/3600);
                 case 'fixed'
                 UTs = struct('ArrayPosition',[51.4;21.15; 1.5],...
                             'UTDirectionOfTravel',[rand(1)*360 - 180;0],...
                             'Mobility',3000/3600);
                 otherwise
                 error('Non-random UE placement - not implemented')
             end   
             
        end
        
        %% Method to place TRPs in geometry "2"
        function TRPs = dropTRPs(obj, UTDirectionOfTravel, Mobility, UTPosition)
            TRPPositions=obj.trp_pos;
            TRPBearings=obj.trp_bearing;
   
            % Codebook of beams equi-spaced in azimuth when using an array defined by
            % txArrayFR2.Size = [4 8 1 1 1]; % Given as [M N P Mg Ng].
            % txArrayFR2.ElementSpacing = [0.5 0.5 0.0 0.0]; % Given as [dV dH dgV dgH].
            phiGCS = repmat(37.5 + 15*(7:-1:0),8,1);% 8X8
            codebookdft = exp(-1i*pi*repmat((0:7).',1,8).*cos(deg2rad(phiGCS)));  %8X8 double 
            %codebookdft = fftshift(dftmtx(8),2);
            codebook(1,:,:) = codebookdft/sqrt(8)/sqrt(4);
            codebook(2,:,:) = codebookdft/sqrt(8)/sqrt(4);
            codebook(3,:,:) = codebookdft/sqrt(8)/sqrt(4);
            codebook(4,:,:) = codebookdft/sqrt(8)/sqrt(4);
            codebook = reshape(codebook,4*8,8);
            
            % Generate a list of TRPs with the desired positions and bearing angles.
            TRPIndex = 1;
            for positionIndex = 1:size(TRPPositions,2)    % TRPPosition: 3 by 12
                for sectorIndex = 1:length(TRPBearings)
%TransmitAntennaArrayOrientation,DelayProfile,CarrierFrequency,MaximumDopplerShift,UTDirectionOfTravel,SampleRate,NumTimeSamples,FirstPathDelay,TRP2UTDir)
                    channel = getDefaultChannel(obj, ...
                        [TRPBearings(sectorIndex); 0; 0],...
                        obj.channel,...      %            'LOS',... obs
                        30e9,...             % CarrierFrequency
                        30e9*Mobility/physconst('lightspeed'),...
                        UTDirectionOfTravel,...
                        obj.Fsamp,...                      %dddddddd???
                        obj.len_channel,...                                          
                        norm(UTPosition-TRPPositions(:,positionIndex))/physconst('lightspeed'),...   % First Path Delay
                        (UTPosition-TRPPositions(:,positionIndex))/norm(UTPosition-TRPPositions(:,positionIndex)));  % TRP2UTDir 3X1 matrix
                    
                    TRPs(TRPIndex) = struct(...
                        'ArrayPosition',TRPPositions(:,positionIndex),...
                        'ArrayBearing',TRPBearings(sectorIndex),...
                        'ArrayBeamWidth',360/length(TRPBearings),...
                        'ArrayCodebook',codebook,...
                        'Channel',channel,...
                        'CellID',TRPIndex-1);
                    TRPIndex = TRPIndex + 1;
                    
                end
            end
            numTRP = length(TRPs); 
        end
        
        %% Method to get a default channel out from the channel model with some basic settings "3"
        function channel = getDefaultChannel(obj, TransmitAntennaArrayOrientation,DelayProfile,CarrierFrequency,MaximumDopplerShift,UTDirectionOfTravel,SampleRate,NumTimeSamples,FirstPathDelay,TRP2UTDir)
            % The following parameters might need to be changed by the function call:
            % DelaySpread.
            
            % Define the tx array. See Table 6.1.1-3 in TR38.901.
            txArrayFR2.Size = [4 8 1 1 1]; % Given as [M N P Mg Ng].
            %txArrayFR2.Size = [4 8 2 1 1]; % Given as [M N P Mg Ng]. Obs
            txArrayFR2.ElementSpacing = [0.5 0.5 0.0 0.0]; % Given as [dV dH dgV dgH].
            txArrayFR2.PolarizationAngles = [0 90]; % In degrees. Given as [theta, rho]. Apply only when P=2.
            txArrayFR2.Orientation = TransmitAntennaArrayOrientation; % In degrees. Given as the bearing angles [alpha; beta; gamma].
            txArrayFR2.Element = '38.855 model 1'; % 'isotropic'|'38.901'|'38.855 model 1' in TR 38.901 7.3 and TR 38.855 6.1.
            txArrayFR2.PolarizationModel = 'Model-2'; % 'Model-1'|'Model-2' in TR 38.901 7.3.2.
            
            rxArrayFR2.Size = [1 1 1 1 1]; % Given as [M N P Mg Ng].
            %rxArrayFR2.Size = [1 1 2 1 1]; % Given as [M N P Mg Ng]. Obs
            rxArrayFR2.ElementSpacing = [0.5 0.5 0.5 0.5]; % Given as [dV dH dgV dgH].
            rxArrayFR2.PolarizationAngles = [0 90]; % In degrees. Given as [theta, rho]. Apply only when P=2.
            rxArrayFR2.Orientation = [0; 0; 0]; % In degrees. Given as the bearing angles [alpha; beta; gamma].
            rxArrayFR2.Element = 'isotropic'; % 'isotropic'|'38.901' in TR 38.901 7.3.
            rxArrayFR2.PolarizationModel = 'Model-2'; % 'Model-1'|'Model-2' in TR 38.901 7.3.2.
            
            channel.DelayProfile = DelayProfile; % One of 'CDL-A', 'CDL-B', 'CDL-C', 'CDL-D', 'CDL-E', 'Custom'.
            channel.PathDelays = 0.0; % In seconds. Only used if 'DelayProfile' is set to 'Custom'.
            channel.AveragePathGains = 0.0; % In dB. Only used if 'DelayProfile' is set to 'Custom'.
            channel.AnglesAoA = 0.0; % In degrees. Only used if 'DelayProfile' is set to 'Custom'.
            channel.AnglesAoD = 0.0; % In degrees. Only used if 'DelayProfile' is set to 'Custom'.
            channel.AnglesZoA = 0.0; % In degrees. Only used if 'DelayProfile' is set to 'Custom'.
            channel.AnglesZoD = 0.0; % In degrees. Only used if 'DelayProfile' is set to 'Custom'.
            channel.HasLOSCluster = false; % true|false. Only used if 'DelayProfile' is set to 'Custom'.
            channel.KFactorFirstCluster = 13.3; % In dB. Only used if 'DelayProfile' is set to 'Custom', and 'HasLOSCluster' is set to 'true'.
            channel.AngleScaling = false; % true|false. Only used if 'DelayProfile' is set to 'CDL-A', 'CDL-B', 'CDL-C', 'CDL-D', or 'CDL-E'.
            channel.AngleSpreads = [5.0 11.0 3.0 3.0]; % Use [CASD, CASA, CZSD, CZSA] as in TR 38.901, Sec. 7.7.1., Step 1 when 'DelayProfile' is set to 'Custom'. Use [ASD ASA ZSD ZSA] as in TR 38.901, Sec. 7.7.5.1 when 'AngleScaling' is set to 'true'.
            channel.MeanAngles = [0.0 0.0 0.0 0.0]; % In degrees. Only used if 'AngleScaling' is set to 'true'.
            channel.XPR = 10.0; % In dB. Only used if 'DelayProfile' is set to 'Custom'.
            channel.DelaySpread = 30e-9; % In seconds. Only used if 'DelayProfile' is set to 'CDL-A', 'CDL-B', 'CDL-C', 'CDL-D', or 'CDL-E'.
            channel.CarrierFrequency = CarrierFrequency; % In Hz. Should be between 0.5 and 100 GHz.
            channel.MaximumDopplerShift = MaximumDopplerShift; % In Hz.
            channel.UTDirectionOfTravel = UTDirectionOfTravel; % In degrees.
            channel.KFactorScaling = false; % true|false. Only used if 'DelayProfile' is set to 'CDL-D' or 'CDL-E'.
            channel.KFactor = 9.0; % In dB. Only used if 'KFactorScaling' is set to 'true'.
            channel.SampleRate = SampleRate; % In Hz.
            channel.TransmitAntennaArray = txArrayFR2;
            channel.ReceiveAntennaArray = rxArrayFR2;
            channel.SampleDensity = 64; % Samples per half wavelength. Related to the coefficient generation sampling rate by Fcg = MaximumDopplerShift x 2 x SampleDensity, if SampleDensity ~= Inf, Fcg = SampleRate, otherwise.
            channel.NormalizePathGains = true; % true|false. Power of paths gains, averaged over time, is 0 dB.
            channel.InitialTime = 0.0; % In seconds.
            channel.NumStrongestClusters = 0; % Number of strongest clusters to split into subclusters.
            channel.ClusterDelaySpread = 3.90625e-9; % In seconds. Only used if 'DelayProfile' is set to 'Custom' and 'NumStrongestClusters' is greater than zero.
            channel.RandomStream = 'mt19937ar with seed'; % Source of random stream.
            channel.Seed = 73; % Initial seed.
            channel.ChannelFiltering = true; % true|false.
            channel.NumTimeSamples = NumTimeSamples; % Sets the duration of the fading process realization.
            channel.NormalizeChannelOutputs = false; % true|false. Normalize channel outputs by the number of receive antennas. Only used if 'ChannelFiltering' is set to 'true'.
            channel.R3Distance = 0; % 3D distance between the TRP and the UT. Only used for LOS cases. OBS
            channel.FirstPathDelay = FirstPathDelay; % In seconds. OBS
            [angleAoDGCS,angleZoDGCS,~] = cart2sph(TRP2UTDir(1),TRP2UTDir(2),TRP2UTDir(3));
            angleAoDGCS = WrapToPlusMinus180(obj, rad2deg(angleAoDGCS));
            angleZoDGCS = 90 - rad2deg(angleZoDGCS);
            if (angleZoDGCS<0 || angleZoDGCS>180)
                disp('WARNING: thetaLCS out of range after `cart2sph'' conversion.');
            end
            [angleAoAGCS,angleZoAGCS,~] = cart2sph(-TRP2UTDir(1),-TRP2UTDir(2),-TRP2UTDir(3));
            angleAoAGCS = WrapToPlusMinus180(obj,rad2deg(angleAoAGCS));
            angleZoAGCS = 90 - rad2deg(angleZoAGCS);
            if (angleZoAGCS<0 || angleZoAGCS>180)
                disp('WARNING: thetaLCS out of range after `cart2sph'' conversion.');
            end
            channel.FirstPathAngleAoA = angleAoAGCS; % In degrees. OBS
            channel.FirstPathAngleAoD = angleAoDGCS; % In degrees. OBS
            channel.FirstPathAngleZoA = angleZoAGCS; % In degrees. OBS
            channel.FirstPathAngleZoD = angleZoDGCS; % In degrees. OBS
        end 
        
        %% Method to choose the closest by distance to one placed UT "2"
        function TRPs = selectClosestTRPs(obj, TRPs , UTs)
            trpselect=obj.trpselect;
            numTRPs=obj.num_trp;
            switch trpselect
                case 'Auto'
                    % Select `numTRPs' closest to UTs(1).
                    d = zeros(length(TRPs),1);
                    for TRPIndex = 1:length(TRPs)
                        
                        UTBearing = UTs(1).ArrayPosition(1:2) - TRPs(TRPIndex).ArrayPosition(1:2);
                        
                        UTBearing = UTBearing/norm(UTBearing);
                        
                        TRPBearing = [cos(deg2rad(TRPs(TRPIndex).ArrayBearing + 90)),sin(deg2rad(TRPs(TRPIndex).ArrayBearing + 90))].';
                        
                        UT2TRPBearing = rad2deg(acos(UTBearing'*TRPBearing));     %remove WrapToPlusMinus180 since the output of acos is within the range (0,180).
                        
                        if (UT2TRPBearing <= TRPs(TRPIndex).ArrayBeamWidth/2)
                            UTVec = UTs(1).ArrayPosition - TRPs(TRPIndex).ArrayPosition;
                            d(TRPIndex) = norm(UTVec);
                        else
                            UTVec = UTs(1).ArrayPosition - TRPs(TRPIndex).ArrayPosition;
                            d(TRPIndex) = norm(UTVec)+max(obj.max_dim)/2;
                        end
                    end
                    [d,I] = sort(d);
                    if (length(d) >= numTRPs)
                        TRPs = TRPs(I(1:numTRPs));
                    end
                    
                case 'Manual'      % 3 subframes cooresponding to 3 sectors of each gnbs
                    d = zeros(length(TRPs),1);
                    for TRPIndex = 1:length(TRPs)
                        UTVec = UTs(1).ArrayPosition - TRPs(TRPIndex).ArrayPosition;
                        d(TRPIndex) = norm(UTVec);
                    end
                    [~,I1] = sort(d);
                    I=zeros(1,numTRPs);
                    for ii=1:3
                        I((1:6)+(ii-1)*6)=I1((0:3:15)+ii);
                    end
                    if (length(TRPs) >= numTRPs)
                        TRPs = TRPs(I(1:numTRPs));
                    end
                    
                    
            end
        end
        
        %% Method to apply the configured channel "1"   
        function sig_o = apply_ch(obj, sig_i,Trp, UTs)
      
            [signalOut,pathGains,sampleTimes] = nrCDLChannel2(sig_i, Trp.Channel);

            d2D = norm(Trp.ArrayPosition(1:2) - UTs.ArrayPosition(1:2));     %7.4.1
            hBS = Trp.ArrayPosition(3);
            hUT = UTs.ArrayPosition(3);
            [PL,sigmaSF] = nrPathloss('InH',Trp.Channel.HasLOSCluster,d2D,hBS,hUT, Trp.Channel.CarrierFrequency);
            PL = PL + randn(1)*sigmaSF;
            sig_o = signalOut/sqrt(power(10,PL/10));
            
        end
        
        %% WrapToPlusMinus180 "2"
        function y = WrapToPlusMinus180(obj, x)
            % Wrap x to [-180,180).
            y = mod(x + 180,2*180) - 180;
            
        end
        
    end
end