function [signalOut,pathGains,sampleTimes,channelOut,averagePathGainsOut,rngStateOut] = nrCDLChannel2(signalIn,channel,rngStateIn)
% nrCDLChannel sends an input signal through a clustered delay line (CDL)
% multi-input multi-output (MIMO) link-level fading channel to obtain the
% channel-impaired signal. The object implements the following aspects of
% TR 38.901 [1]:
%Section 7.7.1: CDL models
%Section 7.7.3: Scaling of delays 
%Section 7.7.5.1: Scaling of angles
%Section 7.7.6: K-factor for LOS channel models
%
% signalIn, a complex vector of size Ns x Nt, where Ns is the number of
% samples, and Nt is the number of transmit antennas. 
% signalOut, a complex vector of size Ns x Nr, where Ns is the number of
% samples, and Nr is the number of receive antennas. 
% pathGains, a complex vector of size Ncs x Np x Nt x Nr, where Ncs is the
% number of channel snapshots, and Np is the number of paths (clusters).
% sampleTimes, sample times of channel snapshots, a real column vector of
% size Ncs x 1.
% rngStateIn, the state Matlab's random number generator used with rand(),
% randi(), and randn() will be restored to.
% channel, input structure plus eventual overrides by the function. For
% example, `channel.pathDelays' can be scaled and translated.
% averagePathGainsOut, a complex vector of size 1 x Np x Nt x Nr.
% rngStateOut, the state of Matlab's random number generator used with rand(),
% randi(), and randn() right after entering the function.
%
% Note that 38l.901 adopts the points of view of the DL. Thus below the BS
% acts as the transmitter, and the UT as the receiver. Of course, the UL
% channel is simple the reciprocal of the DL channel.

Ns=channel.NumTimeSamples;
Nt=channel.TransmitAntennaArray.Size(1)*channel.TransmitAntennaArray.Size(2);





% Handle Matlab's random number generator.
if (nargin > 2)
    rngStateOut = rng(rngStateIn);
else
    rngStateOut = rng;
end

% TODO:
% Figure out if the polarization gain needs to be normalized away.
% Is the implementation of ´rayPhaseCenterOffset' correct?
% `KFactorScaling' should be set to `true' whenever `DelayProfile' is set
% to `CDL-D' or `CDL-E'. Appropriate values of `KFactorScaling' can be
% found in 38.901, Table 7.5-6; see also Section 7.7.6.
% `DelaySpread' is active whenever `DelayProfile' is set to `CDL-A',
% `CDL-B', `CDL-C', `CDL-D' or `CDL-E'. Suitable values of `DelaySpread'
% can be found in TR 38.901 Section 7.7.3 and Tables 7.7.3-1 and 7.7.3-2.
% `ClusterDelaySpread' is active whenever `DelayProfile' is set to `custom'
% and `NumStrongestClusters' to a value greater than zero.
% 
% The 3D distance seems not passed to this function. APPLY EXTERNALLY?
% Check if any of the additional modeling components described in TR 38.901, Sec. 7.6 are needed.
% Need to handle the time aspect in a more efficient way? (Our input
% signals consists of reference signals and is therefore sparse in the
% frequency-time domains.)
% Handle the `Custom' case. Use `DelaySpread', `HasLOSCluster',
% `KFactorFirstCluster', `NumStrongestClusters', `ClusterDelaySpread'.

% For compatibility, define channel properties as in the NR Matlab Toolbox.
txArrayDefault.Size = [2 2 2 1 1]; % Given as [M N P Mg Ng].
txArrayDefault.ElementSpacing = [0.5 0.5 1.0 1.0]; % Given as [dV dH dgV dgH].
txArrayDefault.PolarizationAngles = [45 -45]; % In degrees. Given as [theta, rho]. Apply only when P=2.
txArrayDefault.Orientation = [0; 0; 0]; % In degrees. Given as the bearing angles [alpha; beta; gamma], which correspond to bearing, downtilt, and slant, respectively; see TR 38.901, Sec. 7.1.
txArrayDefault.Element = '38.901'; % 'isotropic'|'38.901'|'38.855 model 1' in TR 38.901 7.3 and TR 38.855 6.1. Or callback function.
txArrayDefault.PolarizationModel = 'Model-2'; % 'Model-1'|'Model-2' in TR 38.901 7.3.2.

rxArrayDefault.Size = [1 1 2 1 1]; % Given as [M N P Mg Ng].
rxArrayDefault.ElementSpacing = [0.5 0.5 0.5 0.5]; % Given as [dV dH dgV dgH].
rxArrayDefault.PolarizationAngles = [0 90]; % In degrees. Given as [theta, rho]. Apply only when P=2.
rxArrayDefault.Orientation = [0; 0; 0]; % In degrees. Given as the bearing angles [alpha; beta; gamma], which correspond to bearing, downtilt, and slant, respectively; see TR 38.901, Sec. 7.1.
rxArrayDefault.Element = 'isotropic'; % 'isotropic'|'38.901'|'38.855 model 1' in TR 38.901 7.3 and TR 38.855 6.1. Or callback function.
rxArrayDefault.PolarizationModel = 'Model-2'; % 'Model-1'|'Model-2' in TR 38.901 7.3.2.

channelDefault.DelayProfile = 'CDL-A'; % One of 'CDL-A', 'CDL-B', 'CDL-C', 'CDL-D', 'CDL-E', 'Custom'.
channelDefault.PathDelays = 0.0; % In seconds. Only used if 'DelayProfile' is set to 'Custom'.
channelDefault.AveragePathGains = 0.0; % In dB. Only used if 'DelayProfile' is set to 'Custom'.
channelDefault.AnglesAoA = 0.0; % In degrees. Only used if 'DelayProfile' is set to 'Custom'.
channelDefault.AnglesAoD = 0.0; % In degrees. Only used if 'DelayProfile' is set to 'Custom'.
channelDefault.AnglesZoA = 0.0; % In degrees. Only used if 'DelayProfile' is set to 'Custom'.
channelDefault.AnglesZoD = 0.0; % In degrees. Only used if 'DelayProfile' is set to 'Custom'.
channelDefault.HasLOSCluster = false; % true|false. Only used if 'DelayProfile' is set to 'Custom'.
channelDefault.KFactorFirstCluster = 13.3; % In dB. Only used if 'DelayProfile' is set to 'Custom', and 'HasLOSCluster' is set to 'true'.
channelDefault.AngleScaling = false; % true|false. Only used if 'DelayProfile' is set to 'CDL-A', 'CDL-B', 'CDL-C', 'CDL-D', or 'CDL-E'.
channelDefault.AngleSpreads = [5.0 11.0 3.0 3.0]; % Use [CASD, CASA, CZSD, CZSA] as in TR 38.901, Sec. 7.7.1., Step 1 when 'DelayProfile' is set to 'Custom'. Use [ASD ASA ZSD ZSA] as in TR 38.901, Sec. 7.7.5.1 when 'AngleScaling' is set to 'true'.
channelDefault.MeanAngles = [0.0 0.0 0.0 0.0]; % In degrees. Only used if 'AngleScaling' is set to 'true'.
channelDefault.XPR = 10.0; % In dB. Only used if 'DelayProfile' is set to 'Custom'.
channelDefault.DelaySpread = 30e-9; % In seconds. Only used if 'DelayProfile' is set to 'CDL-A', 'CDL-B', 'CDL-C', 'CDL-D', or 'CDL-E'.
channelDefault.CarrierFrequency = 4e9; % In Hz. Should be between 0.5 and 100 GHz.
channelDefault.MaximumDopplerShift = 5; % In Hz.
channelDefault.UTDirectionOfTravel = [0; 90]; % In degrees.
channelDefault.KFactorScaling = false; % true|false. Only used if 'DelayProfile' is set to 'CDL-D' or 'CDL-E'.
channelDefault.KFactor = 9.0; % In dB. Only used if 'KFactorScaling' is set to 'true'.
channelDefault.SampleRate = 30.72e6; % In Hz.
channelDefault.TransmitAntennaArray = txArrayDefault;
channelDefault.ReceiveAntennaArray = rxArrayDefault;
channelDefault.SampleDensity = 64; % Samples per half wavelength. Related to the coefficient generation sampling rate by Fcg = MaximumDopplerShift x 2 x SampleDensity, if SampleDensity ~= Inf, Fcg = SampleRate, otherwise.
channelDefault.NormalizePathGains = true; % true|false. Power of paths gains, averaged over time, is 0 dB.
channelDefault.InitialTime = 0.0; % In seconds.
channelDefault.NumStrongestClusters = 0; % Number of strongest clusters to split into subclusters.
channelDefault.ClusterDelaySpread = 3.90625e-9; % In seconds. Only used if 'DelayProfile' is set to 'Custom' and 'NumStrongestClusters' is greater than zero.
channelDefault.RandomStream = 'mt19937ar with seed'; % Source of random stream.
channelDefault.Seed = 73; % Initial seed.
channelDefault.ChannelFiltering = true; % true|false.
channelDefault.NumTimeSamples = 30720; % Sets the duration of the fading process realization.
channelDefault.NormalizeChannelOutputs = false; % true|false. Normalize channel outputs by the number of receive antennas. Only used if 'ChannelFiltering' is set to 'true'.
channelDefault.R3Distance = 0; % 3D distance between the TRP and the UT. Only used for LOS cases. OBS
channelDefault.FirstPathDelay = -Inf; % In seconds. OBS
channelDefault.FirstPathAngleAoA = NaN; % In degrees. OBS
channelDefault.FirstPathAngleAoD = NaN; % In degrees. OBS
channelDefault.FirstPathAngleZoA = NaN; % In degrees. OBS
channelDefault.FirstPathAngleZoD = NaN; % In degrees. OBS

if (nargin<2)
    channel = channelDefault;
end

% Handle predefined delay profiles.
switch lower(channel.DelayProfile)
    case 'cdl-a'
        channel.PathDelays = [0.0000 0.3819 0.4025 0.5868 0.4610 0.5375 0.6708 0.5750 0.7618 1.5375 1.8978 2.2242 2.1718 2.4942 2.5119 3.0582 4.0810 4.4579 4.5695 4.7966 5.0066 5.3043 9.6586];
        channel.AveragePathGains = [-13.4 0 -2.2 -4 -6 -8.2 -9.9 -10.5 -7.5 -15.9 -6.6 -16.7 -12.4 -15.2 -10.8 -11.3 -12.7 -16.2 -18.3 -18.9 -16.6 -19.9 -29.7];
        channel.AnglesAoD = [-178.1 -4.2 -4.2 -4.2 90.2 90.2 90.2 121.5 -81.7 158.4 -83 134.8 -153 -172 -129.9 -136 165.4 148.4 132.7 -118.6 -154.1 126.5 -56.2];
        channel.AnglesAoA = [51.3 -152.7 -152.7 -152.7 76.6 76.6 76.6 -1.8 -41.9 94.2 51.9 -115.9 26.6 76.6 -7 -23 -47.2 110.4 144.5 155.3 102 -151.8 55.2];
        channel.AnglesZoD = [50.2 93.2 93.2 93.2 122 122 122 150.2 55.2 26.4 126.4 171.6 151.4 157.2 47.2 40.4 43.3 161.8 10.8 16.7 171.7 22.7 144.9];
        channel.AnglesZoA = [125.4 91.3 91.3 91.3 94 94 94 47.1 56 30.1 58.8 26 49.2 143.1 117.4 122.7 123.2 32.6 27.2 15.2 146 150.7 156.1];
        channel.AngleSpreadsIntraCluster = [5 11 3 3]; % [CASD, CASA, CZSD, CZSA] as in TR 38.901, Sec. 7.7.1., Step 1.
        channel.XPR = 10;
        channel.HasLOSCluster = false;
        %channel.DelaySpread = ;
    case 'cdl-b'
        channel.PathDelays = [0.0000 0.1072 0.2155 0.2095 0.2870 0.2986 0.3752 0.5055 0.3681 0.3697 0.5700 0.5283 1.1021 1.2756 1.5474 1.7842 2.0169 2.8294 3.0219 3.6187 4.1067 4.2790 4.7834];
        channel.AveragePathGains = [0 -2.2 -4 -3.2 -9.8 -1.2 -3.4 -5.2 -7.6 -3 -8.9 -9 -4.8 -5.7 -7.5 -1.9 -7.6 -12.2 -9.8 -11.4 -14.9 -9.2 -11.3];
        channel.AnglesAoD = [9.3 9.3 9.3 -34.1 -65.4 -11.4 -11.4 -11.4 -67.2 52.5 -72 74.3 -52.2 -50.5 61.4 30.6 -72.5 -90.6 -77.6 -82.6 -103.6 75.6 -77.6];
        channel.AnglesAoA = [-173.3 -173.3 -173.3 125.5 -88.0 155.1 155.1 155.1 -89.8 132.1 -83.6 95.3 103.7 -87.8 -92.5 -139.1 -90.6 58.6 -79.0 65.8 52.7 88.7 -60.4];
        channel.AnglesZoD = [105.8 105.8 105.8 115.3 119.3 103.2 103.2 103.2 118.2 102.0 100.4 98.3 103.4 102.5 101.4 103.0 100.0 115.2 100.5 119.6 118.7 117.8 115.7];
        channel.AnglesZoA = [78.9 78.9 78.9 63.3 59.9 67.5 67.5 67.5 82.6 66.3 61.6 58.0 78.2 82.0 62.4 78.0 60.9 82.9 60.8 57.3 59.9 60.1 62.3];
        channel.AngleSpreadsIntraCluster = [10 22 3 7]; % [CASD, CASA, CZSD, CZSA] as in TR 38.901, Sec. 7.7.1., Step 1.
        channel.XPR = 8;
        channel.HasLOSCluster = false;
        %channel.DelaySpread = ;
    case 'cdl-c'
        channel.PathDelays = [0 0.2099 0.2219 0.2329 0.2176 0.6366 0.6448 0.6560 0.6584 0.7935 0.8213 0.9336 1.2285 1.3083 2.1704 2.7105 4.2589 4.6003 5.4902 5.6077 6.3065 6.6374 7.0427 8.6523];
        channel.AveragePathGains = [-4.4 -1.2 -3.5 -5.2 -2.5 0 -2.2 -3.9 -7.4 -7.1 -10.7 -11.1 -5.1 -6.8 -8.7 -13.2 -13.9 -13.9 -15.8 -17.1 -16 -15.7 -21.6 -22.8];
        channel.AnglesAoD = [-46.6 -22.8 -22.8 -22.8 -40.7 0.3 0.3 0.3 73.1 -64.5 80.2 -97.1 -55.3 -64.3 -78.5 102.7 99.2 88.8 -101.9 92.2 93.3 106.6 119.5 -123.8];
        channel.AnglesAoA = [-101 120 120 120 -127.5 170.4 170.4 170.4 55.4 66.5 -48.1 46.9 68.1 -68.7 81.5 30.7 -16.4 3.8 -13.7 9.7 5.6 0.7 -21.9 33.6];
        channel.AnglesZoD = [97.2 98.6 98.6 98.6 100.6 99.2 99.2 99.2 105.2 95.3 106.1 93.5 103.7 104.2 93.0 104.2 94.9 93.1 92.2 106.7 93.0 92.9 105.2 107.8];
        channel.AnglesZoA = [87.6 72.1 72.1 72.1 70.1 75.3 75.3 75.3 67.4 63.8 71.4 60.5 90.6 60.1 61.0 100.7 62.3 66.7 52.9 61.8 51.9 61.7 58 57];
        channel.AngleSpreadsIntraCluster = [2 15 3 7]; % [CASD, CASA, CZSD, CZSA] as in TR 38.901, Sec. 7.7.1., Step 1.
        channel.XPR = 7;
        channel.HasLOSCluster = false;
        %channel.DelaySpread = ;
    case 'cdl-d'
        % First row is specular path, part of cluster 1; then cluster 1, ...
        channel.PathDelays = [0 0 0.035 0.612 1.363 1.405 1.804 2.596 1.775 4.042 7.937 9.424 9.708 12.525];
        channel.AveragePathGains = [-0.2 -13.5 -18.8 -21 -22.8 -17.9 -20.1 -21.9 -22.9 -27.8 -23.6 -24.8 -30.0 -27.7];
        channel.AnglesAoD = [0 0 89.2 89.2 89.2 13 13 13 34.6 -64.5 -32.9 52.6 -132.1 77.2];
        channel.AnglesAoA = [-180 -180 89.2 89.2 89.2 163 163 163 -137 74.5 127.7 -119.6 -9.1 -83.8];
        channel.AnglesZoD = [98.5 98.5 85.5 85.5 85.5 97.5 97.5 97.5 98.5 88.4 91.3 103.8 80.3 86.5];
        channel.AnglesZoA = [81.5 81.5 86.9 86.9 86.9 79.4 79.4 79.4 78.2 73.6 78.3 87 70.6 72.9];
        channel.AngleSpreadsIntraCluster = [5 8 3 3]; % [CASD, CASA, CZSD, CZSA] as in TR 38.901, Sec. 7.7.1., Step 1.
        channel.XPR = 11;
        channel.HasLOSCluster = true;
        %channel.DelaySpread = ;
        %channel.KFactorScaling = ;
    case 'cdl-e'
        % First row is specular path, part of cluster 1; then cluster 1, ...
        channel.PathDelays = [0.000 0.000 0.5133  0.5440 0.5630 0.5440 0.7112 1.9092 1.9293 1.9589 2.6426 3.7136 5.4524 12.0034 20.6419];
        channel.AveragePathGains = [-0.03 -22.03 -15.8 -18.1 -19.8 -22.9 -22.4 -18.6 -20.8 -22.6 -22.3 -25.6 -20.2 -29.8 -29.2];
        channel.AnglesAoD = [0 0 57.5 57.5 57.5 -20.1 16.2 9.3 9.3 9.3 19 32.7 0.5 55.9 57.6];
        channel.AnglesAoA = [-180 -180 18.2 18.2 18.2 101.8 112.9 -155.5 -155.5 -155.5 -143.3 -94.7 147 -36.2 -26];
        channel.AnglesZoD = [99.6 99.6 104.2 104.2 104.2 99.4 100.8 98.8 98.8 98.8 100.8 96.4 98.9 95.6 104.6];
        channel.AnglesZoA = [80.4 80.4 80.4 80.4 80.4 80.8 86.3 82.7 82.7 82.7 82.9 88 81 88.6 78.3];
        channel.AngleSpreadsIntraCluster = [5 11 3 7]; % [CASD, CASA, CZSD, CZSA] as in TR 38.901, Sec. 7.7.1., Step 1.
        channel.XPR = 8;
        channel.HasLOSCluster = true;
        %channel.DelaySpread = ;
        %channel.KFactorScaling = ;
    case 'los'
        % First row is specular path, part of cluster 1; then cluster 1, ...
        channel.PathDelays = 0.000;
        channel.AveragePathGains = 0.00;
        channel.AnglesAoD = 0;
        channel.AnglesAoA = -180;
        channel.AnglesZoD = 99.6;
        channel.AnglesZoA = 80.4;
        channel.AngleSpreadsIntraCluster = [0 0 0 0]; % [CASD, CASA, CZSD, CZSA] as in TR 38.901, Sec. 7.7.1., Step 1.
        channel.XPR = 8;
        channel.HasLOSCluster = true;
        %channel.DelaySpread = ;
        %channel.KFactorScaling = ;
    case 'custom'
        channel.AngleSpreadsIntraCluster = channel.AngleSpreads; % [CASD, CASA, CZSD, CZSA] as in TR 38.901, Sec. 7.7.1., Step 1.
    otherwise
        disp('Unknown DelayProfile.');
end

% Sanity check.
if ~(channel.CarrierFrequency >= .5e9 &&  channel.CarrierFrequency <= 100e9)
    disp('CarrierFrequency is outside the validity range.');
end

% Scaling of the Ricean K-factor for LOS channel models according to TR
% 38.901, Sec. 7.7.6. This needs to be done before delay scaling according
% to TR 38.901, Sec. 7.3 (see below).
switch lower(channel.DelayProfile)
    case {'cdl-d','cdl-e'}
        if (channel.KFactorScaling)
            KModel = channel.AveragePathGains(1) - pow2db(sum(db2pow(channel.AveragePathGains(2:end))));
            channel.AveragePathGains = channel.AveragePathGains - channel.KFactor + KModel;
            channel.PathDelays = channel.PathDelays/rmsDelay(channel);
        end
end

% Delay scaling according to TR 38.901, Sec. 7.3.
% When computing the rms using the function rmsDelay(channel) defined below
% with the `PathDelays' and `AveragePathGains' given in Tr 38.901, Tables
% 7.7.1-1--7.7.1-5, a normalized value of 1 is obtained. 
% The scaling is accomplished by multiplying `PathDelays' by `DelaySpread'.
switch lower(channel.DelayProfile)
    case {'cdl-a','cdl-b','cdl-c','cdl-d','cdl-e'}
        channel.PathDelays = channel.DelaySpread * channel.PathDelays;
end

% Delay translation. This feature is necessary to relate the channel delays
% to the actual geometry of the simulation.
if (channel.FirstPathDelay >= 0)
    channel.PathDelays = channel.PathDelays - channel.PathDelays(1) + channel.FirstPathDelay;
end
% Take care of LOS angles, if specified by caller.
if (~isnan(channel.FirstPathAngleAoA))
    channel.AnglesAoA(1) = channel.FirstPathAngleAoA;
end
if (~isnan(channel.FirstPathAngleAoD))
    channel.AnglesAoD(1) = channel.FirstPathAngleAoD;
end
if (~isnan(channel.FirstPathAngleZoA))
    channel.AnglesZoA(1) = channel.FirstPathAngleZoA;
end
if (~isnan(channel.FirstPathAngleZoD))
    channel.AnglesZoD(1) = channel.FirstPathAngleZoD;
end

% Scaling of angles according to TR 38.901, Sec. 7.7.5.1.
switch lower(channel.DelayProfile)
    case {'cdl-a','cdl-b','cdl-c','cdl-d','cdl-e'}
        if (channel.AngleScaling)
            for i = 1:10 % Obs
                [muAoD,muAoA,muZoD,muZoA] = meanAngle(channel);
                [asAoD,asAoA,asZoD,asZoA] = rmsAngle(channel);
                channel.AnglesAoD = WrapToPlusMinus180(channelDefault.AngleSpreads(1)/asAoD * (channel.AnglesAoD - muAoD) + channelDefault.MeanAngles(1));
                channel.AnglesAoA = WrapToPlusMinus180(channelDefault.AngleSpreads(2)/asAoA * (channel.AnglesAoA - muAoA) + channelDefault.MeanAngles(2));
                channel.AnglesZoD = ClipBetween0And180(channelDefault.AngleSpreads(3)/asZoD * (channel.AnglesZoD - muZoD) + channelDefault.MeanAngles(3));
                channel.AnglesZoA = ClipBetween0And180(channelDefault.AngleSpreads(4)/asZoA * (channel.AnglesZoA - muZoA) + channelDefault.MeanAngles(4));
            end
        end
end

% Scaling of the path gains. This is an NR Matlab toolbox specific feature.
if (channel.NormalizePathGains)
    channel.AveragePathGains = channel.AveragePathGains - pow2db(sum(db2pow(channel.AveragePathGains)));
end
% Activate the code below (and comment the code above) to normalize away
% the polarimetric channel matrix.
% if (channel.NormalizePathGains)
%     if (channel.HasLOSCluster)
%         kappa = db2pow(channel.XPR);
%         totalPower = 2*(db2pow(channel.AveragePathGains(1)) + sum(db2pow(channel.AveragePathGains(2:end)))*(1+1/kappa));
%     else
%         kappa = db2pow(channel.XPR);
%         totalPower = 2*sum(db2pow(channel.AveragePathGains))*(1+1/kappa);
%     end
%     channel.AveragePathGains = channel.AveragePathGains - pow2db(totalPower);
% end

% Coefficient computation and channel filtering. Are we suppose to filter the input signal?
if (channel.ChannelFiltering)
    if (size(signalIn,1) ~= Ns)
        disp('WARNING: Size of dimension 1 of signalIn does not match channel.NumTimeSamples.');
    end
    if (size(signalIn,2) ~= Nt)
        disp('WARNING: Size of dimension 1 of signalIn does not match number of transmit panels.');
    end
    
    % If so, first compute dimensionality of input/output signals...
    Ns = channel.NumTimeSamples;
    [Nr,Nt,Np] = size(nrCDLChannelCoeff(channel,0));    
    % ...allocate memory space for the output of the filter.
    signalOut = zeros(Ns,Nr);
    % We also need to allocate output signals related to channel coefficients.
    % Units Fcg = [waves/s] x [samples/wave] = [samples/s].
    Fcg = 2 * channel.MaximumDopplerShift * channel.SampleDensity;
    if (Fcg > channel.SampleRate)
        Fcg = channel.SampleRate; % Covers channel.SampleDensity == Inf.
    end
    Ncs = ceil(Ns*Fcg/channel.SampleRate); % # channel snapshots.
    pathGains = zeros(Ncs,Np,Nt,Nr);
    averagePathGainsOut = zeros(1,Np,Nt,Nr);
    sampleTimes = zeros(Ncs,1);
    channelSnapshotIndex = 1;
    % Finally, apply filtering in the time doamin.
    PathDelaysInSamples = floor(channel.PathDelays*channel.SampleRate);
    rays = [];
    for timeSampleIndex = 1:Ns
        sampleTime = channel.InitialTime + (timeSampleIndex-1)/channel.SampleRate;
        if (mod(timeSampleIndex-1,floor(channel.SampleRate/Fcg)) == 0)
            if (timeSampleIndex == 1)
                [coeff,rays,averageCoeff] = nrCDLChannelCoeff(channel,sampleTime);
                averagePathGainsOut(1,:,:,:) = permute(averageCoeff,[3,2,1]);
            else
                coeff = nrCDLChannelCoeff(channel,sampleTime,rays);
            end
            pathGains(channelSnapshotIndex,:,:,:) = permute(coeff,[3,2,1]);
            sampleTimes(channelSnapshotIndex) = sampleTime;
            channelSnapshotIndex = channelSnapshotIndex + 1;
        end
        for pathIndex = 1:Np
            currentPathSampleIndex = timeSampleIndex - PathDelaysInSamples(pathIndex);
            if (currentPathSampleIndex > 0)
                signalOut(timeSampleIndex,:) = signalOut(timeSampleIndex,:) + signalIn(currentPathSampleIndex,:)*coeff(:,:,pathIndex).';
            end
        end
    end
    % The channel outputs might need to be normalized by the number of receive
    % antennas.
    if (channel.ChannelFiltering && channel.NormalizeChannelOutputs)
        M = channel.ReceiveAntennaArray.Size(1);
        N = channel.ReceiveAntennaArray.Size(2);
        signalOut = signalOut/sqrt(M)/sqrt(N);
    end
else
%    % Assume signalIn is just a delta.
%    signalIn = zeros(Ns,Nt);
%    signalIn(1,:) = 1;
    [coeff,~,averageCoeff] = nrCDLChannelCoeff(channel,0);
    signalOut = [];
    pathGains(1,:,:,:) = permute(coeff,[3,2,1]);
    sampleTimes(1) = channel.InitialTime;
    averagePathGainsOut(1,:,:,:) = permute(averageCoeff,[3,2,1]);
end

channelOut = channel;

% Handle Matlab's random number generator.
if (nargin > 2)
    rng(rngStateOut);
end

return;

% We divide the implementation is several parts. In one part, we compute
% the CDL as specified in TR 38.901. In another part, we apply spatial
% filtering. In another part, we apply time or frequency filtering.

function [coeff,rays,averageCoeff] = nrCDLChannelCoeff(channel,t,rays)
% Generate CDL coefficients according to the properties of 'channel'
% according to TR 38.901, Sec. 7.7.1.

% Step 1. Generate departure and arrival angles.
% TR 38.901, Tab. 7.5-3, ray offset angles within a cluster, given for rms
% angle spread normalized to 1.
IntraClusterOffsetAngles = [0.0447 -0.0447 0.1413 -0.1413 0.2492 -0.2492 0.3715 -0.3715 0.5129 -0.5129 0.6797 -0.6797 0.8844 -0.8844 1.1481 -1.1481 1.5195 -1.5195 2.1551 -2.1551];
numRays = length(IntraClusterOffsetAngles);

numClusters = length(channel.AnglesAoD);
if (nargin<3)
    rays.AnglesAoD = (repmat(channel.AnglesAoD,numRays,1) + channel.AngleSpreadsIntraCluster(1)*repmat(IntraClusterOffsetAngles(:),1,numClusters)).';
    rays.AnglesAoA = (repmat(channel.AnglesAoA,numRays,1) + channel.AngleSpreadsIntraCluster(2)*repmat(IntraClusterOffsetAngles(:),1,numClusters)).';
    rays.AnglesZoD = (repmat(channel.AnglesZoD,numRays,1) + channel.AngleSpreadsIntraCluster(3)*repmat(IntraClusterOffsetAngles(:),1,numClusters)).';
    rays.AnglesZoA = (repmat(channel.AnglesZoA,numRays,1) + channel.AngleSpreadsIntraCluster(4)*repmat(IntraClusterOffsetAngles(:),1,numClusters)).';
end

% Step 2. Coupling of rays within a cluster for both azimuth and elevation.
if (nargin<3)
    for n = 1:numClusters
        rays.AnglesAoA(n,:) = rays.AnglesAoA(n,randperm(numRays));
        rays.AnglesZoD(n,:) = rays.AnglesZoD(n,randperm(numRays));
        rays.AnglesZoA(n,:) = rays.AnglesZoA(n,randperm(numRays));
    end
end

% Step 3. Generate cross-polarization power ratios.
if (nargin<3)
    rays.kappa = db2pow(channel.XPR);
end

% Step 4. Coefficient generation.
% Follow the same procedure as in Steps 10 and 11 in Subclause 7.5, with
% the exception that all clusters are treated sa ``weaker clusters'', i.e.
% no further sub-clusters in delay should be generated. Additional
% clusters representing delay spread of the stronger clusters are already
% provided in Tables 7.7.1-1--7.7.1-5.
% >> Step 10. Draw initial random phases.
% This is done for each cluster n, ray m, and polarimetric combination.
if (nargin<3)
    rays.Phases = 2*pi*rand([numClusters,numRays,2,2]) - pi;
end
% >> Step 11. Generate channel coefficients for each cluster n and
% receiver and transmitter element pair u, s.
numRx = prod(channel.ReceiveAntennaArray.Size);
numTx = prod(channel.TransmitAntennaArray.Size);
coeff = zeros(numRx,numTx,numClusters);
averageCoeff2 = zeros(numRx,numTx,numClusters);
for u = 1:numRx
    for s = 1:numTx
        for n = channel.HasLOSCluster+1:numClusters % Treat all clusters as weak clusters.
% Parallelization of the computation of the channel coefficients:
% coeff(u,s,n,t) 
% = sum_m [ sqrt(G(n)/M) * a(u,n,m).' * P(n,m) * a(s,m,n) * p(u,m,n) * q(s,m,n) * d(n,m,t) ] 
% = sqrt(G(n)/M) * sum_m { p(u,m,n) * a(u,n,m).' * P(n,m) * a(s,m,n) * q(s,m,n) * d(n,m,t) } 
% = sqrt(G(n)/M) * sum_m { p(u,m,n) * [a(u,n,m,1), a(u,n,m,2)] * [P(n,m,1,1), P(n,m,1,2); P(n,m,2,1), P(n,m,2,2)] * [a(s,m,n,1); a(s,m,n,2)] * q(s,m,n) * d(n,m,t) } 
% = sqrt(G(n)/M) * sum_m { p(u,m,n) * [a(u,n,m,1)P(n,m,1,1) + a(u,n,m,2)P(n,m,2,1), a(u,n,m,1)P(n,m,1,2) + a(u,n,m,2)P(n,m,2,2)] * [a(s,m,n,1); a(s,m,n,2)] * q(s,m,n) * d(n,m,t) } 
% = sqrt(G(n)/M) * sum_m { p(u,m,n) * [a(u,n,m,1)P(n,m,1,1)a(s,m,n,1) + a(u,n,m,2)P(n,m,2,1)a(s,m,n,1) + a(u,n,m,1)P(n,m,1,2)a(s,m,n,2) + a(u,n,m,2)P(n,m,2,2)a(s,m,n,2)] * q(s,m,n) * d(n,m,t) } 
% = sqrt(G(n)/M) * sum_m { p(u,m,n) * [sum_{i,j} a(u,n,m,i)P(n,m,i,j)a(s,m,n,j)] * q(s,m,n) * d(n,m,t) } 
            a_u_n = panelResponse(channel.ReceiveAntennaArray,rays.AnglesAoA(n,:),rays.AnglesZoA(n,:),u);    %ddddddd
            a_u_n_1 = vec(a_u_n(1,:));
            a_u_n_2 = vec(a_u_n(2,:));
            P_n_1_1 = vec(exp(1i*squeeze(rays.Phases(n,:,1,1))))*1;
            P_n_1_2 = vec(exp(1i*squeeze(rays.Phases(n,:,1,2))))*1/sqrt(rays.kappa);
            P_n_2_1 = vec(exp(1i*squeeze(rays.Phases(n,:,2,1))))*1/sqrt(rays.kappa);
            P_n_2_2 = vec(exp(1i*squeeze(rays.Phases(n,:,2,2))))*1;
            a_s_n = panelResponse(channel.TransmitAntennaArray,rays.AnglesAoD(n,:),rays.AnglesZoD(n,:),s);
            a_s_n_1 = vec(a_s_n(1,:));
            a_s_n_2 = vec(a_s_n(2,:));
            coeff_u_s_n = a_u_n_1.*P_n_1_1.*a_s_n_1 + a_u_n_1.*P_n_1_2.*a_s_n_2 + a_u_n_2.*P_n_2_1.*a_s_n_1 + a_u_n_2.*P_n_2_2.*a_s_n_2;
            p_u_n = vec(rayPhaseCenterOffset(channel.ReceiveAntennaArray,u,rays.AnglesAoA(n,:),rays.AnglesZoA(n,:)));
            q_s_n = vec(rayPhaseCenterOffset(channel.TransmitAntennaArray,s,rays.AnglesAoD(n,:),rays.AnglesZoD(n,:)));
            d_n_t = vec(rayDopplerPhaseOffset(channel,rays.AnglesAoA(n,:),rays.AnglesZoA(n,:),t));
            coeff(u,s,n) = sqrt(db2pow(channel.AveragePathGains(n))/numRays) * sum(p_u_n .* coeff_u_s_n .* q_s_n .* d_n_t);
% Let Q = [1 1/kappa; 1/kappa 1], a2(u,n,m) = abs(a(u,n,m)).^2, and% a2(s,n,m) = abs(a(s,n,m)).^2. 
% Then:
% averageCoeff2(u,s,n) 
% = sum_n [ G(m)/M  * a2(u,n,m).' * Q * a2(s,n,m) ]
% = G(m)/M  * sum_n [ a2(u,n,m).' * Q * a2(s,n,m) ]
% = G(m)/M  * sum_n { [a2(u,n,m,1), a2(u,n,m,2)] * [q11, q12; q21, q22] * [a2(s,n,m,1); a2(s,n,m,1)] }
% = G(m)/M  * sum_n { sum_{i,j} a2(u,n,m,i)*q(i,j)*a2(s,n,m,j) }
            averageCoeff2_u_s_n = abs(a_u_n_1).^2.*abs(a_s_n_1).^2 + abs(a_u_n_1).^2.*abs(a_s_n_2).^2/sqrt(rays.kappa) + abs(a_u_n_2).^2.*abs(a_s_n_1).^2/sqrt(rays.kappa) + abs(a_u_n_2).^2.*abs(a_s_n_2).^2;
            averageCoeff2(u,s,n) = (db2pow(channel.AveragePathGains(n))/numRays) * sum(averageCoeff2_u_s_n);                        
        end        
    end
end
if (channel.HasLOSCluster)
    for u = 1:numRx
        for s = 1:numTx
            % In the LOS case, compute the LOS channel coefficient according to
            % (7.5-29). The LOS parameters are those of cluster number 1.
            coeff(u,s,1) =...
                sqrt(db2pow(channel.AveragePathGains(1)))*...
                panelResponse(channel.ReceiveAntennaArray,channel.AnglesAoA(1),channel.AnglesZoA(1),u).'*...
                [1 0;0 -1]*...%[1 1/sqrt(rays.kappa); 1/sqrt(rays.kappa) 1]*...%[1 0;0 -1]*... OBS
                panelResponse(channel.TransmitAntennaArray,channel.AnglesAoD(1),channel.AnglesZoD(1),s)*...
                rayPhaseCenterOffset(channel.ReceiveAntennaArray,u,channel.AnglesAoA(1),channel.AnglesZoA(1))*...
                rayPhaseCenterOffset(channel.TransmitAntennaArray,s,channel.AnglesAoD(1),channel.AnglesZoD(1))*...
                rayDopplerPhaseOffset(channel,channel.AnglesAoA(1),channel.AnglesZoA(1),t)*...
                ray3DPhaseOffset(channel);
%   (a1c11b1 + a1c12b2/sqrt(kappa) + a2c21b1/sqrt(kappa) + a2c22b2)(a1c11b1 + a1c12b2/sqrt(kappa) + a2c21b1/sqrt(kappa) + a2c22b2)* 
% = |a1|^2|b1|^2 + |a1|^2c11c12*b1b2*/sqrt(kappa) + a1a2*c11c21*|b1|^2/sqrt(kappa) + a1a2*c11c22*b1b2* + |a1|^2|b2|^2/kappa + ... + |a2|^2|b1|^2/kappa + ... + |a2|^2|b2|^2
% ~ |a1|^2|b1|^2 + |a1|^2|b2|^2/kappa + |a2|^2|b1|^2/kappa + |a2|^2|b2|^2 (holds as m\to\infty for fixed gain per cluster) 
% = [|a1|^2 |a2|^2][1 1/kappa;1/kappa 1][|b1|^2;|b2|^2]
            averageCoeff2(u,s,1) =...
                db2pow(channel.AveragePathGains(1))*...
                abs(panelResponse(channel.ReceiveAntennaArray,channel.AnglesAoA(1),channel.AnglesZoA(1),u).').^2*...
                ([1 0;0 -1])*...
                abs(panelResponse(channel.TransmitAntennaArray,channel.AnglesAoD(1),channel.AnglesZoD(1),s)).^2;
        end
    end   
end
averageCoeff = sqrt(averageCoeff2);

function y = ray3DPhaseOffset(channel)
% See (7.5-29) in TR38.901.
lambda0 = physconst('lightspeed')/channel.CarrierFrequency;
y = exp(-1i*2*pi/lambda0*channel.R3Distance);

function y = rayPhaseCenterOffset(antennaArray,elementIndex,angleAzimuth,angleZenith)
% Here, elementLoc is defined in TR38.901, Sec. 7.5 as the location vector of
% antenna element elementIndex for panelIndex = 1,...,Mg*Ng*P*M*N.
[iV,iH,~,igV,igH] = ind2sub(antennaArray.Size,elementIndex); % Sizes order is [M,N,P,Mg,Ng]
dV = antennaArray.ElementSpacing(1);
dH = antennaArray.ElementSpacing(2);
dgV = antennaArray.ElementSpacing(3);
dgH = antennaArray.ElementSpacing(4);
%elementLoc = [+(dV*(iV-1)+dgV*(igV-1)) -(dH*(iH-1)+dgH*(igH-1)) 0].'; % Obs! lambda0 is cancelled out.
elementLoc = [dH*(iH-1) + dgH*(igH-1) 0 dV*(iV-1) + dgV*(igV-1)].'; % Obs! lambda0 is cancelled out.
[phiLCS,thetaLCS] = GCS2LCS(angleAzimuth(:),angleZenith(:),antennaArray.Orientation);
rayDir = [sin(deg2rad(thetaLCS)).*cos(deg2rad(phiLCS)),sin(deg2rad(thetaLCS)).*sin(deg2rad(phiLCS)),cos(deg2rad(thetaLCS))].';
y = exp(1i*2*pi*rayDir.'*elementLoc);

function y = rayDopplerPhaseOffset(channel,angleAzimuth,angleZenith,t)
% Compute the Doppler frequency component in TR 38.901, (7.5-22). The
% Doppler shift can be written as 
% nu(n,m) = <rxDir(n,m),ueVelocity(n,m)> / lambda0
%         = c(n,m) * maximumDopplerShift,
% where <.,.> denote the dot product, c(n,m) = <rxDir(n,m),ueDir(n,m)>
% with |rxDir(n,m)| = |ueDir(n,m)| = 1, 
% and maximumDopplerShift = f*|ueVelocity|/speedOfLight.
myRxDir = [sin(deg2rad(angleZenith(:))).*cos(deg2rad(angleAzimuth(:))),sin(deg2rad(angleZenith(:))).*sin(deg2rad(angleAzimuth(:))),cos(deg2rad(angleZenith(:)))].';
myUTDir = [sin(deg2rad(channel.UTDirectionOfTravel(2)))*cos(deg2rad(channel.UTDirectionOfTravel(1))),sin(deg2rad(channel.UTDirectionOfTravel(2)))*sin(deg2rad(channel.UTDirectionOfTravel(1))),cos(deg2rad(channel.UTDirectionOfTravel(2)))].';
y = exp(1i*2*pi*channel.MaximumDopplerShift*myRxDir.'*myUTDir*t);

function mu = meanDelay(channel)
mu = sum(channel.PathDelays .* db2pow(channel.AveragePathGains))/sum(db2pow(channel.AveragePathGains));

function as = rmsDelay(channel)
mu = meanDelay(channel);
as = sqrt(sum((channel.PathDelays - mu).^2 .* db2pow(channel.AveragePathGains))/sum(db2pow(channel.AveragePathGains)));

function [muAoD,muAoA,muZoD,muZoA] = meanAngle(channel)
% Include ray offsets.
% TR 38.901, Tab. 7.5-3, ray offset angles within a cluster, given for rms
% angle spread normalized to 1.
IntraClusterOffsetAngles = [0.0447 -0.0447 0.1413 -0.1413 0.2492 -0.2492 0.3715 -0.3715 0.5129 -0.5129 0.6797 -0.6797 0.8844 -0.8844 1.1481 -1.1481 1.5195 -1.5195 2.1551 -2.1551];
numRays = length(IntraClusterOffsetAngles);
numClusters = length(channel.AnglesAoD);
anglesAoD = repmat(channel.AnglesAoD,numRays,1) + channel.AngleSpreadsIntraCluster(1)*repmat(IntraClusterOffsetAngles(:),1,numClusters);
anglesAoA = repmat(channel.AnglesAoA,numRays,1) + channel.AngleSpreadsIntraCluster(2)*repmat(IntraClusterOffsetAngles(:),1,numClusters);
anglesZoD = repmat(channel.AnglesZoD,numRays,1) + channel.AngleSpreadsIntraCluster(3)*repmat(IntraClusterOffsetAngles(:),1,numClusters);
anglesZoA = repmat(channel.AnglesZoA,numRays,1) + channel.AngleSpreadsIntraCluster(4)*repmat(IntraClusterOffsetAngles(:),1,numClusters);
averagePathGainsLin = db2pow(repmat(channel.AveragePathGains,numRays,1))/numRays;
% Calculation of mean angle accoring to TR 38.901, Annex A.2.
muAoD = WrapToPlusMinus180(rad2deg(angle(sum(sum(exp(1i*deg2rad(anglesAoD)).*averagePathGainsLin)))));
muAoA = WrapToPlusMinus180(rad2deg(angle(sum(sum(exp(1i*deg2rad(anglesAoA)).*averagePathGainsLin)))));
muZoD = WrapToPlusMinus180(rad2deg(angle(sum(sum(exp(1i*deg2rad(anglesZoD)).*averagePathGainsLin)))));
muZoA = WrapToPlusMinus180(rad2deg(angle(sum(sum(exp(1i*deg2rad(anglesZoA)).*averagePathGainsLin)))));

function [asAoD,asAoA,asZoD,asZoA] = rmsAngle(channel)
% Include ray offsets.
% TR 38.901, Tab. 7.5-3, ray offset angles within a cluster, given for rms
% angle spread normalized to 1.
IntraClusterOffsetAngles = [0.0447 -0.0447 0.1413 -0.1413 0.2492 -0.2492 0.3715 -0.3715 0.5129 -0.5129 0.6797 -0.6797 0.8844 -0.8844 1.1481 -1.1481 1.5195 -1.5195 2.1551 -2.1551];
numRays = length(IntraClusterOffsetAngles);
numClusters = length(channel.AnglesAoD);
anglesAoD = repmat(channel.AnglesAoD,numRays,1) + channel.AngleSpreadsIntraCluster(1)*repmat(IntraClusterOffsetAngles(:),1,numClusters);
anglesAoA = repmat(channel.AnglesAoA,numRays,1) + channel.AngleSpreadsIntraCluster(2)*repmat(IntraClusterOffsetAngles(:),1,numClusters);
anglesZoD = repmat(channel.AnglesZoD,numRays,1) + channel.AngleSpreadsIntraCluster(3)*repmat(IntraClusterOffsetAngles(:),1,numClusters);
anglesZoA = repmat(channel.AnglesZoA,numRays,1) + channel.AngleSpreadsIntraCluster(4)*repmat(IntraClusterOffsetAngles(:),1,numClusters);
averagePathGainsLin = db2pow(repmat(channel.AveragePathGains,numRays,1))/numRays;
% Calculation of mean angle accoring to TR 38.901, Annex A.1.
pow = sum(db2pow(channel.AveragePathGains));
asAoD = rad2deg(sqrt(-2*log(abs(sum(sum(exp(1i*deg2rad(anglesAoD)).*averagePathGainsLin))/pow))));
asAoA = rad2deg(sqrt(-2*log(abs(sum(sum(exp(1i*deg2rad(anglesAoA)).*averagePathGainsLin))/pow))));
asZoD = rad2deg(sqrt(-2*log(abs(sum(sum(exp(1i*deg2rad(anglesZoD)).*averagePathGainsLin))/pow))));
asZoA = rad2deg(sqrt(-2*log(abs(sum(sum(exp(1i*deg2rad(anglesZoA)).*averagePathGainsLin))/pow))));





