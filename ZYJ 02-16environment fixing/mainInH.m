function positioningTestInH6TRP 
% positioningTestInH6TRP
% Implements a simplified version of scenario 1, indoor office for FR1 and
% FR2 (open office), described in TR38.855.
% The test scenario is simplified in that only 6 TRPs are considered, and
% certain assumptions are made on the beamforming capabilities of the TRPs.

% TODO: 
% Define some top level parameters to steer the simulations. Ideally, they
% come directly from Tables 6.1.1-1 to 6.1.1-6 in TR 38.855.
MaxDLNumberRB = 110;
DLNumberRB = 110;
PRSNumberRB = 69;
slotNumber = randi(20,1)-1;
d2Dmin = 0; % In meters. Min. gNB-UE distance (2D).
% TODO:
% Try multiple polarizations.
% Cyclic prefix lengths in nrPRS need to be adjusted.
% CHECK POWERS!
% OK: (Not needed) Plot radiation pattern of codebook.
% OK: CEHCK BEAMFORMING!!! 
% Apply transmit power, antenna gains, and receive noise.
% OK: TRY WITH PURE LOS!
% OK: Remove the cyclic prefix for correlation.
% OK: Correlation plots off by one sample too less?
% OK: Correlera med en symbol bara?
% OK: Something wrong with PL?
% OK: check why the number of taps is 14 for CDL-E.
% OK: ROTATE ARRAY (NOT ONLY ELEMENT). In rayPhaseCenterOffset.
% OK: CHECK WHY TIME EVOLUTION OF CHANNELS DOES NOT SEEM TO WORK.
% Check scaling of channel with the number of antenna elements.
% OK: Apply pathloss.
% OK: Apply LOS phase. NOT NEEDED UNLESS PHASE INDUCED BY MOVEMENT IS TRACKED.
% OK: Implement OFDM waveform. Cell layout: how to assign cell-IDs to the different cells in such a way that intercell interference on the PRS is minimized/maximized?
% OK: Apply beamforming AFTER call to `nrCDLChannel.m'. Figure out which codebook that shoudl be applied.
% OK: Apply delay according to TRP2UT distance.
% Implement callback function for different element responses.
% Consider working in the frequency domain if there is no ISI. Checks!
% Timing considerations of OFDM waveforms from multiple TRPs. What is the general principle?
% How does DTOA really works? Does it need time domain waveforms (in complex baseband representation)?

% Generate UTs and TRPs. Origin (0,0) means (left,bottom) corner of the
% office room.
TRPPositions = [...
    10, 30, 50, 70, 90, 110, 10, 30, 50, 70, 90, 110;...
    15, 15, 15, 15, 15,  15, 35, 35, 35, 35, 35,  35;...
    3,  3,  3,  3,  3,  3,   3,  3,  3,  3,  3,   3];
TRPBearings = [30; 150; -90] - 90; % obs
% OBS, use here d2Dmin.
UTs = struct('ArrayPosition',[rand(1)*120; rand(1)*50; 1.5],'UTDirectionOfTravel',[rand(1)*360 - 180;0],'Mobility',3000/3600);
TRPs = dropTRPs(TRPPositions,TRPBearings,UTs.UTDirectionOfTravel,UTs.Mobility,UTs.ArrayPosition);
figRoom = plotRoom(TRPs,UTs,1.0);
% Select 6 closest TRPs...
TRPs = selectClosestTRPs(TRPs,UTs,6);
% ...and reassign CellIDs.
for TRPIndex = 1:length(TRPs)
    TRPs(TRPIndex).CellID = mod(TRPIndex-1,6);
end
plotRoom(TRPs,UTs,3.0,figRoom,{'1';'2';'3';'4';'5';'6'});

% Loop over TRPs.
for TRPIndex = 1:length(TRPs)
    % Generate the PRS to be transmitted. Symbols between the fourth and
    % the eleventh are allocated for PRS transmission.
    [signalIn,signalInSymbols] = nrPRS(TRPs(TRPIndex).CellID,MaxDLNumberRB,DLNumberRB,PRSNumberRB,slotNumber);
    signalInTRPs(:,TRPIndex) = signalIn;
    % Apply a different precoding codeword per allocated PRS symbol.
    Ns = length(signalIn);
    Nt = length(TRPs(TRPIndex).ArrayCodebook(:,1));
    signalInPrecoded = zeros(Ns,Nt);
    for symbolIndex = 3:10
        symbolStart = signalInSymbols(symbolIndex + 1);
        if (symbolIndex + 1 < length(signalInSymbols))
            symbolEnd = signalInSymbols(symbolIndex + 2) - 1;
        else
            symbolEnd = length(signalIn) - 1;
        end
        codewordIndex = mod(symbolIndex - 3,8);
        signalInPrecoded(symbolStart + 1:symbolEnd + 1,:) = signalIn(symbolStart + 1:symbolEnd + 1)*TRPs(TRPIndex).ArrayCodebook(:,codewordIndex + 1).';
    end
    % 
    [signalOut,pathGains,sampleTimes] = nrCDLChannel(signalInPrecoded,TRPs(TRPIndex).Channel);
    %figure; plot(real(signalIn)); hold on; plot(real(signalOut(:,1))); axis tight; grid on;
    %figure; plot(imag(signalIn)); hold on; plot(imag(signalOut(:,1))); axis tight; grid on;
    d2D = norm(TRPs(TRPIndex).ArrayPosition(1:2) - UTs.ArrayPosition(1:2));
    hBS = TRPs(TRPIndex).ArrayPosition(3);
    hUT = UTs.ArrayPosition(3);
    [PL,sigmaSF] = nrPathloss('InH',TRPs(TRPIndex).Channel.HasLOSCluster,d2D,hBS,hUT,TRPs(TRPIndex).Channel.CarrierFrequency);
    PL = PL + randn(1)*sigmaSF;
    PLTRPs(TRPIndex) = PL;
    signalOut = signalOut/sqrt(power(10,PL/10));
    signalOutTRPs(:,:,TRPIndex) = signalOut;
end

% Time domain signal (real part) before precoding for each of 6 UE closest TRP
figure;
for i=1:6
    y = signalInTRPs(:,i); 
    subplot(2,3,i);%  hold on; 
    plot(real(y)); 
    grid on;
    axis tight;
end

%Time domain signal at UE location (after applied channel) from the 6
%closest TRP
figure; 
for i=1:6
    y = signalOutTRPs(:,1,i); 
    subplot(2,3,i);%  hold on; 
    plot(real(y)); 
    grid on;
    axis tight;
end

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

for i = 1:6
    disp(norm(TRPs(i).ArrayPosition - UTs.ArrayPosition)); 
end

return;

% Loop over TRPs.
%   Compute channel. Apply path loss and shadow fading.
%   Loop over sectors.
%       Compute input waveform, including PRS. Filter with computed channel.
%       Loop TRP/UT over beams.
%           Apply beams. Add waveform.
% Save resultant waveform.
%

function TRPs = dropTRPs(TRPPositions,TRPBearings,UTDirectionOfTravel,Mobility,UTPosition)

% Codebook of beams equispaced in azimuth when using an array defined by
% txArrayFR2.Size = [4 8 1 1 1]; % Given as [M N P Mg Ng].
% txArrayFR2.ElementSpacing = [0.5 0.5 0.0 0.0]; % Given as [dV dH dgV dgH].
phiGCS = repmat(37.5 + 15*(7:-1:0),8,1);
codebookdft = exp(-1i*pi*repmat((0:7).',1,8).*cos(deg2rad(phiGCS)));
%codebookdft = fftshift(dftmtx(8),2);
codebook(1,:,:) = codebookdft/sqrt(8)/sqrt(4);
codebook(2,:,:) = codebookdft/sqrt(8)/sqrt(4);
codebook(3,:,:) = codebookdft/sqrt(8)/sqrt(4);
codebook(4,:,:) = codebookdft/sqrt(8)/sqrt(4);
codebook = reshape(codebook,4*8,8);

% Generate a list of TRPs with the desired positions and bearing angles.
TRPIndex = 1;
for positionIndex = 1:size(TRPPositions,2)
    for sectorIndex = 1:length(TRPBearings)        
        channel = getDefaultChannel(...
            [TRPBearings(sectorIndex); 0; 0],...
            'CDL-E',...%            'LOS',... obs
            30e9,...
            30e9*Mobility/physconst('lightspeed'),...
            UTDirectionOfTravel,...
            30.72e6*8,...
            30720,...
            norm(UTPosition-TRPPositions(:,positionIndex))/physconst('lightspeed'),...
            (UTPosition-TRPPositions(:,positionIndex))/norm(UTPosition-TRPPositions(:,positionIndex)));
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

function channel = getDefaultChannel(TransmitAntennaArrayOrientation,DelayProfile,CarrierFrequency,MaximumDopplerShift,UTDirectionOfTravel,SampleRate,NumTimeSamples,FirstPathDelay,TRP2UTDir)
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
angleAoDGCS = WrapToPlusMinus180(rad2deg(angleAoDGCS));
angleZoDGCS = 90 - rad2deg(angleZoDGCS);
if (angleZoDGCS<0 || angleZoDGCS>180)
    disp('WARNING: thetaLCS out of range after `cart2sph'' conversion.');
end
[angleAoAGCS,angleZoAGCS,~] = cart2sph(-TRP2UTDir(1),-TRP2UTDir(2),-TRP2UTDir(3));
angleAoAGCS = WrapToPlusMinus180(rad2deg(angleAoAGCS));
angleZoAGCS = 90 - rad2deg(angleZoAGCS);
if (angleZoAGCS<0 || angleZoAGCS>180)
    disp('WARNING: thetaLCS out of range after `cart2sph'' conversion.');
end
channel.FirstPathAngleAoA = angleAoAGCS; % In degrees. OBS
channel.FirstPathAngleAoD = angleAoDGCS; % In degrees. OBS
channel.FirstPathAngleZoA = angleZoAGCS; % In degrees. OBS
channel.FirstPathAngleZoD = angleZoDGCS; % In degrees. OBS

function fig = plotRoom(TRPs,UTs,lineWidth,fig,myText)
if (nargin < 3)
    lineWidth = 1.0;
end
if (nargin < 4)
    fig = figure;
else
    figure(fig);
    hold on;
end
R = 3;
hold on;
line([0 120 120 0 0],[0 0 50 50 0]);
for TRPIndex = 1:length(TRPs)
    plot(TRPs(TRPIndex).ArrayPosition(1),TRPs(TRPIndex).ArrayPosition(2),'bo','LineWidth',lineWidth);
    h = line(TRPs(TRPIndex).ArrayPosition(1)+R*[0,cos(deg2rad(TRPs(TRPIndex).ArrayBearing + 90))],TRPs(TRPIndex).ArrayPosition(2)+R*[0,sin(deg2rad(TRPs(TRPIndex).ArrayBearing + 90))]);
    set(h,'Color','r');
    set(h,'LineWidth',lineWidth);
    phiGCS = 37.5 + 15*(7:-1:0);
    if (lineWidth>1)
        for AngleIdx = 1:length(phiGCS)
            hold on;
            h = line(TRPs(TRPIndex).ArrayPosition(1)+2*R*[0,cos(deg2rad(TRPs(TRPIndex).ArrayBearing + phiGCS(AngleIdx)))],TRPs(TRPIndex).ArrayPosition(2)+2*R*[0,sin(deg2rad(TRPs(TRPIndex).ArrayBearing + phiGCS(AngleIdx)))]);
            set(h,'Color','k');
            set(h,'LineWidth',lineWidth/2);
        end
    end
    if (nargin>4)
        text(TRPs(TRPIndex).ArrayPosition(1)+3,TRPs(TRPIndex).ArrayPosition(2),myText{TRPIndex});
    end
end
plot(UTs(1).ArrayPosition(1,:),UTs.ArrayPosition(2,:),'gx');
h = line(UTs.ArrayPosition(1)+R*[0,cos(deg2rad(UTs.UTDirectionOfTravel(1)))],UTs.ArrayPosition(2)+R*[0,sin(deg2rad(UTs.UTDirectionOfTravel(1)))]);
set(h,'Color','k');
set(h,'LineWidth',lineWidth);
axis equal; 
grid on;

function TRPs = selectClosestTRPs(TRPs,UTs,numTRPs)
% Select `numTRPs' closest to UTs(1).
d = zeros(length(TRPs),1);
for TRPIndex = 1:length(TRPs)
    UTBearing = UTs(1).ArrayPosition(1:2) - TRPs(TRPIndex).ArrayPosition(1:2);
    UTBearing = UTBearing/norm(UTBearing);
    TRPBearing = [cos(deg2rad(TRPs(TRPIndex).ArrayBearing + 90)),sin(deg2rad(TRPs(TRPIndex).ArrayBearing + 90))].';
    UT2TRPBearing = WrapToPlusMinus180(rad2deg(acos(UTBearing'*TRPBearing)));
    if (abs(UT2TRPBearing) <= TRPs(TRPIndex).ArrayBeamWidth/2)
        UTVec = UTs(1).ArrayPosition - TRPs(TRPIndex).ArrayPosition;
        d(TRPIndex) = norm(UTVec);
    else
        d(TRPIndex) = Inf;
    end
end
[d,I] = sort(d);
if (length(d) > numTRPs)
    TRPs = TRPs(I(1:numTRPs));
end

function y = WrapToPlusMinus180(x)
% Wrap x to [-180,180).
y = mod(x + 180,2*180) - 180;
