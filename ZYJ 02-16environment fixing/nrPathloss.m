function [PL,sigmaSF] = nrPathloss(scenario,hasLOS,d2D,hBS,hUT,fc)
% nrPathloss computes the distance dependence pathloss. The function
% implements the following aspects of TR 38.901 [1]:
%Section 7.4.1: pathloss
%
% d, a real scalar  
% signalOut, a complex vector of size Ns x Nr, where Ns is the number of
% samples, and Nr is the number of receive antennas. 
% pathGains, a complex vector of size Ncs x Np x Nt x Nr, where Ncl is the
% number of clusters, and Np is the number of paths.
% sampleTimes, sample times of channel snapshots, a real column vector of
% size Ncs x 1.
%
% Note that 38l.901 adopts the points of view of the DL. Thus below the BS
% acts as the transmitter, and the UT as the receiver. Of course, the UL
% channel is simple the reciprocal of the DL channel.

switch lower(scenario)
    case 'rma'
        h = 5; % 5 <= h <= 50 is the avg. building height.
        W = 20; % 5 <= W <= 50 is the avg. street width.
        if ((hUT < 1) || (hUT > 10))
            disp('WARNING: `hUT'' outside range of validity in call to `nrPathloss''');
        end
        if ((hBS < 10) || (hBS > 150))
            disp('WARNING: `hBS'' outside range of validity in call to `nrPathloss''');
        end
        if ((fc < .5e9) || (fc > 30e9))
            disp('WARNING: `fc'' outside range of validity in call to `nrPathloss''');
        end
        if (hasLOS)
            if ((d2D < 10) || (d2D > 10000))
                disp('WARNING: `d2D'' outside range of validity in call to `nrPathloss''');
                PL = NaN;
                sigmaSF = NaN;
            else
                [PL,sigmaSF] = nrPathlossRMaLOS(d2D,hBS,hUT,fc,h);
            end
        else
            if ((d2D < 10) || (d2D > 5000))
                disp('WARNING: `d2D'' outside range of validity in call to `nrPathloss''');
                PL = NaN;
                sigmaSF = NaN;
            else
                [PLLOS,~] = nrPathlossRMaLOS(d2D,hBS,hUT,fc,h);
                [PLNLOS,sigmaSF] = nrPathlossRMaNLOS(d2D,hBS,hUT,fc,h,W);
                PL = max(PLLOS,PLNLOS);
            end
        end
    case 'uma'
        if ((hUT < 1.5) || (hUT > 22.5))
            disp('WARNING: `hUT'' outside range of validity in call to `nrPathloss''');
        end
        if (hBS ~= 25)
            disp('WARNING: `hBS'' outside range of validity in call to `nrPathloss''');
        end
        if ((fc < .5e9) || (fc > 100e9))
            disp('WARNING: `fc'' outside range of validity in call to `nrPathloss''');
        end
        if (hasLOS)
            if ((d2D < 10) || (d2D > 5000))
                disp('WARNING: `d2D'' outside range of validity in call to `nrPathloss''');
                PL = NaN;
                sigmaSF = NaN;
            else
                PL = nrPathlossUMaLOS(d2D,hBS,hUT,fc);
                sigmaSF = 4;
            end
        else
            if ((d2D < 10) || (d2D > 5000))
                disp('WARNING: `d2D'' outside range of validity in call to `nrPathloss''');
                PL = NaN;
                sigmaSF = NaN;
            else
                PLLOS = nrPathlossUMaLOS(d2D,hBS,hUT,fc);
                PLNLOS = nrPathlossUMaNLOS(d2D,hBS,hUT,fc);
                PL = max(PLLOS,PLNLOS);
                sigmaSF = 6;                
                % PL = 32.4 + 20*log10(fc/1e9) + 30*log10(d3D); % optional.
                % sigmaSF = 7.8; % optional.
            end
        end
    case 'umi'
        if ((hUT < 1.5) || (hUT > 22.5))
            disp('WARNING: `hUT'' outside range of validity in call to `nrPathloss''');
        end
        if (hBS ~= 10)
            disp('WARNING: `hBS'' outside range of validity in call to `nrPathloss''');
        end
        if ((fc < .5e9) || (fc > 100e9))
            disp('WARNING: `fc'' outside range of validity in call to `nrPathloss''');
        end
        if (hasLOS)
            if ((d2D < 10) || (d2D > 5000))
                disp('WARNING: `d2D'' outside range of validity in call to `nrPathloss''');
                PL = NaN;
                sigmaSF = NaN;
            else
                PL = nrPathlossUMiLOS(d2D,hBS,hUT,fc);
                sigmaSF = 4;
            end
        else
            if ((d2D < 10) || (d2D > 5000))
                disp('WARNING: `d2D'' outside range of validity in call to `nrPathloss''');
                PL = NaN;
                sigmaSF = NaN;
            else
                PLLOS = nrPathlossUMiLOS(d2D,hBS,hUT,fc);
                PLNLOS = nrPathlossUMiNLOS(d2D,hBS,hUT,fc);
                PL = max(PLLOS,PLNLOS);
                sigmaSF = 7.82;
                % PL = 32.4 + 20*log10(fc/1e9) + 31.9*log10(d3D); % optional.
                % sigmaSF = 8.2; % optional.
            end
        end
    case 'inh'
        if ((fc < .5e9) || (fc > 100e9))
            disp('WARNING: `fc'' outside range of validity in call to `nrPathloss''');
        end
        d3D = sqrt(d2D^2 + (hBS-hUT)^2);
        if ((d3D < 1) || (d3D > 150))
            disp('WARNING: `d3D'' outside range of validity in call to `nrPathloss''');
            PL = NaN;
            sigmaSF = NaN;
        else
            if (hasLOS)
                PL = nrPathlossInHLOS(d3D,fc);
                sigmaSF = 3;
            else
                PLLOS = nrPathlossInHLOS(d3D,fc);
                PLNLOS = nrPathlossInHNLOS(d3D,fc);
                PL = max(PLLOS,PLNLOS);
                sigmaSF = 8.03;
                % PL = 32.4 + 17.3*log10(d3D) + 20*log10(fc/1e9); % optional.
                % sigmaSF = 8.29; % optional.
            end
        end
    otherwise
        disp('WARNING: Unkown scenario type in call to `nrPathloss''');
        PL = NaN;
        sigmaSF = NaN;
end

function [PL,sigmaSF] = nrPathlossRMaLOS(d2D,hBS,hUT,fc,h)
    d3D = sqrt(d2D^2 + (hBS-hUT)^2);
    dBP = nrBreakPointDistanceRMa(hBS,hUT,fc);
    if ((d2D >= 10) && (d2D <= dBP))
        PL = 20*log10(40*pi*d3D*fc/1e9/3) + min(.03*h^1.72,10)*log10(d3D) - min(.044*h^1.72,14.77)+ .002*log10(h)*d3D;
        sigmaSF = 4;
    else
        PL = 20*log10(40*pi*d3D*fc/1e9/3) + min(.03*h^1.72,10)*log10(d3D) - min(.044*h^1.72,14.77)+ .002*log10(h)*d3D;
        PL = PL + 40*log10(d3D/dBP);
        sigmaSF = 6;
    end

function [PL,sigmaSF] = nrPathlossRMaNLOS(d2D,hBS,hUT,fc,h,W)
    d3D = sqrt(d2D^2 + (hBS-hUT)^2);
    PL = 161.04 - 7.1*log10(W) + 7.5*log10(h) - (24.37-3.7*(h/hBS)^2)*log10(hBS) + (43.42-3.1*log10(hBS))*(log10(d3D)-3) + 20*log10(fc/1e9) - (3.2*(log10(11.75*hUT)^2)-4.97);
    sigmaSF = 8;

function dBP = nrBreakPointDistanceRMa(hBS,hUT,fc)
dBP = 2*pi*hBS*hUT*fc/physconst('lightspeed');

function PL = nrPathlossUMaLOS(d2D,hBS,hUT,fc)
    d3D = sqrt(d2D^2 + (hBS-hUT)^2);
    dBP = nrBreakPointDistanceUMa(hBS,hUT,fc);
    if ((d2D >= 10) && (d2D <= dBP))
        PL = 28.0 + 22*log10(d3D) + 20*log10(fc/1e9);
    else
        PL = 28.0 + 40*log10(d3D) + 20*log10(fc/1e9) - 9*log10(dBP^2+(hBS-hUT)^2);
    end

function PL = nrPathlossUMaNLOS(d2D,hBS,hUT,fc)
    d3D = sqrt(d2D^2 + (hBS-hUT)^2);
    PL = 13.54 + 39.08*log10(d3D) + 20*log10(fc/1e9) - .6*(hUT-1.5);

function dBP = nrBreakPointDistanceUMa(hBS,hUT,fc)
if (d2D <= 18)
    g = 0;
else
    g = 5/4*(d2D/100)^3*exp(-d2D/150);
end
if (hUT < 13)
    C = 0;
else
    C = ((hUT-13)/10)^1.5*g;
end
if (rand(1) < 1/(1+C))
    hE = 1;
else
    hEPoints = 12:3:hUT-1.5;
    hE = hEPoints(randi(length(hEPoints)));
end
dBP = 4*(hBS-hE)*(hUT-hE)*fc/physconst('lightspeed');

function PL = nrPathlossUMiLOS(d2D,hBS,hUT,fc)
    d3D = sqrt(d2D^2 + (hBS-hUT)^2);
    dBP = nrBreakPointDistanceUMi(hBS,hUT,fc);
    if ((d2D >= 10) && (d2D <= dBP))
        PL = 32.4 + 21*log10(d3D) + 20*log10(fc/1e9);
    else
        PL = 32.4 + 40*log10(d3D) + 20*log10(fc/1e9) - 9.5*log10(dBP^2+(hBS-hUT)^2);
    end

function PL = nrPathlossUMiNLOS(d2D,hBS,hUT,fc)
    d3D = sqrt(d2D^2 + (hBS-hUT)^2);
    PL = 22.4 + 35.3*log10(d3D) + 21.3*log10(fc/1e9) - .3*(hUT-1.5);
   
function dBP = nrBreakPointDistanceUMi(hBS,hUT,fc)
hE = 1;
dBP = 4*(hBS-hE)*(hUT-hE)*fc/physconst('lightspeed');

function PL = nrPathlossInHLOS(d3D,fc)
    PL = 32.4 + 17.3*log10(d3D) + 20*log10(fc/1e9);

function PL = nrPathlossInHNLOS(d3D,fc)
    PL = 17.3 + 38.3*log10(d3D) + 24.9*log10(fc/1e9);
