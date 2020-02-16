function [signalOut,signalOutSymbols] = nrPRS(NPRSid,MaxDLNumberRB,DLNumberRB,PRSNumberRB,slotNumber)   %seems like slotNumber doesn't show up in the code...
% nrPRS generates positioning reference signals (PRS) to be sent out on
% antenna port 6. The function implements the following aspects of 
% TR 36.211:
%Section 6.10.4: Positioning reference signals
%Section 7.2: Pseudo-random sequence generation
%
% NPRSid, is the PRS id = 0,...,4095, which equals the cell id unless
% configured by higher layers.
% signalOut, a complex vector of size Ns x Nsymb, where Ns is the number of
% subcarriers, rounded to the nearest non-smaller power of 2, Nsymb is 14,
% the number of symbols per slot. signalOut is normalized to havinng unit
% variance, when accounting for the duty cycle.
% signalOutSymbols is an Nsymb x 1 vector whose k:th entry corresponds to
% the first sample 0..Ns-1 of symbol k.
%
% The layout used for NR is the same as specified for LTE for the case of 1
% or 2 PBCH antenna ports. 

if ((NPRSid < 0) || (NPRSid > 4095))
    disp('WARNING: `NPRSid'' out of validity range in call to `nrPRS''');
end
if ((slotNumber < 0) || (slotNumber > 19))
    disp('WARNING: `slotNumber'' out of validity range in call to `nrPRS''');
end

% Implement TR 36.211, Sec. 6.10.4.2. Assume NCP, 1 or 2 PBCH antenna
% ports, and no frequency hoping.
Ns = power(2,ceil(log2(12*MaxDLNumberRB)));                                %Zhang: Ns<12*MaxDLNumberRB and Ns is 2 to power of an integer. Is MaxDLNumberRB a constant? (110)
m = 0:2*PRSNumberRB-1;              
mprim = m + MaxDLNumberRB - PRSNumberRB;
signalOut = zeros(Ns,14);
signalOutSymbols = zeros(14,1);
symbolNumberList = [3,4,5,6,7,8,9,10];    
%symbolNumberList = [3,5,6,8,9,10,12,13];
nuShift = mod(NPRSid,6); 
for symbolNumber = symbolNumberList                               
    ell = mod(symbolNumber,7);                                           
    k = round(6*(m + Ns/12 - PRSNumberRB) + mod(6 - ell - nuShift,6));     % zhang: Ns/12=MaxDLNumberRB ?
%    k = 6*(m + DLNumberRB - PRSNumberRB) + mod(6 - ell - nuShift,6);
    r = nrPRBS(slotNumber,symbolNumber,NPRSid,MaxDLNumberRB);           
    signalOut(k+1,symbolNumber+1) = r(mprim);
end

% Prepare for IFFT...
signalOut = ifftshift(signalOut,1);                      
% ... and avoid DC subcarrier.
signalOut(1+1:6*PRSNumberRB+1,:) = signalOut(1:6*PRSNumberRB,:);  
signalOut(1,:) = 0;

% Do IFFT...
signalOutTmp = ifft(signalOut,[],1);
% ... and add CP.
signalOut = [];
CPLength = [160,144,144,144,144,144,144,160,144,144,144,144,144,144];   
sampleIndex = 0;                                                          
for symbolIndex = 0:13                                                 
    signalOutSymbols(symbolIndex+1) = sampleIndex+1;
    signalOut(sampleIndex+1:sampleIndex+CPLength(symbolIndex+1)) = signalOutTmp(end-CPLength(symbolIndex+1)+1:end,symbolIndex+1); 
    sampleIndex = sampleIndex + CPLength(symbolIndex+1);
    signalOut(sampleIndex+1:sampleIndex+2048) = signalOutTmp(:,symbolIndex+1);
    sampleIndex = sampleIndex + 2048;                                      %zhang: ....
end

signalOut = signalOut/sqrt(var(signalOut)*14/length(symbolNumberList));    %zhang: Normalization. seems like the some of the symbols are excluded (14/8)

%var(signalOut)
%figure; plot(real(signalOut)); axis tight; grid on;
signalOut = signalOut(:);

return;

function r = nrPRBS(slotNumber,symbolNumber,NPRSid,MaxDLNumberRB)
% Implement TR 36.211, Sec. 6.10.4.1.
NCP = 1; % Obs! Only normal CP supported.
ns = 2*slotNumber + floor(symbolNumber/7); % Map NR to LTE subframe number.
ell = mod(symbolNumber,7); % Map NR to LTE symbol number
cinit = 2^28*floor(NPRSid/512) + 2^10*(7*(ns+1)+ell+1)*(2*mod(NPRSid,512)+1) + 2*mod(NPRSid,512) + NCP;  % Zhang: looks different from 36.211
c = ltePRBS(cinit, 4*MaxDLNumberRB);                                       %(PRBS) generator
c = reshape(c,2,2*MaxDLNumberRB).';
r = 1/sqrt(2)*(1-2*c(:,1)) + 1i/sqrt(2)*(1-2*c(:,2));
