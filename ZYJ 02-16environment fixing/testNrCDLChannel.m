function testNrCDLChannel

signalIn=[]; 

[signalOut,pathGains,sampleTimes] = nrCDLChannel(signalIn);

figure; plot(abs(signalOut(:,1)));

return;

