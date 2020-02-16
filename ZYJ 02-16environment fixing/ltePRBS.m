function [seq, cinit] = ltePRBS(cinit, n)
% [seq, cinit] = ltePRBS(cinit, n) returns the first n elements of the
% pseudorandom binary sequence (PRBS) generator when initialized with
% cinit. For uniformity with the channel specific PRBS functions, ltePRBS
% also returns the initialization value cinit.
% The generation procedure is described in 3GPP 36.211 v14 section 7.2.

Nc = 1600;

% Initialization of x1(0:30), x2(0:30).
x1(0+1) = 1; 
%x1(1+1:30+1,1) = zeros(30,1);
x1(1+1:30+1) = zeros(30,1);
%x2 = flipud(str2num(dec2bin(cinit)'));
x2 = de2bi(cinit);
if (length(x2)>31)
    x2 = x2(0+1:30+1);
else
    x2(end+1:30+1) = zeros(31-length(x2),1);
end

% Generate x1(31:Nc+n-1), x2(31:Nc+n-1).
for nn = 31-31:Nc+n-1-31
    x1(nn+31+1) = mod(x1(nn+3+1) + x1(nn+0+1),2);
    x2(nn+31+1) = mod(x2(nn+3+1) + x2(nn+2+1) + x2(nn+1+1) + x2(nn+0+1),2);
end

% Generate a length-31 Gold sequence, seq(0:n-1).
seq(0+1:n-1+1) = mod(x1(Nc+0+1:Nc+n-1+1)+x2(Nc+0+1:Nc+n-1+1),2);

return;
