function B = de2bi(D)
% de2bi Convert decimal numbers to binary numbers.
% B = de2bi(D) converts a nonnegative integer decimal vector D to a binary
% matrix B. Each row of the binary matrix B corresponds to one element of
% D. The default orientation of the binary output is Right-MSB; the first
% element in B represents the lowest bit.

B = flipud(str2double(dec2bin(D)'));