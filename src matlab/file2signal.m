%Bin file to binary signal

function [signal errmsg] = File2signal(FILENAME)
% input:
% FILENAME - name of file to open
%
% output:
% errmsg - returns a system-dependent error message if fopen fails to open the file. Otherwise, errmsg is an empty string
% signal - binary signal [-1 1 -1 1]
signal = 0;
[fileID,errmsg] = fopen(FILENAME,'r');
if length(errmsg) ~= 0
    return;
end

a = fread(fileID,'uint8')';
fclose(fileID);

b = dec2bin(a,8)';  %Convert decimal to binary number in string. Produces a binary representation with at least n bits.
c = strcat(b(:));   %Concatenate strings horizontally
d = bin2dec(c);     %Convert binary number string to decimal number
signal = d*2-1;     %[-1 1 -1 1]
