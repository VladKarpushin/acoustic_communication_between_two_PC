%Bin signal to binary file

function [errmsg] = signal2file(FILENAME, signal)
% input:
% FILENAME  - name of file to open
% signal    - binary signal [-1 1 -1 1]
%
% output:
% errmsg - returns a system-dependent error message if fopen fails to open the file. Otherwise, errmsg is an empty string

%c_est = num2str(d);        %Convert number to string

d = boolean(signal + 1);      %[0 1 0 1]
c_est = dec2bin(d);         %Convert decimal to binary number in string
b_est = sscanf(c_est', '%8u');
b_est_str = num2str(b_est); %Convert number to string
a_est = bin2dec(b_est_str); %Convert binary number string to decimal number

%[a' a_est]

[fileID,errmsg] = fopen(FILENAME, 'w');
if length(errmsg) ~= 0
    return;
end
fwrite(fileID, a_est);   %fwrite(fileID,A) write the elements of array A as 8-bit unsigned integers to a binary file in column order
fclose(fileID);
