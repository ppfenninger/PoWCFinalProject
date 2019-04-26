constant = ones(100000, 1)*.5;

tmp = zeros(1,length(constant)*2);
tmp(1:2:end) = real(constant);
tmp(2:2:end) = imag(constant);

f1 = fopen('txFinder.dat', 'wb');
% write the values as a float32
fwrite(f1, tmp, 'float32');
% close the file
fclose(f1);
