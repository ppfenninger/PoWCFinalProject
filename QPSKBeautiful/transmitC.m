knownWNF = wgn(64, 10, 0);
knownWNF = 0.5.*(sign(real(knownWNF)) + 1i*sign(imag(knownWNF)));
knownWN = ifft(knownWNF);
knownWN = [knownWN(end-15:end, :); knownWN];

ts = wgn(64, 1, 0, 'complex');
ts = sign(ts); 
lts = [ts; ts; ts].';
lts = 0.35.*lts; 

save('lts', 'lts'); 
% 
save('knownWN', 'knownWN'); 


dataBitsRaw = sign(wgn(128, 100, 0)); 
dataBits = zeros(64,100);
for m = 1:100
    for n = 1:64
        if isequal(dataBitsRaw((n-1)*2+1:n*2,m), [1;1])
            dataBits(n,m) = 1+1i;
        elseif isequal(dataBitsRaw((n-1)*2+1:n*2,m), [-1;1])
            dataBits(n,m) = -1+1i;
        elseif isequal(dataBitsRaw((n-1)*2+1:n*2,m), [-1;-1])
            dataBits(n,m) = -1-1i;
        elseif isequal(dataBitsRaw((n-1)*2+1:n*2,m), [1;-1])
            dataBits(n,m) = 1-1i;
        end
    end
end
dataBits = dataBits*sqrt(2);
%Add pilots
dataBits(7,:)  =  sqrt(2) + sqrt(2)*1i;
dataBits(21,:) = -sqrt(2) + sqrt(2)*1i;
dataBits(44,:) = -sqrt(2) - sqrt(2)*1i;
dataBits(58,:) =  sqrt(2) - sqrt(2)*1i;
dataBits = dataBits.*0.5;
data = ifft(dataBits);
data = [data(end-15:end, :); data];

save('dataBitsRaw', 'dataBitsRaw'); 

tx = [lts, reshape(knownWN, [1, 800]), reshape(data, [1, 8000])];
tx = tx./rms(tx);
tx = 0.25.*tx; 
tx = [zeros(1, 6400), tx];

save('tx', 'tx'); 

tmp = zeros(1,length(tx)*2);
tmp(1:2:end) = real(tx);
tmp(2:2:end) = imag(tx);

f1 = fopen('tx.dat', 'wb');
% write the values as a float32
fwrite(f1, tmp, 'float32');
% close the file
fclose(f1);