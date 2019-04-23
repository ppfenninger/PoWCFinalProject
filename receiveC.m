% Open the file containing the received samples
f2 = fopen('rx.dat', 'rb');
 
% read data from the file
tmp = fread(f2, 'float32');
rx = tmp(1:2:end)+1i*tmp(2:2:end);
rx = rx.';

% close the file
fclose(f2);

% f1 = fopen('tx.dat', 'rb');
% tx = fread(f1, 'float32'); 

load('tx'); 
disp(length(rx));

load('lts'); 

load('knownWN'); 

[r, lags] = xcorr(rx, lts); 
[~, startIndex] = max(abs(r)); 
startLag = lags(startIndex);
plot(lags, abs(r)); 

lts = rx(startLag:startLag + 191); 
ts1 = lts(64:127); 
ts2 = lts(128:191); 

fDeltaSum = 0; 

for m = 1:64
    fDeltaSum = fDeltaSum + angle(ts2(m)/ts1(m));
end

fDelta = fDeltaSum / (64^2); 
rxLong = rx((startLag+192):(startLag + 8799 + 192)); 

for k = 1:length(rxLong)
    rxLong(k) = rxLong(k)*(exp(-1i*fDelta*k)); 
end

rxR = reshape(rxLong, [80, 110]); 
rxR = rxR(17:80, :); 

rxData = rxR(:, 11:end); 

WN = rxR(:, 1:10);
knownWN = knownWN(17:80, :); 
H = fft(WN)./fft(knownWN); 
H = mean(H, 2); 

dataF = fft(rxData); 
% load('bestH'); 
dataF = dataF./H;

% Use pilots for finer angle offset correction
pilot7  = angle( dataF(7,:) /( 1+1j) );
pilot21 = angle( dataF(21,:)/(-1+1j) );
pilot44 = angle( dataF(44,:)/(-1-1j) );
pilot58 = angle( dataF(58,:)/( 1-1j) );
angleOffset = (pilot7+pilot21+pilot44+pilot58)./4;

dataF = dataF./exp(1j*angleOffset);

dataHat = sign(real(dataF)); 

load('dataBits'); 

% totalWrong = reshape(dataHat, [1, 6400]) ~= reshape(dataBits, [1, 6400]);
dataHat  = dataHat ([1:6,8:20,22:27,39:43,45:57,59:64], 1:95); 
dataBits = dataBits([1:6,8:20,22:27,39:43,45:57,59:64], 1:95); 
% dataHat = dataHat(:, 1:50); 
% dataBits = dataBits(:, 1:50); 
totalWrong = (dataHat ~= dataBits); 
stem(totalWrong)
% plot(reshape(totalWrong, [1, 5300]), '.'); 

bitErrorRate = sum(sum(totalWrong)) / (49*95); 

% stem(reshape(dataHat, [1, 64000]))
% hold on
% stem(reshape(dataBits, [1, 64000]))