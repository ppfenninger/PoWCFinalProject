clear;
clf;

% Open the file containing the received samples
f2 = fopen('rx.dat', 'rb');
%  
% % read data from the file
tmp = fread(f2, 'float32');
rx = tmp(1:2:end)+1i*tmp(2:2:end);
rx = rx.';
rx = rx(10000:end);

% close the file
fclose(f2);

% f1 = fopen('tx.dat', 'rb');
% tx = fread(f1, 'float32'); 

load('tx'); 
% rx = nonflat_channel_timing_error(tx.').'; 
% rx = nonflat_channel_timing_error(tx.').';
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
disp(fDelta);
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

dataHat    = zeros(64*4,100);

dataF = dataF*2/sqrt(2); %Should scale to maximum of 1

for m = 1:100
    for n = 1:64
        if real(dataF(n,m))>2/3
            if imag(dataF(n,m))>2/3
                dataHat((n-1)*4+1:n*4,m)=[1;1;1;1];
            elseif imag(dataF(n,m))>0
                dataHat((n-1)*4+1:n*4,m)=[1;1;1;-1];
            elseif imag(dataF(n,m))>-2/3
                dataHat((n-1)*4+1:n*4,m)=[1;1;-1;-1];
            else
                dataHat((n-1)*4+1:n*4,m)=[1;1;-1;1];
            end
        elseif real(dataF(n,m))>0
            if imag(dataF(n,m))>2/3
                dataHat((n-1)*4+1:n*4,m)=[1;-1;1;1];
            elseif imag(dataF(n,m))>0
                dataHat((n-1)*4+1:n*4,m)=[1;-1;1;-1];
            elseif imag(dataF(n,m))>-2/3
                dataHat((n-1)*4+1:n*4,m)=[1;-1;-1;-1];
            else
                dataHat((n-1)*4+1:n*4,m)=[1;-1;-1;1];
            end
        elseif real(dataF(n,m))> -2/3
            if imag(dataF(n,m))>2/3
                dataHat((n-1)*4+1:n*4,m)=[-1;-1;1;1];
            elseif imag(dataF(n,m))>0
                dataHat((n-1)*4+1:n*4,m)=[-1;-1;1;-1];
            elseif imag(dataF(n,m))>-2/3
                dataHat((n-1)*4+1:n*4,m)=[-1;-1;-1;-1];
            else
                dataHat((n-1)*4+1:n*4,m)=[-1;-1;-1;1];
            end
        else
            if imag(dataF(n,m))>2/3
                dataHat((n-1)*4+1:n*4,m)=[-1;1;1;1];
            elseif imag(dataF(n,m))>0
                dataHat((n-1)*4+1:n*4,m)=[-1;1;1;-1];
            elseif imag(dataF(n,m))>-2/3
                dataHat((n-1)*4+1:n*4,m)=[-1;1;-1;-1];
            else
                dataHat((n-1)*4+1:n*4,m)=[-1;1;-1;1];
            end
        end
    end
end

load('dataBitsRaw'); 

% totalWrong = reshape(dataHat, [1, 6400]) ~= reshape(dataBits, [1, 6400]);
% dataHat  = dataHat       ([1:12, 15:40, 43:54,77:86,89:114,117:128], 1:90); 
% dataBitsRaw = dataBitsRaw([1:12, 15:40, 43:54,77:86,89:114,117:128], 1:90);

bitsWeWant = ones(64,1);
bitsWeWant([7,21,28:38,44,58]) = 0;
bitsWeWant = upsample(bitsWeWant,4);
bitMask = ones(4,1);
bigBitsWeWant = conv(bitsWeWant,bitMask);
bigBitsWeWant = bigBitsWeWant(1:256);

dataHat = dataHat.*bigBitsWeWant;
dataBitsRaw = dataBitsRaw.*bigBitsWeWant;
dataHat( all(~dataHat,2), : ) = [];
dataHat = dataHat(:, 1:80); 
dataBitsRaw( all(~dataBitsRaw,2), : ) = [];
dataBitsRaw = dataBitsRaw(:, 1:80); 



% Saved BPSK data chopping
% dataBits = dataBits([1:6,8:20,22:27,39:43,45:57,59:64], 1:95);

% dataHat = dataHat(:, 1:50); 
% dataBits = dataBits(:, 1:50); 
totalWrong = (dataHat ~= dataBitsRaw); 
x = linspace(1, length(totalWrong), length(totalWrong));
stem(x./4, totalWrong)
% plot(reshape(totalWrong, [1, 5300]), '.'); 

bitErrorRate = sum(sum(totalWrong))/ (49*80*4); 

figure;
plot(abs(H))
disp(bitErrorRate*100)
% plot(abs(H([1:6, 7:20, 22:27, 38:43, 45:57, 59:64])))


% stem(reshape(dataHat, [1, 64000]))
% hold on
% stem(reshape(dataBits, [1, 64000]))