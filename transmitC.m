% knownWNF = wgn(64, 10, 0);
% knownWN = ifft(knownWNF); 
% knownWN = [knownWN(end-15:end, :); knownWN];
% 
% ts = wgn(64, 1, 0, 'complex');
% ts = sign(ts); 
% lts = [ts; ts; ts].';
% lts = 0.45.*lts; 
% 
% save('lts', 'lts'); 
% % 
% save('knownWN', 'knownWN'); 
% 
% dataBits = sign(wgn(64, 100, 0)); 
% data = ifft(dataBits);
% data = [data(end-15:end, :); data];
% 
% save('dataBits', 'dataBits'); 
% 
% tx = [lts, reshape(knownWN, [1, 800]), reshape(data, [1, 8000])];
% tx = tx./rms(tx);
% tx = 0.5.*tx; 
% tx = [zeros(1, 6400), tx];
% 
% save('tx', 'tx'); 
% 
% tmp = zeros(1,length(tx)*2);
% tmp(1:2:end) = real(tx);
% tmp(2:2:end) = imag(tx);
% 
% f1 = fopen('tx.dat', 'wb');
% % write the values as a float32
% fwrite(f1, tmp, 'float32');
% % close the file
% fclose(f1);
