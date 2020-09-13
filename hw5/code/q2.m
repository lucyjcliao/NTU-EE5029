close all;
clear;
clc

%% Simulation
for Q = 1:3
    tic
    fprintf('[Quantization level = %d]\n',2^Q)
    fprintf('SNR')
    frame = 0;
    BER = [];
    
    for snr = 2:0.5:7           % SNR increment
        fprintf('->%0.1f',snr)
        err_sum = 0;
        
        while err_sum < 1000
            frame = frame+1;
            u = randi(2,1,1000)-1;
            c = Encoder(u);
            s = Modulation(c);
            r = AWGNchannel(s,snr);
            q = Quantization(r,Q);
            u_decode = Decoder(q,12);
            bit_error = ErrorCalculate(u,u_decode);
            err_sum = err_sum+bit_error;
        end
        BER(end+1) = err_sum/(1000*frame);
    end
    fprintf('\n')
    total_BER(Q,:) = BER;
    toc
    fprintf('\n')
end


%% plot
SNR_axis=repmat(2:0.5:7,3,1);
semilogy(2:0.5:7,total_BER,'-o')
axis([2 7 10^(-9) 1])
title('(2,1,4) convolution code with Q =2,4,8 and truncation length=12');
legend('Q=2','Q=4','Q=8');
grid on
xlabel('Eb/No')
ylabel('BER')
set(gcf, 'Position',  [400, 400, 700, 400])
saveas(gcf,'BER(b).jpg')


%% print
fileID = fopen('result(b).txt','w');
fprintf(fileID,'(2,1,4) convolution code with Q =2,4,8 and truncation length=12\n');
for i = 1:size(total_BER,2)
    format_2 = ['SNR = %.1f(dB), BER = ' repmat('%0.15f/',1,size(total_BER,1)) '\n'];
    fprintf(fileID,format_2,(i+3)/2,total_BER(:,i));
end

