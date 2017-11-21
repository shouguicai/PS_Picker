%%
Smooth_window = 50;
env = conv(amp,ones(1,Smooth_window)/Smooth_window);
env1 = conv(env,ones(1,Smooth_window)/Smooth_window);
differ = diff(env1);
signer = sign(differ);
%%
plot(amp);
hold on
plot(env1,'r');
% plot(differ*1e5,'m');
%%
F = 32;
[Y,F1,T1,P1] = spectrogram(differ,F,F/2,F,fs,'yaxis');
%%
figure;
surf(T1,F1,10*log10(abs(P1)),'EdgeColor','none');
axis xy; axis tight; 
colormap(jet); 
view(0,90);
colorbar
xlabel('Time (sec)','FontSize',12);
ylabel('Frequency (Hz)','FontSize',12);
%%
figure;
surf(T1,F1,10*log10(abs(P1)),'EdgeColor','none');
axis xy; axis tight; 
colormap(jet); 
view(0,0);
colorbar
xlabel('Time (sec)','FontSize',12);
ylabel('Frequency (Hz)','FontSize',12);
%%
zz = hilbert(differ*1e5);
figure(2)
subplot(2,1,1)
plot(abs(zz))
ylabel('Amp')
grid on
title('Instantaneous Amp')
instfreq = fs/(2*pi)*diff(unwrap(angle(zz)));
figure(2)
subplot(2,1,2)
plot(instfreq)
ylabel('Hz')
grid on
title('Instantaneous Frequency')
%%
alarm = envelop_hilbert_v2(differ,round(50*fs),1,round(20*fs),1);