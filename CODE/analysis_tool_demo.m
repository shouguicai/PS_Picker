% clear all;
warning off all
%%
root = '../after/';
file_list = dir(root);
%%
sel_num = 132;
sel_index = 5;
imlist = file_list(2+sel_num).name
%%
file = [root,imlist];
X=rdsac(file);
X_d = X.d;
%% calu arrival time for P and S wave
t_start = X.t(1);
t0 = X.t(1) - X.HEADER.B/86400;
x_t = (X.t-t0)*86400-X.HEADER.B;
%%
fs = 1/X.HEADER.DELTA;  % Hz
%%
z = hilbert(X_d);
amp = abs(z);

threshold1 = 5*mean(amp); % need to be adaptive
trigger = 1e5*(sign(amp-threshold1)+1);
trigger_index = find(trigger>0);  
if size(trigger_index,1) == 0
    error('no wave found');   
end
tmp = trigger;
differ = diff(trigger_index);
tmp(trigger_index(1)) = 1;
threshold2 = 5e3;  % signal palse: 50s ,fixed
tmp(trigger_index(2:end)) = sign(differ-threshold2)+1;
tmp_index = find(tmp>0); 
%%

start = max(1,floor(tmp_index(sel_index)-45/X.HEADER.DELTA));
ends = min(start + 4500*2,length(X_d));
XX_d = X_d(start:ends);
XX_t = x_t(start:ends);
%%
figure(sel_num)
subplot(3,1,1)
plot(XX_t,XX_d)
xlim = [min(XX_t),max(XX_t)];
set(gca,'XLim',xlim)
grid on 
xlabel('Time(sec)')
ylabel('Count')
%%
% % F = 32;
% % [Y,F1,T1,P1] = spectrogram(XX_d,F,F/2,F,fs,'yaxis');
% % %%
% % figure;
% % surf(T1,F1,10*log10(abs(P1)),'EdgeColor','none');
% % axis xy; axis tight; 
% % colormap(jet); 
% % view(0,90);
% % colorbar
% % % caxis([-120,-50]);
% % xlabel('Time (sec)','FontSize',12);
% % ylabel('Frequency (Hz)','FontSize',12);
%%
% % figure;
% % surf(T1,F1,10*log10(abs(P1)),'EdgeColor','none');
% % axis xy; axis tight; 
% % colormap(jet); 
% % view(0,0);
% % colorbar
% % % caxis([-120,-50]);
% % xlabel('Time (sec)','FontSize',12);
%%
z = hilbert(XX_d);
amp = abs(z);
figure(sel_num)
subplot(3,1,2)
plot(XX_t,amp)
xlim = [min(XX_t),max(XX_t)];
set(gca,'XLim',xlim)
xlabel('Time(sec)')
ylabel('Amp')
grid on
title('Instantaneous Amp')
instfreq = fs/(2*pi)*diff(unwrap(angle(z)));
figure(sel_num)
subplot(3,1,3)
plot(XX_t(2:end),instfreq)
xlabel('Time(sec)')
ylabel('Hz')
xlim = [min(XX_t),max(XX_t)];
set(gca,'XLim',xlim)
grid on
title('Instantaneous Frequency')
%% Run PphasePicker with optional picking parameters
type = 'na'; % no bandpass filtering  
pflag = 'Y'; % To plot waveform and P-phase onset

% USE Tn = 0.1 for records with low sampling rate lower than 100 sps 
Tn = 0.1;       % undamped natural period of oscillator in second 
xi = 0.6;       % damping ratio 
nbins = 200;    % histogram bin size
o = 'to_peak' ;  % 'to_peak' to take segment of waveform from beginning to
                % absolute peak value
Hfilter = noisefilter;
% XX_d = filter(Hfilter,XX_d);
[loc, snr_db] = PphasePicker(XX_d, X.HEADER.DELTA, type, pflag, Tn, xi, nbins, o);

%% STA
    % define the window
    wintype = 'rectwin';
    winlen = 51;
    winamp = [0.5,1]*(1/winlen);
    % find the zero-crossing rate
    E = energy(XX_d.',wintype,winamp(2),winlen);
    % time index for the ST-ZCR and STE after delay compensation
    N = length(XX_d); % signal length
    out = (winlen-1)/2:(N+winlen-1)-(winlen-1)/2;
    t = (out-(winlen-1)/2)*(1/fs);
    hold on;
    plot(t,2*max(XX_d)*E(out)/max(E(out)),'b-.','Linewidth',2); xlabel('t, seconds');
    if loc ~= -1
        loc_S = pickS(loc,E,out,X.HEADER.DELTA);
%         loc_S = s_catch(XX_d);
        line([loc_S loc_S]*X.HEADER.DELTA,ylim,'Color','m','LineWidth',2);
        ind = aic_pick(XX_d,'whole');
        line([ind ind]*X.HEADER.DELTA,[ylim],'Color','b','LineWidth',2);
        hold off
    end  
    figure(sel_num)
    subplot(3,1,1)
    line(XX_t(1)+[loc_S loc_S]*X.HEADER.DELTA,ylim,'Color','m','LineWidth',2);
    line(XX_t(1)+[loc loc],ylim,'Color','r','LineWidth',2);
    line(XX_t(1)+[ind ind]*X.HEADER.DELTA,[ylim],'Color','b','LineWidth',2);
    figure(sel_num)
    subplot(3,1,2)
    line(XX_t(1)+[loc_S loc_S]*X.HEADER.DELTA,ylim,'Color','m','LineWidth',2);
    line(XX_t(1)+[loc loc],ylim,'Color','r','LineWidth',2);
    line(XX_t(1)+[ind ind]*X.HEADER.DELTA,[ylim],'Color','b','LineWidth',2);
    figure(sel_num)
    subplot(3,1,3)
    line(XX_t(1)+[loc_S loc_S]*X.HEADER.DELTA,ylim,'Color','m','LineWidth',2);
    line(XX_t(1)+[loc loc],ylim,'Color','r','LineWidth',2);
    line(XX_t(1)+[ind ind]*X.HEADER.DELTA,[ylim],'Color','b','LineWidth',2);

%********************************************************************************%
%%
sel_num = sel_num -1;
imlist = file_list(2+sel_num).name
%%
file = [root,imlist];
X=rdsac(file);
X_d = X.d;
%% calu arrival time for P and S wave
t_start = X.t(1);
t0 = X.t(1) - X.HEADER.B/86400;
x_t = (X.t-t0)*86400-X.HEADER.B;
%%
fs = 1/X.HEADER.DELTA;  % Hz
%%
XX_d = X_d(start:ends);
XX_t = x_t(start:ends);
%%
figure(sel_num)
subplot(3,1,1)
plot(XX_t,XX_d)
xlim = [min(XX_t),max(XX_t)];
set(gca,'XLim',xlim)
grid on 
xlabel('Time(sec)')
ylabel('Count')
%%
z = hilbert(XX_d);
amp = abs(z);
figure(sel_num)
subplot(3,1,2)
plot(XX_t,amp)
xlim = [min(XX_t),max(XX_t)];
set(gca,'XLim',xlim)
xlabel('Time(sec)')
ylabel('Amp')
grid on
title('Instantaneous Amp')
instfreq = fs/(2*pi)*diff(unwrap(angle(z)));
figure(sel_num)
subplot(3,1,3)
plot(XX_t(2:end),instfreq)
xlabel('Time(sec)')
ylabel('Hz')
xlim = [min(XX_t),max(XX_t)];
set(gca,'XLim',xlim)
grid on
title('Instantaneous Frequency')
%% Run PphasePicker with optional picking parameters
type = 'na'; % no bandpass filtering  
pflag = 'Y'; % To plot waveform and P-phase onset

% USE Tn = 0.1 for records with low sampling rate lower than 100 sps 
Tn = 0.1;       % undamped natural period of oscillator in second 
xi = 0.6;       % damping ratio 
nbins = 200;    % histogram bin size
o = 'to_peak' ;  % 'to_peak' to take segment of waveform from beginning to
                % absolute peak value
Hfilter = noisefilter;
% XX_d = filter(Hfilter,XX_d);
[loc, snr_db] = PphasePicker(XX_d, X.HEADER.DELTA, type, pflag, Tn, xi, nbins, o);

%% STA
    % define the window
    wintype = 'rectwin';
    winlen = 51;
    winamp = [0.5,1]*(1/winlen);
    % find the zero-crossing rate
    E = energy(XX_d.',wintype,winamp(2),winlen);
    % time index for the ST-ZCR and STE after delay compensation
    N = length(XX_d); % signal length
    out = (winlen-1)/2:(N+winlen-1)-(winlen-1)/2;
    t = (out-(winlen-1)/2)*(1/fs);
    hold on;
    plot(t,2*max(XX_d)*E(out)/max(E(out)),'b-.','Linewidth',2); xlabel('t, seconds');
    if loc ~= -1
        loc_S = pickS(loc,E,out,X.HEADER.DELTA);
%         loc_S = s_catch(XX_d);
        line([loc_S loc_S]*X.HEADER.DELTA,ylim,'Color','m','LineWidth',2);
        ind = aic_pick(XX_d,'whole');
        line([ind ind]*X.HEADER.DELTA,ylim,'Color','b','LineWidth',2);
        hold off
    end  
    figure(sel_num)
    subplot(3,1,1)
    line(XX_t(1)+[loc_S loc_S]*X.HEADER.DELTA,ylim,'Color','m','LineWidth',2);
    line(XX_t(1)+[loc loc],ylim,'Color','r','LineWidth',2);
    line(XX_t(1)+[ind ind]*X.HEADER.DELTA,ylim,'Color','b','LineWidth',2);
    figure(sel_num)
    subplot(3,1,2)
    line(XX_t(1)+[loc_S loc_S]*X.HEADER.DELTA,ylim,'Color','m','LineWidth',2);
    line(XX_t(1)+[loc loc],ylim,'Color','r','LineWidth',2);
    line(XX_t(1)+[ind ind]*X.HEADER.DELTA,ylim,'Color','b','LineWidth',2);
    figure(sel_num)
    subplot(3,1,3)
    line(XX_t(1)+[loc_S loc_S]*X.HEADER.DELTA,ylim,'Color','m','LineWidth',2);
    line(XX_t(1)+[loc loc],ylim,'Color','r','LineWidth',2);
    line(XX_t(1)+[ind ind]*X.HEADER.DELTA,ylim,'Color','b','LineWidth',2);
    
    %********************************************************************************%
%%
sel_num = sel_num -1;
imlist = file_list(2+sel_num).name
%%
file = [root,imlist];
X=rdsac(file);
X_d = X.d;
%% calu arrival time for P and S wave
t_start = X.t(1);
t0 = X.t(1) - X.HEADER.B/86400;
x_t = (X.t-t0)*86400-X.HEADER.B;
%%
fs = 1/X.HEADER.DELTA;  % Hz
%%
XX_d = X_d(start:ends);
XX_t = x_t(start:ends);
%%
figure(sel_num)
subplot(3,1,1)
plot(XX_t,XX_d)
xlim = [min(XX_t),max(XX_t)];
set(gca,'XLim',xlim)
grid on 
xlabel('Time(sec)')
ylabel('Count')
%%
z = hilbert(XX_d);
amp = abs(z);
figure(sel_num)
subplot(3,1,2)
plot(XX_t,amp)
xlim = [min(XX_t),max(XX_t)];
set(gca,'XLim',xlim)
xlabel('Time(sec)')
ylabel('Amp')
grid on
title('Instantaneous Amp')
instfreq = fs/(2*pi)*diff(unwrap(angle(z)));
figure(sel_num)
subplot(3,1,3)
plot(XX_t(2:end),instfreq)
xlabel('Time(sec)')
ylabel('Hz')
xlim = [min(XX_t),max(XX_t)];
set(gca,'XLim',xlim)
grid on
title('Instantaneous Frequency')
%% Run PphasePicker with optional picking parameters
type = 'na'; % no bandpass filtering  
pflag = 'Y'; % To plot waveform and P-phase onset

% USE Tn = 0.1 for records with low sampling rate lower than 100 sps 
Tn = 0.1;       % undamped natural period of oscillator in second 
xi = 0.6;       % damping ratio 
nbins = 200;    % histogram bin size
o = 'to_peak' ;  % 'to_peak' to take segment of waveform from beginning to
                % absolute peak value
Hfilter = noisefilter;
% XX_d = filter(Hfilter,XX_d);
[loc, snr_db] = PphasePicker(XX_d, X.HEADER.DELTA, type, pflag, Tn, xi, nbins, o);

%% STA
    % define the window
    wintype = 'rectwin';
    winlen = 51;
    winamp = [0.5,1]*(1/winlen);
    % find the zero-crossing rate
    E = energy(XX_d.',wintype,winamp(2),winlen);
    % time index for the ST-ZCR and STE after delay compensation
    N = length(XX_d); % signal length
    out = (winlen-1)/2:(N+winlen-1)-(winlen-1)/2;
    t = (out-(winlen-1)/2)*(1/fs);
    hold on;
    plot(t,2*max(XX_d)*E(out)/max(E(out)),'b-.','Linewidth',2); xlabel('t, seconds');
    if loc ~= -1
        loc_S = pickS(loc,E,out,X.HEADER.DELTA);
%         loc_S = s_catch(XX_d);
        line([loc_S loc_S]*X.HEADER.DELTA,ylim,'Color','m','LineWidth',2);
        ind = aic_pick(XX_d,'whole');
        line([ind ind]*X.HEADER.DELTA,[ylim],'Color','b','LineWidth',2);
        hold off
    end  
    figure(sel_num)
    subplot(3,1,1)
    line(XX_t(1)+[loc_S loc_S]*X.HEADER.DELTA,ylim,'Color','m','LineWidth',2);
    line(XX_t(1)+[loc loc],ylim,'Color','r','LineWidth',2);
    line(XX_t(1)+[ind ind]*X.HEADER.DELTA,[ylim],'Color','b','LineWidth',2);
    figure(sel_num)
    subplot(3,1,2)
    line(XX_t(1)+[loc_S loc_S]*X.HEADER.DELTA,ylim,'Color','m','LineWidth',2);
    line(XX_t(1)+[loc loc],ylim,'Color','r','LineWidth',2);
    line(XX_t(1)+[ind ind]*X.HEADER.DELTA,[ylim],'Color','b','LineWidth',2);
    figure(sel_num)
    subplot(3,1,3)
    line(XX_t(1)+[loc_S loc_S]*X.HEADER.DELTA,ylim,'Color','m','LineWidth',2);
    line(XX_t(1)+[loc loc],ylim,'Color','r','LineWidth',2);
    line(XX_t(1)+[ind ind]*X.HEADER.DELTA,[ylim],'Color','b','LineWidth',2);
%%
F = 16;
[Y,F1,T1,P1] = spectrogram(XX_d,F,F/2,F,fs,'yaxis');
% figure;
% surf(T1,F1,(abs(P1)),'EdgeColor','none');
% axis xy; axis tight; 
% colormap(jet); 
% view(0,90);
% colorbar
% % caxis([-120,-50]);
% xlabel('Time (sec)','FontSize',12);
% ylabel('Frequency (Hz)','FontSize',12);

% max_p = max(P1.');
% norm_P = P1./(repmat(max_p.',1,length(T1)));
% tmp = norm_P(1:4,:);
% sum_tmp = sum(tmp,1);
% figure;
% plot(T1,sum_tmp);
%%
max_p = max(P1);
% max_p = P1(1,:);
figure(9);
plot(T1,max_p);hold on
% line(xlim,4*[mean(max_p),mean(max_p)],'Color','g','LineWidth',2);
% sta = 1;
% loc_s = pickS(sta,max_p,sta:length(max_p),1);
delta_T = T1(2) - T1(1);
min_pos = ceil(ind*X.HEADER.DELTA/delta_T);
[~,max_pos] = max(max_p);
if max_pos <= min_pos
    max_pos = max_pos + ceil(5/delta_T);
end
sta = min_pos;
loc_s = aic_pick(max_p(sta:max_pos),'whole')+sta;
line([T1(loc_s),T1(loc_s)],ylim,'Color','y','LineWidth',2);
figure(sel_num)
subplot(3,1,1)
line(XX_t(1)+[T1(loc_s) T1(loc_s)],ylim,'Color','y','LineWidth',2);
%% Test for Kurtosis©\Based method
FB = [1 5;6 10;11 15;16 20];
T = [1 2 3];
n_smooth = 10;
start_interst = 1;
end_interst = floor(3*length(XX_d)/4);
[ind_P, ind_S] = Pick_CF(XX_d,fs,FB,T,n_smooth,start_interst,end_interst);
figure(66);
line([ind_P ind_P]*X.HEADER.DELTA,ylim,'Color','b','LineWidth',2);
line([ind_S ind_S]*X.HEADER.DELTA,ylim,'Color','g','LineWidth',2);