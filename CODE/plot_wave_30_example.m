% clear all;
warning off all
%%
root = '../example30/example30/';
file_list = dir(root);
sel_num = 3;
% 15 30
imlist = file_list(2+sel_num).name
%%
file = [root,imlist];
X=rdsac(file);
%% calu arrival time for P and S wave
t_start = X.t(1);
t0 = X.t(1) - X.HEADER.B/86400;
p_arrival = t0 + X.HEADER.A/86400;
s_arrival = t0 + X.HEADER.T0/86400;
%%
fs = 1/X.HEADER.DELTA;  % Hz
x_t = (X.t-t0)*86400-X.HEADER.B;
x_t_p_arrival = (p_arrival-t0)*86400-X.HEADER.B;
x_t_s_arrival = (s_arrival-t0)*86400-X.HEADER.B;
figure(1)
subplot(3,1,1)
plot(x_t,X.d)
hold on
plot([x_t_p_arrival,x_t_p_arrival],[min(X.d),max(X.d)],'r',[x_t_s_arrival,x_t_s_arrival],[min(X.d),max(X.d)],'m','LineWidth',2.0)
hold off
xlim = [min(x_t),max(x_t)];
set(gca,'XLim',xlim)
grid on
% xlim = [min(X.t),max(X.t)];
% set(gca,'XLim',xlim)
% datetick('x','keeplimits');
% xlabel(sprintf('%s to %s',datestr(xlim(1)),datestr(xlim(2))))    
xlabel('Time(sec)')
ylabel('Count')
%%
F = 16;
[Y,F1,T1,P1] = spectrogram(X.d,F,F/2,F,fs,'yaxis');
%%
figure;
surf(T1,F1,(abs(P1)),'EdgeColor','none');
axis xy; axis tight; 
colormap(jet); 
view(0,90);
colorbar
% caxis([-120,-50]);
xlabel('Time (sec)','FontSize',12);
ylabel('Frequency (Hz)','FontSize',12);
%%
% % figure;
% % surf(T1,F1,(abs(P1)),'EdgeColor','none');
% % axis xy; axis tight; 
% % colormap(jet); 
% % view(0,0);
% % colorbar
% % % caxis([-120,-50]);
% % xlabel('Time (sec)','FontSize',12);
% % % %%
% % % figure;
% % % plot(T1,max(P1))
%%
z = hilbert(X.d);
amp = abs(z);
figure(1)
subplot(3,1,2)
plot(x_t,amp)
hold on
plot([x_t_p_arrival,x_t_p_arrival],[min(amp),max(amp)],'r',[x_t_s_arrival,x_t_s_arrival],[min(amp),max(amp)],'m','LineWidth',2.0)
hold off
xlim = [min(x_t),max(x_t)];
set(gca,'XLim',xlim)
xlabel('Time(sec)')
ylabel('Amp')
grid on
title('Instantaneous Amp')
instfreq = fs/(2*pi)*diff(unwrap(angle(z)));
figure(1)
subplot(3,1,3)
plot(x_t(2:end),instfreq)
hold on
plot([x_t_p_arrival,x_t_p_arrival],[min(instfreq),max(instfreq)],'r',[x_t_s_arrival,x_t_s_arrival],[min(instfreq),max(instfreq)],'m','LineWidth',2.0)
hold off
xlabel('Time(sec)')
ylabel('Hz')
xlim = [min(x_t),max(x_t)];
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
XX_d = X.d;
% XX_d = filter(Hfilter,XX_d);
[loc, snr_db] = PphasePicker(XX_d, X.HEADER.DELTA, type, pflag, Tn, xi, nbins, o);
ind = aic_pick(XX_d,'whole');
[Loc_p, ~]=P_cat(XX_d,2000); % AAE : average amplitude ratio 
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
    loc_S_1 = pickS(loc,E,out,X.HEADER.DELTA);
%     loc_S = s_catch(XX_d);
%     line([loc_S loc_S]*X.HEADER.DELTA,ylim,'Color','k','LineWidth',1);
    line([loc_S_1 loc_S_1]*X.HEADER.DELTA,ylim,'Color','b','LineWidth',1);    % pickS
    line([x_t_s_arrival x_t_s_arrival],ylim,'Color','m','LineWidth',1);
    hold off
    
    figure(1)
    subplot(3,1,1)
%     line([loc_S loc_S]*X.HEADER.DELTA,ylim,'Color','k','LineWidth',2);
    hold on;
    line([loc_S_1 loc_S_1]*X.HEADER.DELTA,ylim,'Color',[0,1,1],'LineWidth',2);
    line([loc loc],ylim,'Color','g','LineWidth',2);                  % phase picker
    line([ind ind]*X.HEADER.DELTA,ylim,'Color','b','LineWidth',2);   % aic
    line([Loc_p Loc_p],ylim,'Color',[0,0.5,0.5],'LineWidth',2);      % imporved aic
%%
% % max_p = max(P1.');
% % norm_P = P1./(repmat(max_p.',1,length(T1)));
% % % figure;
% % % surf(T1,F1,(abs(norm_P)),'EdgeColor','none');
% % % axis xy; axis tight;
% % % colormap(jet);
% % % view(0,90);
% % %%
% % tmp = norm_P(1:3,:);
% % sum_tmp = sum(tmp,1);
% % figure;
% % plot(T1,sum_tmp);hold on
% % % plot(T1(2:end),diff(sum_tmp),'r');
% % sta = 1;
% % loc_s = pickS(sta,sum_tmp,sta:length(sum_tmp),1);
% % line([T1(loc_s),T1(loc_s)],[ylim],'Color','y','LineWidth',2);
% % figure(1)
% % subplot(3,1,1)
% % line([T1(loc_s),T1(loc_s)],[ylim],'Color','y','LineWidth',2);
%%
max_p = max(P1);
% max_p = P1(1,:);
figure(9);
plot(T1,max_p);hold on
% line(xlim,4*[mean(max_p),mean(max_p)],'Color','g','LineWidth',2);
% sta = 1;
% loc_s = pickS(sta,max_p,sta:length(max_p),1);
delta_T = T1(2) - T1(1);
min_pos = ceil(Loc_p/delta_T);
[~,max_pos] = max(max_p);
% therold_max_pos = min_pos + ceil(5/delta_T);
% if max_pos < therold_max_pos
%     [~,max_pos] = max(max_p(therold_max_pos:end));
%     max_pos = max_pos + therold_max_pos;
% end
if max_pos <= min_pos
    max_pos = max_pos + ceil(5/delta_T);
end
sta = min_pos;
loc_s = aic_pick(max_p(sta:max_pos),'whole')+sta;
line([T1(loc_s),T1(loc_s)],ylim,'Color','y','LineWidth',2);
figure(1)
subplot(3,1,1)
line([T1(loc_s),T1(loc_s)],ylim,'Color','y','LineWidth',2);
%%
% alarm = envelop_hilbert_v2(max_p,round(fs),1,round(fs),1);
%% Test for Kurtosis©\Based method
% sta1 = floor((Loc_p)/X.HEADER.DELTA);
% data = XX_d(sta1:end);
% [ind_S] = pickS_CF(data,fs,[1 5;6 10;11 15;16 20],[0.5 1 1.5 2 3],10,1,length(data)-1);
% figure(66);
% ind_S = ind_S + sta1;
% line([ind_S ind_S]*X.HEADER.DELTA,ylim,'Color','g','LineWidth',2);

FB = [1 5;6 10;11 15;15 20];
T = [1 2 3];
n_smooth = 10;
start_interst = 1;
end_interst = floor(3*length(XX_d)/4);
[ind_P, ind_S] = Pick_CF(XX_d,fs,FB,T,n_smooth,start_interst,end_interst);
figure(66);
line([ind_P ind_P]*X.HEADER.DELTA,ylim,'Color','b','LineWidth',2);
line([ind_S ind_S]*X.HEADER.DELTA,ylim,'Color','g','LineWidth',2);