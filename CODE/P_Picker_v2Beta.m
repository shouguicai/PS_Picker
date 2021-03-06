clear all;
root = '../after/';
file_list = dir(root);
%%
N_file = length(file_list)-2;
code_index = [];
for ii=1:N_file
    imlist = file_list(2+ii).name;
    if imlist(end) == 'Z'
        code_index = [code_index,ii];
    end
end
%%
% code_skip_after = [17,90,91,151,152,157,164,191,239,249,330,359:376,432:436];
% code_skip_before = [58,187,217];
for code_ii = 19:length(code_index)
% if ~isempty(find(code_ii ==code_skip_after,1)) 
%     continue;
% end
sel_num = code_index(code_ii);
imlist = file_list(2+sel_num).name;
%%
file = [root,imlist];
X=rdsac(file);
X_t = X.t;
% X_d = detrend(X.d);
X_d = X.d;
%%
station_name = X.HEADER.KSTNM;
% path = ['result/',station_name,'.csv'];
path = ['result/DAY818APS.csv'];
% path = ['result/',imlist,'.csv'];
flag_cre_file = 0;

%%
t_start = X_t(1);
t0 = X_t(1) - X.HEADER.B/86400;
x_t = (X_t-t0)*86400-X.HEADER.B;

% % figure(1)
% % subplot(3,1,1)
% % % 
% % plot(x_t,X_d)
% % xlim = [min(x_t),max(x_t)];
% % set(gca,'XLim',xlim)
% % grid on
% % xlabel('Time(sec)')
% 
% plot(X_t,X_d)
% xlim = [min(X_t),max(X_t)];
% set(gca,'XLim',xlim)
% datetick('x','keeplimits');
% xlabel(sprintf('%s to %s',datestr(xlim(1)),datestr(xlim(2))))
% 
% ylabel('Count');
%%
fs = 1/X.HEADER.DELTA;  % Hz
%%
z = hilbert(X_d);
amp = abs(z);

threshold1 = 5*mean(amp); % need to be adaptive
trigger = 1e5*(sign(amp-threshold1)+1);
trigger_index = find(trigger>0);  
if size(trigger_index,1) == 0
    continue ;
end
tmp = trigger;
differ = diff(trigger_index);
tmp(trigger_index(1)) = 1;
threshold2 = 5e3;  % signal palse: 50s ,fixed
tmp(trigger_index(2:end)) = sign(differ-threshold2)+1;
tmp_index = find(tmp>0); 

%% Run PphasePicker with optional picking parameters
type = 'na'; % no bandpass filtering   na sm 
pflag = 'Y'; % To plot waveform and P-phase onset

% USE Tn = 0.1 for records with low sampling rate lower than 100 sps 
Tn = 0.1;       % undamped natural period of oscillator in second 
xi = 0.6;       % damping ratio 
nbins = 200;    % histogram bin size
o = 'full';   % 'to_peak' to take segment of waveform from beginning to
                % absolute peak value
Hfilter = noisefilter;
%%
for ii = 1:length(tmp_index)
    start = max(1,floor(tmp_index(ii)-45/X.HEADER.DELTA));
    ends = min(start + 4500*2,length(X_d));
    XX_d = X_d(start:ends);
%     XX_d = filter(Hfilter,XX_d);
    if max(abs(XX_d))==0
        continue;
    end

    [loc, snr_db] = PphasePicker(XX_d.', X.HEADER.DELTA, type, pflag, Tn, xi, nbins, o);
    
%     if ind ~= 0
%         line([ind ind]*X.HEADER.DELTA,[ylim],'Color','b','LineWidth',2);
%     end
    
    if loc ~= -1.00
%         if flag_cre_file ~=1
%             fid=fopen(path,'a+');
%             flag_cre_file = 0;
%         end
        ind = aic_pick(XX_d,'whole');
        line([ind ind]*X.HEADER.DELTA,ylim,'Color','b','LineWidth',2);
        % STA
        % define the window
        wintype = 'rectwin';
        winlen = 51;
        winamp = [0.5,1]*(1/winlen);
        
        %*** Z ****%
        % find the zero-crossing rate
% %         E = energy(XX_d.',wintype,winamp(2),winlen);
% %         % time index for the ST-ZCR and STE after delay compensation
% %         N = length(XX_d); % signal length
% %         out = (winlen-1)/2:(N+winlen-1)-(winlen-1)/2;
% %         loc_S_Z = pickS(loc,E,out,X.HEADER.DELTA);
            F = 16;
            [~,~,T1,P1] = spectrogram(XX_d,F,F/2,F,fs,'yaxis');
            max_p = max(P1);
            delta_T = T1(2) - T1(1);
            min_pos = min(ceil(ind*X.HEADER.DELTA/delta_T),length(max_p)-10);
            [~,max_pos] = max(max_p);
            if max_pos <= min_pos
                max_pos = min(min_pos + ceil(5/delta_T),length(max_p));
            end
            sta = min_pos;
            loc_s = aic_pick(max_p(sta:max_pos),'whole')+sta-1;
            loc_S_Z = T1(loc_s);

        %*** N ****%
        sel_num_N = sel_num -1;
        imlist_N = file_list(2+sel_num_N).name;
        file_N = [root,imlist_N];
        X_N = rdsac(file_N);
        X_d_N = X_N.d;
        XX_d_N = X_d_N(min(start,length(X_d_N)):min(ends,length(X_d_N)));
% %         % find the zero-crossing rate
% %         E_N = energy(XX_d_N.',wintype,winamp(2),winlen);
% %         % time index for the ST-ZCR and STE after delay compensation
% %         N = length(XX_d_N); % signal length
% %         out = (winlen-1)/2:(N+winlen-1)-(winlen-1)/2;
% %         loc_S_N = pickS(loc,E_N,out,X.HEADER.DELTA);
            F = 16;
            [~,~,T1,P1] = spectrogram(XX_d_N,F,F/2,F,fs,'yaxis');
            max_p = max(P1);
            delta_T = T1(2) - T1(1);
            min_pos = min(ceil(ind*X.HEADER.DELTA/delta_T),length(max_p)-10);
            [~,max_pos] = max(max_p);
            if max_pos <= min_pos
                max_pos = min(min_pos + ceil(5/delta_T),length(max_p));
            end
            sta = min_pos;
            loc_s = aic_pick(max_p(sta:max_pos),'whole')+sta-1;
            loc_S_N = T1(loc_s);
            
        
        %*** E ****%
        sel_num_E = sel_num -2;
        imlist_E = file_list(2+sel_num_E).name;
        file_E = [root,imlist_E];
        X_E=rdsac(file_E);
        X_d_E = X_E.d;
        XX_d_E = X_d_E(min(start,length(X_d_E)):min(ends,length(X_d_N)));
% %         % find the zero-crossing rate
% %         E_E = energy(XX_d_E.',wintype,winamp(2),winlen);
% %         % time index for the ST-ZCR and STE after delay compensation
% %         N = length(XX_d_E); % signal length
% %         out = (winlen-1)/2:(N+winlen-1)-(winlen-1)/2;
% % %         t = (out-(winlen-1)/2)*(1/fs);
% %         loc_S_E = pickS(loc,E_E,out,X.HEADER.DELTA);
            F = 16;
            [Y,F1,T1,P1] = spectrogram(XX_d_E,F,F/2,F,fs,'yaxis');
            max_p = max(P1);
            delta_T = T1(2) - T1(1);
            min_pos = min(ceil(ind*X.HEADER.DELTA/delta_T),length(max_p)-10);
            [~,max_pos] = max(max_p);
            if max_pos <= min_pos
                max_pos = min(min_pos + ceil(5/delta_T),length(max_p));
            end
            sta = min_pos;
            loc_s = aic_pick(max_p(sta:max_pos),'whole')+sta-1;
            loc_S_E = T1(loc_s);
            
        
        loc_S_ENZ = [loc_S_E,loc_S_N,loc_S_Z];
        loc_S_ENZ_SORT = sort(loc_S_ENZ);
        if (isempty(loc_S_E)+isempty(loc_S_N)+isempty(loc_S_Z))>0
                loc_S = loc_S_ENZ_SORT(1);
        else
                loc_S = loc_S_ENZ_SORT(2);
        end

%         figure(66);
%         hold on;
%         plot(t,max(XX_d)*E(out)/max(E(out)),'b-.','Linewidth',2); xlabel('t, seconds');
        line([loc_S loc_S],ylim,'Color','m','LineWidth',2);
%         hold off;
        
%          Loc_shift_secs = start*X.HEADER.DELTA+loc; % sec
%          time_submissionp  = datestr(t0+(Loc_shift_secs+8*3600)/86400,'yyyymmddHHMMSS.FFF'); % BEIJING 时间

%         fprintf(fid,'%s,%s,P,%s\n',station_name,time_submissionp(1:end-1),snr_db);
        time_submissions  = datestr(t0+((start)*X.HEADER.DELTA+loc_S+8*3600)/86400,'yyyymmddHHMMSS.FFF'); % BEIJING 时间
        time_submissionp_aic  = datestr(t0+((start+ind)*X.HEADER.DELTA+8*3600)/86400,'yyyymmddHHMMSS.FFF'); % BEIJING 时间
        fid=fopen(path,'a+');
        fprintf(fid,'%s,%s,S,%s\n',station_name,time_submissions(1:end-1),snr_db);
        fprintf(fid,'%s,%s,P,%s\n',station_name,time_submissionp_aic(1:end-1),snr_db);      
        fclose(fid);
    end
end

% fclose(fid);
end