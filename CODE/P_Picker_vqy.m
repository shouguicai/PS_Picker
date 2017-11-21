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
for code_ii = 1:length(code_index)
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
path = ['result/DAY817AP.csv'];
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
fid=fopen(path,'a+');
for ii = 1:length(tmp_index)
    [Loc_p, AAE]=P_cat(X_d,tmp_index(ii)); % AAE : average amplitude ratio 
    if Loc_p~=-1    
        time_submissionp  = datestr(t0+(Loc_p+8*3600)/86400,'yyyymmddHHMMSS.FFF'); % BEIJING Ê±¼ä
        fprintf(fid,'%s,%s,P,%s\n',station_name,time_submissionp(1:end-1),AAE);
    end
end
fclose(fid);
end