clear all;close all;
root = '../example30/example30/';
file_list = dir(root);
%%   找到Z结尾的文件的index
N_file = length(file_list)-2;
code_index = [];
for ii=1:N_file  
    imlist = file_list(2+ii).name;
    if imlist(end-4) == 'Z'
        code_index =[code_index,ii];
    end
end
%%
% code_skip_after = [17,90,91,151,152,157,164,191,239,249,330,359:376,432:436];
% code_skip_before = [58,187,217];
for code_ii =1:length(code_index) 
    sel_num = code_index(code_ii);
imlist = file_list(2+sel_num).name;
%% 读取sac 文件
file = [root,'\',imlist];
X=rdsac(file);
X_t = X.t;
% X_d = detrend(X.d);
X_d = X.d;
station_name = X.HEADER.KSTNM;
path = ['result/DAY822.csv'];
flag_cre_file = 0;

%%
t_start = X_t(1);
t0 = X_t(1) - X.HEADER.B/86400;
x_t = (X_t-t0)*86400-X.HEADER.B;

WW=2000;
%  fid=fopen(path,'a+');
for kkk=1:length(WW)

[Loc_p, AAE]=P_cat(X_d,WW(kkk)); % AAE : average amplitude ratio 
if Loc_p~=-1    
     figure(66);plot(x_t,X_d);hold on;
     line([Loc_p Loc_p],[ylim],'Color','k','LineWidth',2);
     pp=X.HEADER.A-X.HEADER.B;
     line([pp pp],[ylim],'Color','r','LineWidth',2);
%       time_submissionp  = datestr(t0+(Loc_p+8*3600)/86400,'yyyymmddHHMMSS.FFF'); % BEIJING 时间
%       fprintf(fid,'%s,%s,P,%s\n',station_name,time_submissionp(1:end-1),AAE);
end
end
end




