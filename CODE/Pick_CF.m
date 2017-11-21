function [ ind_P, ind_S] = Pick_CF( data,fs,FB,T,n_smooth,first,last )
% XX_d,fs,[1 5;6 10;15 20],[1 2 3],10,1,length(XX_d)-1
[~,CF_matrix]=trace2FWkurto(data,fs,FB,T,n_smooth,first,last);
figure;
plot(CF_matrix);
tmp = CF_matrix(:,1);
% tmp = sum(CF_matrix,2);
differ = diff(tmp);
win = 200;
[~,ind_first] = max(differ);
differ_cut = differ(ind_first+win:end);
[~,ind_tmp] = max(differ_cut);
ind_second = ind_tmp + win +ind_first;
ind_P = min(ind_first,ind_second);
ind_S = max(ind_first,ind_second);
end

