% [mean_M,CF_matrix]=trace2FWkurto(XX_d,fs,[1 5;6 10;15 20],[1 2 3],10,1,length(XX_d)-1);
[mean_M,CF_matrix]=trace2FWkurto(XX_d,fs,[1 5;6 10;11 15;15 20],[1 2 3],10,1,length(XX_d)-1);
[ind_pick,vals_kurto]=follow_extrem2(CF_matrix,10,1);
%%
plot(CF_matrix);
%%
% tmp = CF_matrix(:,1);
tmp = sum(CF_matrix,2);
differ = diff(tmp);
% [PKS,LOCS]= findpeaks(differ,'SORTSTR','descend'); 
win = 200;
[~,ind_P] = max(differ);
differ_cut = differ(ind_P+win:end);
[~,ind_tmp] = max(differ_cut);
ind_S = ind_tmp + win +ind_P;
figure;
plot(differ);hold on;
plot(tmp,'r');
line([ind_P ind_P],ylim,'Color','r','LineWidth',2);
line([ind_S ind_S],ylim,'Color','m','LineWidth',2);