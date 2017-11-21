function [ind_S] = pickS_CF( data,fs,FB,T,n_smooth,first,last )
% XX_d,fs,[1 5;6 10;15 20],[1 2 3],10,1,length(XX_d)-1
[~,CF_matrix]=trace2FWkurto(data,fs,FB,T,n_smooth,first,last);
tmp = CF_matrix(:,1);
% tmp = sum(CF_matrix,2);
differ = diff(tmp);
[~,ind_S] = max(differ);
end

