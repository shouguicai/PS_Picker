function [ loc_S ] = pickS(loc,E,out,delta)
%PICKSS Summary of this function goes here
%   Detailed explanation goes here
loc_P = loc/delta;
E_out = E(1,out);
[~,loc_MAX] = max(E_out(loc_P:end));
loc_MAX = loc_MAX + loc_P;
search_start =floor((loc_P+loc_MAX)/2);
% search_start =loc_P;%+1/delta;
search_end = loc_MAX;
E_search = E_out(search_start:search_end);
N_search = length(E_search);
search_ii = floor(N_search/2);
while(max(E_search(1:search_ii-1))>=min(E_search(search_ii+1:end)))
    search_ii = floor((search_ii+N_search)/2);
    if search_ii-floor((search_ii+N_search)/2) == 0
        break;
    end
end

search_jj = search_ii;

% while(E_search(search_jj)>mean(E_search(1:search_jj-1)))
%     search_jj = search_jj - 1;
% end

loc_S = search_start + search_jj;

end

