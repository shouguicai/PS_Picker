function loc_s=s_catch(z)

% clear all ;clc;close all;
% file='01.JMG.BHZ.SAC';
% file1=file;file1(end-4)='E';
% file2=file;file2(end-4)='N';
% XE=rdsac(file1);
% XN=rdsac(file2);
% XZ=rdsac(file);
% z=XZ.d;

for i=1:length(z)-1
    Q(i)=((z(i+1)-z(i))^2)/0.01;%%(x(i+1)-x(i))^2+(y(i+1)-y(i))^2+
end
% figure;plot(Q);
N1=51;% length of short window   
N2=200;% length of long  window 
for ii=1:length(Q)-N2-N1-1
   R(ii)=mean(Q(ii:ii+N1))/mean(Q(ii+N1+1:ii+N1+N2+1));
end
% figure;plot(R);
[~,maxId]=max(R);
loc_s=(maxId);

figure(66)
hold on
plot([1:length(R)]/100,R*2e3,'r-.');hold on
% line([maxId maxId]*X.HEADER.DELTA,ylim,'Color','y','LineWidth',1);hold on
% 
% search_start =1;
% search_end = maxId;
% E_search = R(search_start:search_end);
% N_search = length(E_search);
% search_ii = floor(N_search/2);
% while(max(E_search(1:search_ii-1))>=min(E_search(search_ii+1:end)))
%     search_ii = floor((search_ii+N_search)/2);
%     if search_ii-floor((search_ii+N_search)/2) == 0
%         break;
%     end
% end
% search_jj = search_ii;
% loc_s = search_start + search_jj;

end

%%



    
