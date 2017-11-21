% clear all;close all;clc;


% at first we must determine a rough pisition where an earthquake happened.
% then we can according to the pisition and cut a period of signal (about
% 40s).The picked signal as an input of this function .

%%%% read file .sac
% N1=0;N2=0;
% Th=40;
% root = 'F:\余震检测\workspace\example30';
% file_list = dir(root);
% N_file = length(file_list)-2;
% code_index = [];
% for ii=1:N_file  
%     imlist = file_list(2+ii).name;
%     if imlist(end-4) == 'Z'
%         code_index = [code_index,ii];
%     end
% end
% for jj=13%:length(code_index)
% file=[root,'\',file_list(code_index(jj)+2).name];
% x=XZ.d;
% figure ;plot(x);
% y=x(2064000:2068000);
% figure ;plot(y);
% file='XX.PWU.2008185000000.BHZ';
% X=rdsac(file);
% X_t = X.t;
% X_d = detrend(X.d);
% X_d = X.d;
% load WW.mat
function [Loc_p, ER]=P_cat(X_d,WW)

% INPUT :
%        WW  : a rough position about P/S
%        X_d : raw data 
 start = max(1,floor(WW-20/0.01));
 ends = min(start + 6000,length(X_d));
 x=X_d(start:ends);
 clear XX
 v=1;
% for j=1:length(x)
%     if (x(j))>-1000&& (x(j))<100;
%         XX(v)=x(j);v=v+1;
%     end
% end
x=x-mean(x);% 消除偏置 使STA/LTA比值更合理
%%%% ESM algorithm
signal=hilbert(x);
CF=abs(signal);% characteristic function
%%%% sta/lta algorithm 
L1=40;% length of short window 
L2=500; % length of long window ,usually L2=(5~10)*L1;
for ii=1:length(CF)-L2-1
   RR(ii)=mean(CF(ii+L2-L1:ii+L2-1))/mean(CF(ii:ii+L2-1));
end
R=[zeros(1,L2),RR];
% figure;plot(R);
SR=conv(ones(1,40),R)/40;
% hold on;plot(SR(11:end),'r');

Q=SR(11:end);%% after a smooth operation ,there a approimately 10 point offset ,so correct here
% equal_id=1;
% while(exist('equal_id','var') )
% clear equal_id
% q=1;
% for j=1:length(Q)-1
%     if Q(j)==Q(j+1)
%         equal_id(jj)=j;
%         Q(j)=Q(j)-0.001;
%         q=q+1;
%     end
% end
% end
% 


for j=2:length(Q);% do not consider the first point 
    if (Q(j)-Q(j-1)>0)&&(Q(j+1)-Q(j))<=0
        aa=j;
    break;
    end
end

cnt=1;
% K=aa+5000;
for j=2:length(Q)-1;% do not consider the first point 
    if (Q(j)-Q(j-1)>0)&&(Q(j+1)-Q(j))<=0
        extreme_index(cnt)=j;
        extreme_value(cnt)=Q(j);
        cnt =cnt+1;
    end
end
% figure;plot(Q)%;hold on;plot(extreme_index,extreme_value,'*');
% 
%% new method 
T_p=2.4;
count=1;
for i=1:length(Q)
    if Q(i)>T_p
        index1=i;
        break;
    end
end

if exist('index1','var')
    for i=index1:length(Q)
        if Q(i)<T_p
            index2=i;
            break;
        end
    end
end

if exist('index1','var')&&exist('index2','var')
% figure ;plot(x);
    [ER,IDP1]=max(Q(index1:index2));   
    IDP=IDP1+index1-1;
    P_R1=max(IDP-200,1);
    P_R2=min(IDP+20,6000);
    Sig_P=R(P_R1:P_R2);
    indp = aic_pick(Sig_P,'whole')+P_R1;
    Loc_p=(indp+start-1)*0.01;
%     figure(66);plot(x);
%     line([indp indp],[ylim],'Color','k','LineWidth',2);
end

if ~exist('index1','var')
    Loc_p=-1;
    ER=0;
end





% T_s=2.4;
% count=1;
% for i=1:length(Q)
%     if Q(i)>T_s
%         index3(count)=i;count=count+1;
%         if count>
%         break;
%     end
% end
% 
% if exist('index3','var')
%     for i=index3:length(Q)
%         if Q(i)<T_s
%             index4=i;
%             break;
%         end
%     end
% end
% 
% if exist('index4','var')&&exist('index3','var')
% % figure ;plot(x);
%     [ER,IDP1]=max(Q(index3:index4));   
%     IDP=IDP1+index1-1;
%     P_R1=max(IDP-200,1);
%     P_R2=min(IDP+20,6000);
%     Sig_P=R(P_R1:P_R2);
%     indp = aic_pick(Sig_P,'whole')+P_R1;
%     Loc_p=(indp+start-1)*0.01;
% %     figure;plot(x);
% %    line([indp indp],[ylim],'Color','k','LineWidth',2);
% end
% 
% if ~exist('index1','var')
%     Loc_p=-1;
%     ER=0;
% end

end

        
  


% %
% EX.index=extreme_index;EX.value=extreme_value;
% [EXnew.value,index] = sort([EX(:).value],'descend');
% EXnew.index=EX.index(index);
% %%%%% avoid the interference of some extreme point around max ,delete them
% % A=sort(EXnew.index(1:3));
% B1=EXnew.index(1);
% B2=EXnew.index(2);
% % C=A(end);
% cnt1=1;
% for k=1:length(EXnew.index)
%     if((EXnew.index(k)-B1)<=100&&(EXnew.index(k)-B1)>-100)||((EXnew.index(k)-B2)<=100&&(EXnew.index(k)-B2)) % 100 need adjust
%         if(EXnew.index(k)~=B1 &&EXnew.index(k)~=B2)
%             pesudo_index(cnt1)=k;cnt1=cnt1+1;
%         end
%     end
% end
% 
% 
% XX=EXnew.index;
% YY=EXnew.value;
% if exist( 'pesudo_index', 'var')
% XX(pesudo_index)=[];
% YY(pesudo_index)=[];
% end
% % 
% % figure;plot(Q);hold on;plot(XX,YY,'*');hold on;plot(x)
% 
% % D=sort(YY(1:3));
% % IDP=D(1);IDS=D(2); %%%% 取第三大是有问题的
% % if  abs>1200; 
% %     IDS=D(2);
% % end
% 
% MM=sort(XX);
% M=MM(1);
% if XX(1)<XX(2)
%     IDP=XX(1);
%     SV1=YY(2);%SV2=YY(3);
%     IDS1=XX(2);
%     TT=find(IDS1==extreme_index)-1;
%     if TT>0;
%         IDS2=extreme_index(TT);
%         SV2=extreme_value(TT);
%     else
%         IDS2=extreme_index(TT+1);
%         SV2=extreme_value(TT+1);
%     end
%     RatioS=1.3;
%     if SV1/SV2>RatioS
%         IDS=IDS1;
%     else if IDS1<IDS2
%             IDS=IDS1;
%         else
%             IDS=IDS2;
%         end
%     end
% else
%     IDS=XX(1);
%     if (XX(3)<XX(2))&&(XX(3)~=M)
%         IDP=XX(3);
%     else
%         IDP=XX(2);
%     end
% %%%%%%%%  useless 
% % % % % % % % %;SV1=YY(2);        
% % % % % % % % %     TT=find(IDP1==extreme_index)-1;
% % % % % % % % % if TT>0;
% % % % % % % % % IDP2=extreme_index(TT);
% % % % % % % % % SV2=extreme_value(TT);
% % % % % % % % % else
% % % % % % % % %  IDP2=extreme_index(TT+1);
% % % % % % % % %  SV2=extreme_value(TT+1);
% % % % % % % % % end
% % % % % % % % % 
% % % % % % % % % RatioP=1.5;% maybe need adjust 
% % % % % % % % % if SV1/SV2>RatioP
% % % % % % % % %     IDP=IDP1;
% % % % % % % % % else if IDP1<IDP2
% % % % % % % % %         IDP=IDP1;
% % % % % % % % %     else
% % % % % % % % %         IDP=IDP2;
% % % % % % % % %     end
% % % % % % % % % end
% end  
% %%%% determine the search region  [IDP-120 IDP+100]  [IDS-120 IDS+100]
% %%%% adjust 
% P_R1=IDP-120;
% P_R2=IDP+20;
% S_R_1=IDS-120;
% S_R_2=IDS+20;
% Sig_P=R(P_R1:P_R2);
% Sig_S=R(S_R_1:S_R_2);
% ss=XZ.HEADER.T0-XZ.HEADER.B;
% pp=XZ.HEADER.A-XZ.HEADER.B;
% % line([ss ss]/XZ.HEADER.DELTA,[ylim],'Color','b','LineWidth',2);
% line([pp pp]/XZ.HEADER.DELTA,[ylim],'Color','r','LineWidth',2);


% % ps=R(920:1150);
% % s1=R(1250:1450);

% inds = aic_pick(Sig_S,'whole')+S_R_1;

% line([inds inds],[ylim],'Color','g','LineWidth',2);

% % if abs(indp-pp/0.01)<=Th
% %     N1=N1+1;
% %     AC1(N1)=abs(indp*0.01-pp);
% % end
% % if abs(inds-ss/0.01)<=Th
% %     N2=N2+1;
% %     AC2(N2)=abs(indp*0.01-pp);
% % end
% % 
% clear  XX YY EXnew  EX extreme_index   extreme_value  pesudo_index ;

% end




% hitrate=(length(AC1)+length(AC2))/56;
% RMSE=sqrt((sum((AC1).^2)+sum((AC2).^2)+58-(length(AC1)+length(AC2)))/56);