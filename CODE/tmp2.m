% amp = abs(hilbert(XX_d));
% Smooth_window = 30;
% env = conv(amp,ones(1,Smooth_window)/Smooth_window);
% env1 = conv(env,ones(1,Smooth_window)/Smooth_window);
% plot(x_t,XX_d);
% hold on;
% plot(x_t,env1(Smooth_window:length(amp)+Smooth_window-1),'r-.');


tmp = max(P1);