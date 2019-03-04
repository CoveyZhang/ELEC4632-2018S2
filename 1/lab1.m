clear all
close all
clc

% %pre.1
% t = 0:0.1:100;
% t=t';
% 
% %pre.2
% y1 = sin(0.02 * pi .* t) + 0.4 * rand(length(t),1) - 0.2;
% y2 = cos(0.02 * pi .* t) + 0.4 * rand(length(t),1) - 0.2;
% 
% %pre.3
% t_new = t(1:800);
% y1_new = y1(201:1000);
% y2_new = y2(201:1000);
% 
% %pre.4
% data_new = [t_new y1_new y2_new];
% 
% %pre.5/6
% figure;
% subplot(2,1,1);
% plot(t,y1,'r',t,y2,'b');
% axis([0 140 -1.5 1.5]);
% grid on;
% legend('sin(0.02\pit)','cos(0.02\pit)');
% title('Original Data');
% xlabel({'Time(sec)','(a)'});
% ylabel('Data');
% 
% subplot(2,1,2);
% plot(t_new,y1_new,'g')
% hold on
% plot(t_new,y2_new,'color',[0.6 0.7 0.8]);
% axis([0 140 -1.5 1.5]);
% grid on;
% legend('sin(0.02\pit)','cos(0.02\pit)');
% title('Cut-off Data');
% xlabel({'Time(sec)','(b)'});
% ylabel('Data');

%1.a
load 'SysIdenData_StudentVersion.mat';

t = LogData.time;
y_act = LogData.signals(1).values(:,2);
y_actm = LogData.signals(1).values(:,1);
u_act = LogData.signals(2).values;

%1.b
figure;
subplot(2,1,1);
plot(t,y_act,'b',t,y_actm,'r');
axis([0 700 1 4]);
grid on;
legend('Noise-Reduced Output','Measured Output');
title('Actual Output Signal');
xlabel('Time(sec)');
ylabel('Water Level(V)');

subplot(2,1,2);
plot(t,u_act,'b');
axis([0 700 1 4]);
grid on;
legend('Actual Input');
title('Actual Input Signal');
xlabel('Time(sec)');
ylabel('Pump Voltage(V)');

%1.c
i = 2;
while u_act(i) == u_act(i-1)
    i = i + 1;
end
y_offset = mean(y_act(1:(i-1)));
u_offset = mean(u_act(1:(i-1)));
y = y_act - y_offset;
u = u_act - u_offset;

figure;
subplot(2,1,1);
plot(t,y,'r');
axis([0 700 -2 1]);
grid on;
legend('Actual Output');
title('Actual Offset-Free Output Signal');
xlabel('Time(sec)');
ylabel('Water Level(V)');

subplot(2,1,2);
plot(t,u,'b');
axis([0 700 -0.5 0.5]);
grid on;
legend('Actual Input');
title('Actual Offset-Free Input Signal');
xlabel('Time(sec)');
ylabel('Pump Voltage(V)');

%2.a
k = 3:(round(length(y)/2)+3);
Y = y(k);
phi = [];
for i=k
    temp=[y(i-1),y(i-2),u(i-1),u(i-2)];
    phi=[phi;temp];
end

%2.b
theta = ((phi'*phi)^-1)*phi'*Y;

%2.c
Ts = t(2)-t(1);
g=[0,1;theta(2),theta(1)];
h=[0;1];
c=[theta(4),theta(3)];
d=0;
sys=ss(g,h,c,d,Ts)
G=tf([theta(3),theta(4)],[1,-theta(1),-theta(2)],Ts)

%3
kk=(round(length(y)/2)+3):length(t);
Y1=lsim(sys,u(kk),t(kk));
Y2=lsim(sys,u,t);

figure;
subplot(2,1,1);
plot(t(1:length(kk)),Y1,'b--',t(1:length(kk)),y(kk),'r');
axis([0 350 -2 2]);
grid on;
legend('Simulated Output','Actual Output');
title('Offset-Free Model Vertification(2nd half)');
xlabel('Time(sec)');
ylabel('Water Level(V)');

subplot(2,1,2);
plot(t,Y2,'b--',t,y,'r');
axis([0 700 -2 2]);
grid on;
legend('Simulated Output','Actual Output');
title('Offset-Free Model Vertification(Entire)');
xlabel('Time(sec)');
ylabel('Pump Voltage(V)');

%post.1
phi1=[y(1:451) u(1:451)];
theta1=(phi1'*phi1)^(-1)*phi1'*y(2:452);
a1_new=-theta1(1);
b1_new=theta1(2);
G1=-a1_new;
H1=1;
C1=b1_new;
D1=0;
sys1=ss(G1,H1,C1,D1,Ts);

%post.2
Y3=lsim(sys1,u,t);
figure
plot(t,y,'r');
hold on
plot(t,Y2,'b--');
hold on
plot(t,Y3,'g:.');
hold on
axis([0 700 -1.5 1]);
grid on;
legend('Actual Output','1st Order Model Response','2nd Order Model Response');
title('Offset-Free Model Vertification(Entire)');
xlabel('Time(sec)');
ylabel('Pump Voltage(V)');
%post.3
MSE1=immse(y,Y3);
MSE2=immse(y,Y2);
MSE1=['MSE1 = ' num2str(MSE1)];
MSE2=['MSE2 = ' num2str(MSE2)];
text(0,0.6,MSE1);
text(0,0.4,MSE2);