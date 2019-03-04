clear all
close all
clc

%1.a
load 'SysIdenData_1.mat';

t = LogData.time(1066:end);
y_act = LogData.signals(1).values(1066:end,2);
y_actm = LogData.signals(1).values(1066:end,1);
u_act = LogData.signals(2).values(1066:end);
t = t-t(1);

%1.b
figure;
subplot(2,1,1);
plot(t,y_act,'b',t,y_actm,'r');
axis([0 750 1 4]);
grid on;
legend('Noise-Reduced Output','Measured Output');
title('Actual Output Signal');
xlabel('Time(sec)');
ylabel('Water Level(V)');

subplot(2,1,2);
plot(t,u_act,'b');
axis([0 750 1 4]);
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
axis([0 750 -2.5 1.5]);
grid on;
legend('Actual Output');
title('Actual Offset-Free Output Signal');
xlabel('Time(sec)');
ylabel('Water Level(V)');

subplot(2,1,2);
plot(t,u,'b');
axis([0 750 -0.5 0.5]);
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
axis([0 400 -2 2]);
grid on;
legend('Simulated Output','Actual Output');
title('Offset-Free Model Vertification(2nd half)');
xlabel('Time(sec)');
ylabel('Water Level(V)');

subplot(2,1,2);
plot(t,Y2,'b--',t,y,'r');
axis([0 750 -2 2]);
grid on;
legend('Simulated Output','Actual Output');
title('Offset-Free Model Vertification(Entire)');
xlabel('Time(sec)');
ylabel('Pump Voltage(V)');