clear
close
clc

load 'SysIdenData_1.mat';

t = LogData.time(1066:end);
y_act = LogData.signals(1).values(1066:end,2);
y_actm = LogData.signals(1).values(1066:end,1);
u_act = LogData.signals(2).values(1066:end);
t = t-t(1);

i = 2;
while u_act(i) == u_act(i-1)
    i = i + 1;
end
y_offset = mean(y_act(1:(i-1)));
u_offset = u_act(1);
y = y_act - y_offset;
u = u_act - u_offset;

k = 3:(round(length(y)/2)+3);
Y = y(k);
phi = [];
for i=k
    temp=[y(i-1),y(i-2),u(i-1),u(i-2)];
    phi=[phi;temp];
end

theta = ((phi'*phi)^-1)*phi'*Y;

Ts = t(2)-t(1);
g=[0,1;theta(2),theta(1)];
h=[0;1];
c=[theta(4),theta(3)];
d=0;
sys=ss(g,h,c,d,Ts)

Kp=0.65;
Ki=0.025;
sim('lab5_simulink.slx');

y2 = ScopeData(:,2) ;
u2 = ScopeData(:,3) ;

load('SFControlData_1.mat')
treal = SFLogData.time;
yref = SFLogData.signals(1).values(:,1);
yreal = SFLogData.signals(1).values(:,2); 
ureal = SFLogData.signals(2).values;
y2 = y2(1:701)+y_offset ;
figure(1);
subplot(211);
plot(treal, y2 ,'b');
hold on;
plot(treal, yreal,'r');
plot(treal, yref,'g');
xlim([0, 600]);
ylim([-1, 5]);
xlabel ('Time(sec)');
ylabel ('Water Level(V)');
title('Output Feedback Control Result: Output Signal');
grid on;
legend ('Simulated Output','Actual Output','Reference Output');
u2 = u2 + u_offset;
subplot(212);
plot(treal,ureal,'r')
hold on;
plot(treal,u2(1:701),'b')
xlim([0, 600]);
ylim([1, 3]);
plot ([0 550],[2.5 2.5],'g:','linewidth',2);
plot ([0 550],[1.5 1.5],'g:','linewidth',2);
UMIN='V_m_i_n = 1.5';
UMAX='V_m_a_x = 2.5';
text(100,2.6,UMAX);
text(100,1.4,UMIN);
xlabel ('Time(sec)');
ylabel ('Pump Voltage(V)');
title('Control Input Signal');
grid on;
legend ('Simulated Control Input','Actual Control Input');