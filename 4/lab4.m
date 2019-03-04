clear
close
clc

load 'SysIdenData_4.mat';
UMIN=['u_m_i_n = ' num2str(0.55)];
UMAX=['u_m_a_x = ' num2str(0.55)];
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

gg=[0,1;theta(2),theta(1)]';
hh=[0;1]';
cc=[theta(4),theta(3)]';
dd=0;
sys_o=ss(gg,cc,hh,dd,Ts)

%1.a
x0 = [0; 0];

%1.b
L_db = acker(gg,cc,[0 0])
L_ndb = acker(gg,cc,[0.4  0.9])

%1.c
%Deadbeat
sysdb=ss((gg-cc*L_db),cc,hh,dd,Ts);
[Y1,T1,X1]=initial(sysdb,x0,0:Ts:50);
U1=-L_db*X1';
%non
sysndb=ss((gg-cc*L_ndb),cc,hh,dd,Ts);
[Y2,T2,X2]=initial(sysndb,x0,0:Ts:50);
U2=-L_ndb*X2';

y_ref=[0 0.7 0.2 -0.5 0];
gain_db=dcgain(sysdb);
gain_ndb=dcgain(sysndb);
Period=140;
one=ones(1,Period);
temp=[];
Y_ref=[];
for i=1:5
    temp=one.*y_ref(i);
    Y_ref=[Y_ref,temp];
end
len=length(Y_ref);
L_ndb1=acker(gg,cc,[0.92 0.92]);
sys_ndb1=ss(gg-cc*L_ndb1,cc,hh,dd,Ts);
gain_ndb1=dcgain(sys_ndb1);

[y_ndb1,T,x_ndb1]=lsim(sys_ndb1,Y_ref*1/gain_ndb1);
u_ndb1=-L_ndb1*x_ndb1'+1/gain_ndb1*Y_ref;

y_sp=zeros(1,len);
x_sp=zeros(2,len);
u_sp=zeros(1,len);
for k=(1:len)
    u_sp(k)=Y_ref(k)/gain_ndb1-L_ndb1*x_sp(:,k);
    y_sp(k)=hh*x_sp(:,k);
    if k ~=len
        x_sp(:,k+1)=gg*x_sp(:,k)+cc*u_sp(k);
    end
end
figure;
subplot(2,1,1);
stairs(T, y_sp,'r');
hold on;stairs( T, Y_ref,'g','linewidth',2);
hold off;
set(gca,'YTick',[-1 -0.2 0 0.5 0.7]);
ylim([-1,1]);grid on;
title({'Set-Point Control Results: Simulation','Output Signal'});
legend('Simulated Output','Reference Output');
xlabel({'Time(sec)','(a)'});ylabel({'Offset-Free','Water Level(V)'});
subplot(2,1,2);
stairs(T,u_sp,'b');
hold on;
plot ([0 525],[0.5 0.5],'g:','linewidth',2);
plot ([0 525],[-0.5 -0.5],'g:','linewidth',2);
hold off;
ylim([-1,1]);
title('Control Input Signal');
legend('Simulated Control Input');
xlabel({'Time(sec)','(b)'});ylabel({'Offset-Free','Pump Voltage(V)'});
text(100,0.7,UMAX);text(100,-0.7,UMIN);
grid on

K_db=(acker(gg',hh',[0 0]))';
L_ndb3=acker(gg, cc,[0.9 0.9]);
sys_ndb3=ss(gg-cc*L_ndb3,cc,hh,dd,Ts);
gain_ndb3=dcgain(sys_ndb3);
y_pm=zeros(1,len);
y_ob=zeros(1,len);
x_pm=zeros(2,len);
x_ob=zeros(2,len);
u_ndb3=zeros(1,len);
x_pm(:,1)=x0;
for k=(1:len)
    u_ndb3(k)=Y_ref(k)/gain_ndb3-L_ndb3*x_ob(:,k);
    y_pm(k)=hh*x_pm(:,k);
    y_ob(k)=hh*x_ob(:,k);
    if k~=len
        x_pm(:,k+1)=gg*x_pm(:,k)+cc*u_ndb3(k);
        x_ob(:,k+1)=gg*x_ob(:,k)+cc*u_ndb3(k)+K_db*(y_pm(k)-y_ob(k));
    end
end
Error=x_pm-x_ob;
fig=figure;
set(fig,'position',[400 150 600 600]);
subplot(3,1,1);
stairs(T, y_pm,'r');
hold on;
stairs( T, Y_ref,'g','linewidth',2);
hold off;
set(gca,'YTick',[-1 -0.2 0 0.5 0.7]);
ylim([-1,1]);
grid on;
title({'Set-Point Control Results: Simulation','Output Signal'});
legend('Simulated Output','Reference Output');
xlabel({'Time(sec)','(a)'});ylabel({'Offset-Free','Water Level(V)'});
subplot(3,1,2);
stairs(T,u_ndb3,'b');
hold on;
plot ([0 525],[0.5 0.5],'g:','linewidth',2);
plot ([0 525],[-0.5 -0.5],'g:','linewidth',2);
hold off;
ylim([-1,1]);
title('Control Input Signal');
legend('Simulated Control Input');
xlabel({'Time(sec)','(b)'});ylabel({'Offset-Free','Pump Voltage(V)'});
text(100,0.7,UMIN);text(100,-0.7,UMAX);
grid on;
subplot(3,1,3);
stairs(T,Error(1,:),'r');
hold on;
stairs(T,Error(2,:));
hold off
axis([0 3 0 0.4]);grid on;
title('State estimation error');
legend('x_1(k)-x_1\^(k)','x_2(k)-x_2\^(k)');
xlabel({'Time(sec)','(c)'});ylabel('Estimation Error');


load('SFControlData_1.mat');
treal = SFLogData.time;
yref = SFLogData.signals(1).values(:,1);
yreal = SFLogData.signals(1).values(:,2);
ureal = SFLogData.signals(2).values;
figure;subplot(2,1,1);
stairs(T,y_pm+y_offset,'b');
hold on;
stairs(T,Y_ref+y_offset,'g');
stairs(treal,yreal,'r')
hold off;
legend('Simulated Output','Reference Output','Actual Output');
xlabel({'Time(sec)','(a)'});
ylabel('Water Level(V)');
subplot(2,1,2);
stairs(T,u_sp+u_offset,'b'); 
hold on;
stairs(treal,ureal,'r')
plot ([0 550],[2.5 2.55],'g:','linewidth',2);
plot ([0 550],[1.5 1.45],'g:','linewidth',2);
hold off; 
ylim([1,3]);
UMIN='V_m_i_n = 1.55';
UMAX='V_m_a_x = 2.45';
text(100,2.6,UMAX);
text(100,1.4,UMIN);
legend('Simulated Control Input','Actual Control Input');
xlabel({'Time(sec)','(b)'});
ylabel('Pump Voltage(V)');
