close all
clear
clc

load 'SysIdenData_1.mat';

%pre.1
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
%G=tf([theta(3),theta(4)],[1,-theta(1),-theta(2)],Ts)

p=[1,-theta(1),-theta(2)];
eigen=roots(p)
error=0;
for i=1:length(eigen)
    if (eigen(i)>1)
        error=1;
        disp('This system is NOT stable')
    end
end
if (error==0)
    disp('This system is Stable')
end

z=zero(sys)
if (abs(z)<1)
    disp('This system is a Minimun Phase System')
end

%pre.2
rc=[h,g*h];
r_rc=rank(rc);
if (r_rc==2)
    disp('This system is Reachable and Controllable')
end

od=[c;c*g];
r_od=rank(od);
if (r_od==2)
    disp('This system is Obserbale and Detectable')
end

%pre.3
gg=[0,1;theta(2),theta(1)]';
hh=[0;1]';
cc=[theta(4),theta(3)]';
dd=0;
sys_o=ss(gg,cc,hh,dd,Ts)

%1.a
x0 = [0; 0.3];

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
figure(1);
subplot(211);
hold on;
axis([0 50 -1 1]);
stairs(T1,Y1,'r');
stairs(T2,Y2,'b');
title('Regularion Response by State Feedback y(k)');
legend ('Deadbated Response','Non-Deabeat Response y(k)');
grid on;
xlabel ({'Time(sec)','(a)'})
ylabel ('Offeset-Free Water Level(V)')

subplot(212);
hold on;
stairs(T1,U1,'r');
stairs(T1,U2,'b');
u_min = ones(1, length(U1))*0.5;
u_max = ones(1, length(U1))*-0.5;
plot(T1,u_min,'g:','linewidth',2);
plot(T1,u_max,'g:','linewidth',2);
axis([0 50 -1 1])
xlabel ({'Time(sec)','(b)'});
ylabel ('Offeset-Free Pump Voltage(V)');
grid on;
UMIN=['u_m_i_n = ' num2str(u_min(1))];
UMAX=['u_m_a_x = ' num2str(u_max(1))];
text(5,0.65,UMAX);
text(5,-0.65,UMIN);

%2.a
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
%2.b
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
% 
%3.a
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

% %post1
% dist=[zeros(1,Period) ones(1,len-Period)];
% W=[0.5;0];
% y_sp2=zeros(1,len);
% x_sp2=zeros(2,len);
% u_sp2=zeros(1,len);
% for k=(1:len)
%     u_sp2(k)=Y_ref(k)/gain_ndb1-L_ndb1*x_sp2(:,k);
%     y_sp2(k)=hh*x_sp2(:,k);
%     if k ~=len
%         x_sp2(:,k+1)=gg*x_sp2(:,k)+cc*u_sp2(k)+W*dist(:,k);
%     end
% end
% figure;
% subplot(2,1,1);
% stairs(T, y_sp2,'r');
% hold on;
% stairs(T, Y_ref,'g');
% hold off;
% grid on;
% title({'Set-Point Control Results: Simulation','Output Signal'});
% legend('Simulated Output','Reference Output');
% xlabel({'Time(sec)','(a)'});ylabel({'Offset-Free','Water Level(V)'});
% subplot(2,1,2);
% stairs(T, u_sp2,'b');
% hold on;
% plot ([0 550],[0.5 0.5],'g:','linewidth',2);
% plot ([0 550],[-0.5 -0.5],'g:','linewidth',2);
% hold off;
% legend('Simulated Control Input');
% xlabel({'Time(sec)','(b)'});ylabel({'Offset-Free','Pump Voltage(V)'});
% %text(100,0.6,UMIN);text(100,-0.6,UMAX);
% grid on;
% aa=[0.001*(rand-0.5) 0.01*(rand-0.5) 0.1*(rand-0.5)];
% for i=1:3
%     %a=0.01*(rand-0.5);
%     a=aa(i);
%     ggg=[0,1;theta(2)-a,theta(1)-a]';
%     ccc=[theta(4)+a,theta(3)+a]';
%     sys_rob=ss(ggg-ccc*L_ndb1,ccc,hh,dd,Ts);
%     gain_rob=dcgain(sys_rob);
%     [y_rob,T,x_rob]=lsim(sys_rob,Y_ref*1/gain_rob);
%     u_rob=-L_ndb1*x_rob'+1/gain_rob*Y_ref;
%     uncertainly=abs(a)*100;
%     str=['uncertainly=',num2str(uncertainly),'%'];
%     figure;
%     subplot(2,1,1);
%     plot(T, y_rob,'r', T, y_sp,'b',T, Y_ref, 'g');
%     set(gca,'YTick',[-1 -0.5 0 0.2 0.7]);
%     grid on;
%     text(15,0.85,str);
%     title({'Set-Point Control Results: Simulation','Output Signal'});
%     legend('Robust Behaviour Output','Set-point Output','Reference Output');
%     xlabel({'Time(sec)','(a)'});ylabel({'Offset-Free','Water Level(V)'});
%     subplot(2,1,2);plot(T, u_rob,'r',T, u_sp,'b');
%     hold on;
%     plot ([0 550],[0.5 0.5],'g:','linewidth',2);
%     plot ([0 550],[-0.5 -0.5],'g:','linewidth',2);
%     hold off;
%     title('Control Input Signal');
%     legend('Robust Behaviour Input','Simulated Control Input');
%     xlabel({'Time(sec)','(b)'});ylabel({'Offset-Free','Pump Voltage(V)'});
%     text(100,0.6,UMIN);text(100,-0.6,UMAX);
%     grid on;
% end