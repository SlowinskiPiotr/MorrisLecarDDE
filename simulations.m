close all
clear

addpath('./System_Files/'); % add path to files with the system

phi=0.125;
gca=4.0;
gl=2.0;
gk=8.0;
v1=-2.8;
v2=26.0;
v3=12.0;
v4=17.4;
vl=-60.0;
vca=120.0;
vk=-91.893652410657;
C=0.5;
I=-37;
Kappa=60;
v_s=0.0;
Kappa_s=5.0;
par=[phi, gca, gl, gk, v1, v2, v3, v4, vl, vca, vk, C, I, Kappa, v_s, Kappa_s];

clc

tau=40;

period=40;
tspan=linspace(0,40,101);
t_start=0;
t_end=t_start+5000;

% we set the stimualtion as a history
history(1).h=zeros(2,101)-80; 
history(1).h(2,:)=zeros(1,101)+0.01;
history(1).h(1,(1:4)+5)=50;
history(1).h(1,(1:4)+19)=50;

history_fun = @(t) interp1(tspan'-period, history(1).h', -mod(-t,period))';  %defining the required inputs
tic,
sol(1).s = dde23(@(t,xx,yy)ml_dde_rhs(t,xx,yy,par), tau, history_fun, [t_start, t_end]);
toc,

subplot(1,2,1)
t=sol(1).s.x(:);
x=sol(1).s.y(1,:);

plot(t,x,'color',lines(1))
hold off

xlim([0,500])
xlabel('Time (ms)');
ylabel('V');

history(2).h=zeros(2,101)-80;
history(2).h(2,:)=zeros(1,101)+0.01;
history(2).h(1,(1:4)+5)=50;
history(2).h(1,(1:4)+20)=50;


history_fun = @(t) interp1(tspan'-period, history(2).h', -mod(-t,period))';  %defining the required inputs
tic,
sol(2).s = dde23(@(t,xx,yy)ml_dde_rhs(t,xx,yy,par), tau, history_fun, [t_start, t_end]);
toc,

subplot(1,2,2)
t=sol(2).s.x(:);
x=sol(2).s.y(1,:);

plot(t,x,'color',lines(1))
hold off

xlim([0,500])
xlabel('Time (ms)');
ylabel('V');

[~,idx]=findpeaks(sol(1).s.y(1,:),'MinPeakHeight',10);
ddd=diff(sol(1).s.x(idx));      
ddd(end-10:end)
[~,idx]=findpeaks(sol(2).s.y(1,:),'MinPeakHeight',10);
ddd=diff(sol(2).s.x(idx));
ddd(end-10:end)