%close all
clear

addpath('./System_Files/'); % add path to files with the system

phi=0.0125;
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
C=5;
I=-37;
Kappa=60;
v_s=0.0;
Kappa_s=5.0;
par=[phi, gca, gl, gk, v1, v2, v3, v4, vl, vca, vk, C, I, Kappa, v_s, Kappa_s];

clc

tau=400;

period=tau;
tspan=linspace(0,tau,101);
t_start=0;
t_end=t_start+5000;

% we set the stimualtion as a history
history(1).h=zeros(2,101)-80;
history(1).h(2,:)=zeros(1,101)+0.01;
history(1).h(1,(1:4)+1)=50;
history(1).h(1,(1:4)+15)=50;

history_fun = @(t) interp1(tspan'-period, history(1).h', -mod(-t,period))';  %defining the required inputs
tic,
sol(1).s = dde23(@(t,xx,yy)ml_dde_rhs(t,xx,yy,par), tau, history_fun, [t_start, t_end]);
toc,

h1=subplot(2,1,1);
t=sol(1).s.x(:);
x=sol(1).s.y(1,:);

plot((tspan-400),history(1).h(1,:),'color','r')
hold on
plot(t,x,'color',lines(1))
hold off

xlim([-405,2000])
set(h1,'Xtick',-400:400:2000,'XtickLabel',-400:400:2000)
xlabel('Time (ms)');
ylabel('V');

history(2).h=zeros(2,101)-80;
history(2).h(2,:)=zeros(1,101)+0.01;
history(2).h(1,(1:4)+1)=50;
history(2).h(1,(1:4)+31)=50;

history_fun = @(t) interp1(tspan'-period, history(2).h', -mod(-t,period))';  %defining the required inputs
tic,
sol(2).s = dde23(@(t,xx,yy)ml_dde_rhs(t,xx,yy,par), tau, history_fun, [t_start, t_end]);
toc,

h2=subplot(2,1,2);
t=sol(2).s.x(:);
x=sol(2).s.y(1,:);

plot((tspan-400),history(2).h(1,:),'color','r')
hold on
plot(t,x,'color',lines(1))
hold off

xlim([-405,2000])
set(h2,'Xtick',-400:400:2000,'XtickLabel',-400:400:2000)
xlabel('Time (ms)');
ylabel('V');

[~,idx]=findpeaks(sol(1).s.y(1,:),'MinPeakHeight',10);
ddd=diff(sol(1).s.x(idx));
ddd(end-10:end)
[~,idx]=findpeaks(sol(2).s.y(1,:),'MinPeakHeight',10);
ddd=diff(sol(2).s.x(idx));
ddd(end-10:end)
%% long time simulation by concatenating short ones

for i_sol=1:2
    ddd=[];
    sols_1k(1,i_sol).sols=sol(i_sol).s;
    sol_h=sols_1k(1,i_sol).sols;
    
    [~,idx]=findpeaks(sol_h.y(1,:),'MinPeakHeight',10);
    i_last_peak=idx(end-1);
    t_last_peak=sol_h.x(i_last_peak-100);
    
    history_fun = @(t) interp1(sol_h.x(i_last_peak-600:i_last_peak-100)'-t_last_peak, sol_h.y(:,i_last_peak-600:i_last_peak-100)', -mod(-t,480),'linear');  %defining the required inputs
    
    for k=1:4000
        k,
        
        tic,
        sol_end = dde23(@(t,xx,yy)ml_dde_rhs(t,xx,yy,par), tau, history_fun, [0, 10000]);
        toc
        
        sols_1k(k+1,i_sol).sols=sol_end;
        
%                 if k>1
%                     plot(sol_h.x(i_last_peak-600:i_last_peak-100),sol_h.y(1,i_last_peak-600:i_last_peak-100))
%                     hold on
%                     plot(sol_end.x(1:500)+t_last_peak,sol_end.y(1,1:500))
%                     hold off,
%                     drawnow,
%                     pause(1)
%                 end
        
        [~,idx]=findpeaks(sol_end.y(1,:),'MinPeakHeight',10);
        i_last_peak=idx(end-1);
        t_last_peak=sol_end.x(i_last_peak-100);
        
        sol_h=sol_end;
        history_fun = @(t) interp1(sol_h.x(i_last_peak-600:i_last_peak-100)'-t_last_peak, sol_h.y(:,i_last_peak-600:i_last_peak-100)', -mod(-t,480),'linear');  %defining the required inputs
        
        
        [~,idx]=findpeaks(sol_end.y(1,:),'MinPeakHeight',10);
        ddd=[ddd diff(sol_end.x(idx))];
        ddd(end-10:end)
    end
end
%%