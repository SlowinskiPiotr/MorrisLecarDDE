close all
clear

% set structure with rhs and derivatives functions for continuation
% functions generated using maple script: ml_dde_gen_my_sys.mw
addpath('./System_Files/'); % add path to files with the system
ml_dde_funcs=set_funcs('sys_rhs', @sys_rhs,...
                       'sys_tau', @()17,...
                       'sys_deri', @sys_deri,...
                       'sys_dtau', @sys_dtau,...
                       'sys_mfderi', @sys_mfderi,...
                       'sys_ntau', @()1);
                   
%% change default method values
method=df_mthod(ml_dde_funcs,'stst',1);
method.continuation.plot=0;
method.stability.minimal_real_part=-2;
method.stability.max_number_of_eigenvalues=10000;
method.stability.root_accuracy=1e-8;
method.stability.interpolation_order=6;
method.stability.minimal_time_step=0.001;
method.stability.minimal_time_step=0.01;

my_methods=method;

%% setting parameters
phi=0.125;  %1
gca=4.0;    %2
gl=2.0;     %3 
gk=8.0;     %4 
v1=-2.8;    %5 
v2=26.0;    %6 
v3=12.0;    %7
v4=17.4;    %8
vl=-60.0;   %9
vca=120.0;  %10
vk=-91.893652410657; %11
C=20.0;     %12
I=-60.0;    %13
Kappa=1.8;  %14 %arbitrary choice so that the Hopf point is near I=15
v_s=0.0;    %15
Kappa_s=5.0;%16
Tau=10.0;   %17
parameters=[phi, gca, gl, gk, v1, v2, v3, v4, vl, vca, vk, C, I, Kappa, v_s, Kappa_s, Tau];

ind_C=12;
ind_I=13;
ind_Kappa=14;
ind_Tau=17;

% 1st 1D continuation in I
x0=[-100;0]; %initial guess of an equilibrium
stst_I=SetupStst(ml_dde_funcs,'parameter',parameters,'x',x0,'contpar',ind_I,'step',0.01);

stst_I.parameter.max_step=[13,  0.1];
stst_I.parameter.min_bound=[13,  -60];
stst_I.parameter.max_bound=[13,   60];

stst_I.method=my_methods;

tic,
stst_I=br_contn(ml_dde_funcs,stst_I,10000);
stst_I=br_rvers(stst_I);
stst_I=br_contn(ml_dde_funcs,stst_I,10000);
toc,
tic,
stst_I=br_stabl(ml_dde_funcs,stst_I,0,0);
toc,

%% plot the computed branch
figure(1)
[x_I,x_var]=df_measr(0,stst_I);
br_splot(stst_I,x_I,x_var)

%% contiunation of the Hopf branch
nunst=GetStability(stst_I,'funcs',ml_dde_funcs);
indhopf=find(abs(diff(nunst))==2);

hopf_tk=SetupHopf(ml_dde_funcs,stst_I,indhopf,'contpar',[ind_Kappa,ind_Tau],'dir',14,'step',0.01,'plot',0);

hopf_tk.parameter.free=[ind_Kappa ind_Tau];
hopf_tk.parameter.max_step=[ind_Kappa 0.1;ind_Tau 0.1];    
hopf_tk.parameter.min_bound=[ind_Kappa 0; ind_Tau -20];
hopf_tk.parameter.max_bound=[ind_Kappa 300; ind_Tau 300];

tic,
hopf_tk=br_contn(ml_dde_funcs,hopf_tk,35000);
hopf_tk=br_rvers(hopf_tk);
hopf_tk=br_contn(ml_dde_funcs,hopf_tk,35000);
toc,

%% plot the computed branch
[x_kappa,x_tau]=df_measr(0,hopf_tk);
br_plot(hopf_tk,x_tau,x_kappa)

%% initial PO from hopf for kappa=60
br_PO_k60_t=[];

kappa_val=60;

ind=find(diff(sign(arrayfun(@(x) x.parameter(14), hopf_tk.point)-kappa_val))~=0);
points=arrayfun(@(x) x, hopf_tk.point(ind));
intervals=64;
degree=3;
k=1;

for i=length(ind):-1:1
    if points(i).parameter(17)>0 && points(i).parameter(17)<140
        try  % just in case any the branches couldn's be started  (e.g. stariting with tau)
            points(i).parameter(14)=kappa_val; %we fix kappas
            
            method=df_mthod(ml_dde_funcs,'hopf',1); % get hopf calculation method parameters:
            [points(i),~]=p_correc(ml_dde_funcs,points(i),ind_Tau,[],method.point); % correct hopf for kappa=60
            [psol,stepcond]=p_topsol(ml_dde_funcs,points(i),1e-2,degree,intervals); % hopf to psol
            
            %correct periodic solution guess:
            method=df_mthod(ml_dde_funcs,'psol');
            psol=p_correc(ml_dde_funcs,psol,17,stepcond,method.point);
            
            br_PO_k60_t(k).branch=df_brnch(ml_dde_funcs,17,'psol'); % empty branch:
            % make degenerate periodic solution with amplitude zero at hopf point:
            deg_psol=p_topsol(ml_dde_funcs,points(i),0,degree,intervals);
            
            % use deg_psol and psol as first two points on branch:
            br_PO_k60_t(k).branch.point=deg_psol;
            br_PO_k60_t(k).branch.point(2)=psol;
            br_PO_k60_t(k).branch.method.continuation.plot=0;
            
            br_PO_k60_t(k).branch.parameter.min_bound=[ind_Tau  0];
            br_PO_k60_t(k).branch.parameter.max_bound=[ind_Tau  200];
            br_PO_k60_t(k).branch.parameter.max_step =[ind_Tau  0.2];
            
            tic,
            br_PO_k60_t(k).branch=br_contn(ml_dde_funcs,br_PO_k60_t(k).branch,550);
            toc,
            
            [~,p_ampl]=df_measr(0,br_PO_k60_t(k).branch);
            
            br_plot3(br_PO_k60_t(k).branch,x_tau,x_kappa,p_ampl)
            drawnow,
            pause(1),
            disp(k);
            k=k+1;
        catch
        end
    end
end
%% visual inspection to choose branches and tau values
for i=1:15
    br_plot(br_PO_k60_t(i).branch,x_tau,p_ampl,'k')
end
% selected in colour (we choose the 1st branches and and a point with a
% quite low amplitude)
for i=[2 4 6 8 10 12 13 14 15]
    br_plot(br_PO_k60_t(i).branch,x_tau,p_ampl) 
end

%% continue the solution in C and I to get to the values we want to analyse
% solution and points have to be seelected by hand (not all the points can
% be continued all the way to C=0.5
x_C=x_kappa;
x_C.col=12;
p_period=p_ampl;
p_period.field='period';
p_period.col=1;

t=  [11.5 26 42 57 72 86 102 116 132];
brn=[2    4  6  8  10 12 13  14  15];

move_PO_C=[];
move_PO_I=[];

for k=numel(t):-1:1
    i_brn=brn(k);
    ind=find(diff(sign(arrayfun(@(x) x.parameter(17), br_PO_k60_t(i_brn).branch.point)-t(k)))~=0);
    if ~isempty(ind)
        ind=ind(1);
        try
            disp(k),
            disp(ind),
            disp(i_brn),
            disp(t(k)),
            move_PO_C(k).branch=br_PO_k60_t(i_brn).branch;
            move_PO_C(k).branch.method.point.adapt_mesh_before_correct=1;
            move_PO_C(k).branch.method.point.adapt_mesh_after_correct=1;
            
            move_PO_C(k).branch.point=[];
            move_PO_C(k).branch.parameter.free=ind_C;
            move_PO_C(k).branch.parameter.max_step=[ind_C 0.1];
            move_PO_C(k).branch.parameter.min_bound=[ind_C 0.5]; % this is the value of C that we want
            move_PO_C(k).branch.parameter.max_bound=[ind_C 80];
            
            p_ini=p_remesh(br_PO_k60_t(i_brn).branch.point(ind),4,64);
            
            [~,p_ampl]=df_measr(0,br_PO_k60_t(i_brn).branch);
            
            tic,
            move_PO_C(k).branch=correct_ini(ml_dde_funcs,move_PO_C(k).branch,p_ini,12,-0.1,1);
            move_PO_C(k).branch=br_contn(ml_dde_funcs,move_PO_C(k).branch,5000);
            toc,
            
            figure(1) %to see the branches that are being computed
            subplot(2,2,1)
            br_plot(br_PO_k60_t(i_brn).branch,x_tau,p_ampl)
            subplot(2,2,2)
            br_plot(move_PO_C(k).branch,x_C,p_ampl)
            subplot(2,2,3)
            br_plot(br_PO_k60_t(i_brn).branch,x_tau,p_period)
            subplot(2,2,4)
            br_plot(move_PO_C(k).branch,x_C,p_period)
            drawnow,
            pause(1)
        catch
        end
    end
    try
        move_PO_I(k).branch=move_PO_C(k).branch;
        move_PO_I(k).branch.method.point.adapt_mesh_before_correct=1;
        move_PO_I(k).branch.method.point.adapt_mesh_after_correct=1;
        
        move_PO_I(k).branch.point=[];
        move_PO_I(k).branch.parameter.free=ind_I;
        move_PO_I(k).branch.parameter.max_step=[ind_I 0.2];
        move_PO_I(k).branch.parameter.min_bound=[ind_I -37]; % this is the final value of I that we want
        move_PO_I(k).branch.parameter.max_bound=[ind_I 140];
        
        p_ini=move_PO_C(k).branch.point(end);

        tic,
        move_PO_I(k).branch=correct_ini(ml_dde_funcs,move_PO_I(k).branch,p_ini,13,-0.01,1);
        move_PO_I(k).branch=br_contn(ml_dde_funcs,move_PO_I(k).branch,2500);
        toc,
        figure(2)
        br_plot3(move_PO_I(k).branch,x_C,x_I,p_ampl)
        drawnow
        pause(1),
    catch
    end
end
%% now we continue the branches with I=-37 and C=0.5  in tau
PO_t=[];

for k=1:9
    PO_t(k).branch=move_PO_I(k).branch;
    % we increase accuracy of the computations
    
    PO_t(k).branch.method.point.adapt_mesh_before_correct=1;
    PO_t(k).branch.method.point.adapt_mesh_after_correct=1;
    PO_t(k).branch.method.point.minimal_accuracy=1e-9;
    PO_t(k).branch.method.point.halting_accuracy=1e-11;

    PO_t(k).branch.point=[];
    PO_t(k).branch.parameter.free=ind_Tau;
    PO_t(k).branch.parameter.max_step=[ind_Tau 0.01];
    
    PO_t(k).branch.parameter.min_bound=[ind_Tau 0];
    PO_t(k).branch.parameter.max_bound=[ind_Tau max([55 t(k)+1])];
    
    p_ini=move_PO_I(k).branch.point(end);
    
    tic,
    PO_t(k).branch=correct_ini(ml_dde_funcs,PO_t(k).branch,p_ini,17,0.01,1);
    PO_t(k).branch=br_contn(ml_dde_funcs,PO_t(k).branch,2500);
    toc,
    tic,
    PO_t(k).branch=br_rvers(PO_t(k).branch);
    PO_t(k).branch=br_contn(ml_dde_funcs,PO_t(k).branch,2500);
    toc,
end

%% we compute stability and trim to tau=57 
% (to save time, originally we needed higher values of tau to start the continution)

for k=1:9 
    idx1=find(arrayfun(@(x) x.parameter(17), PO_t(k).branch.point)<57,1,'first');
    idx2=find(arrayfun(@(x) x.parameter(17), PO_t(k).branch.point)<57,1,'last');
    PO_t(k).branch.point=PO_t(k).branch.point(idx1:idx2);
    
    tic
    % we increase accuracy of the stability of computations
    PO_t(k).branch.method.stability.max_number_of_eigenvalues=10000;
    PO_t(k).branch.method.stability.minimal_modulus=1e-8;
    PO_t(k).branch=br_stabl(ml_dde_funcs,PO_t(k).branch,0,1);
    toc
end
%% plot the branches with stability information
clf
for k=1:9
    figure(3)
    subplot(1,2,1)
    br_splot(PO_t(k).branch,x_tau,p_ampl),
    axis square
    xlim([0 55])
    subplot(1,2,2)
    br_splot(PO_t(k).branch,x_tau,p_period),
    axis square
    axis([0 55 0 55])
end

%% the 1st branch of the folds is tricky to compute we do it later
for k=2:9
        [tau_v,fold_PO_idx]=min(arrayfun(@(x) x.parameter(17), PO_t(k).branch.point));
    
    tic,
    
    [ml_dde_foldfuncs,fold_PO_tk(k).branch]=SetupPOfold(ml_dde_funcs,PO_t(k).branch,fold_PO_idx,...
        'contpar',[ind_Kappa, ind_Tau],...
        'dir',ind_Kappa,'step',-0.1,...
        'min_bound',[ind_Kappa,  0;   ind_Tau,0],...
        'max_bound',[ind_Kappa, 70;  ind_Tau,55],...
        'max_step', [ind_Kappa,0.5; ind_Tau,0.5]);
    
    fold_PO_tk(k).branch.method.continuation.plot=0;
    fold_PO_tk(k).branch.method.continuation.plot_progress=0;
    fold_PO_tk(k).branch.method.point.print_residual_info=0;
    fold_PO_tk(k).branch.method.point.adapt_mesh_before_correct=1;
    fold_PO_tk(k).branch.method.point.adapt_mesh_after_correct=1;
    
    fold_PO_tk(k).branch=br_contn(ml_dde_foldfuncs,fold_PO_tk(k).branch,2500);
    fold_PO_tk(k).branch=br_rvers(fold_PO_tk(k).branch);
    fold_PO_tk(k).branch=br_contn(ml_dde_foldfuncs,fold_PO_tk(k).branch,2500);
    disp(k)    
    toc,
    
    figure(1)
    br_plot(fold_PO_tk(k).branch,x_tau, x_kappa)
    drawnow
    pause(1)
end
%% to check that the 1st fold branch is indeed a border of a region of stable solutions we use it initiate PO 
% and to find a better point to start the computations 

taus=5:5:35;

br_plot3(PO_t(1).branch,x_tau,x_kappa,p_period) %a modified plotting function

for kk=7:-1:1
    ind=find(diff(sign(arrayfun(@(x) x.parameter(17), PO_t(1).branch.point)-taus(kk)))~=0); 
    % find points that are very close to the taus values
    
    % define new branch
    PO_k(kk).branch=PO_t(1).branch;
    PO_k(kk).branch.point=[];
    PO_k(kk).branch.parameter.free=ind_Kappa;

    p_ini=PO_t(k).branch.point(ind);
    
    tic,
    PO_k(kk).branch=correct_ini(ml_dde_funcs,PO_t(1).branch,p_ini,ind_Kappa,0.01,1);
    
    % we increase accuracy of the computations
    PO_k(kk).branch.point(1)=p_remesh(PO_k(kk).branch.point(1),4,64);
    PO_k(kk).branch.point(2)=p_remesh(PO_k(kk).branch.point(2),4,64);
    
    PO_k(kk).branch.method.point.adapt_mesh_before_correct=1;
    PO_k(kk).branch.method.point.adapt_mesh_after_correct=1;
    PO_k(kk).branch.method.point.minimal_accuracy=1e-10;
    PO_k(kk).branch.method.point.halting_accuracy=1e-12;
    PO_k(kk).branch.method.stability.max_number_of_eigenvalues=10000;
    PO_k(kk).branch.method.stability.minimal_modulus=1e-8;
    
    PO_k(kk).branch.parameter.max_step=[ind_Kappa 0.01];
    PO_k(kk).branch.parameter.min_bound=[ind_Kappa 0];
    PO_k(kk).branch.parameter.max_bound=[ind_Kappa 70];
    
    tic,
    PO_k(kk).branch=br_contn(ml_dde_funcs,PO_k(kk).branch,2500);
    toc,
    tic,
    PO_k(kk).branch=br_rvers(PO_k(kk).branch);
    PO_k(kk).branch=br_contn(ml_dde_funcs,PO_k(kk).branch,1500);
    toc,
    tic,
    PO_k(kk).branch=br_stabl(ml_dde_funcs,PO_k(kk).branch,0,0);
    toc,
    br_splot3(PO_k(kk).branch,x_tau,x_kappa,p_period) %a modified plotting function
    drawnow,
    pause(0.5)
end
%% the same but now we continue in kappa
kappas=50:5:70;

br_plot3(PO_t(1).branch,x_tau,x_kappa,p_period)
br_splot3(PO_k(1).branch,x_tau,x_kappa,p_period)

for kk=5:-1:1
    ind=find(diff(sign(arrayfun(@(x) x.parameter(14), PO_k(2).branch.point)-kappas(kk)))~=0);
    PO_t_extra(kk).branch=PO_k(2).branch;
    PO_t_extra(kk).branch.point=[];
    
    PO_t_extra(kk).branch.parameter.free=ind_Tau;
    PO_t_extra(kk).branch.parameter.max_step=[ind_Tau 0.1];
    PO_t_extra(kk).branch.parameter.min_bound=[ind_Tau 0];
    PO_t_extra(kk).branch.parameter.max_bound=[ind_Tau 55];
    
    % we increase accuracy of the computations
    PO_t_extra(kk).branch.method.point.adapt_mesh_before_correct=1;
    PO_t_extra(kk).branch.method.point.adapt_mesh_after_correct=1;
    PO_t_extra(kk).branch.method.point.minimal_accuracy=1e-9;
    PO_t_extra(kk).branch.method.point.halting_accuracy=1e-11;
    PO_t_extra(kk).branch.method.stability.max_number_of_eigenvalues=10000;
    PO_t_extra(kk).branch.method.stability.minimal_modulus=1e-8;
    
    p_ini=PO_k(2).branch.point(ind(1));
    
    PO_t_extra(kk).branch=correct_ini(ml_dde_funcs,PO_t_extra(kk).branch,p_ini,ind_Tau,0.0001,1);
    PO_t_extra(kk).branch.point(1)=p_remesh(PO_t_extra(kk).branch.point(1),4,64);
    PO_t_extra(kk).branch.point(2)=p_remesh(PO_t_extra(kk).branch.point(2),4,64);
    
    tic,
    PO_t_extra(kk).branch=br_contn(ml_dde_funcs,PO_t_extra(kk).branch,2500);
    toc,
    tic,
    PO_t_extra(kk).branch=br_rvers(PO_t_extra(kk).branch);
    PO_t_extra(kk).branch=br_contn(ml_dde_funcs,PO_t_extra(kk).branch,500);
    toc,
    tic,
    PO_t_extra(kk).branch=br_stabl(ml_dde_funcs,PO_t_extra(kk).branch,0,0);
    toc,
    
    br_splot3(PO_t_extra(kk).branch,x_tau,x_kappa,p_period)
    zlim([0 55])
    drawnow,
    pause(0.5)
end
%% here we actually compute the 1st fold

k=1
[~,fold_PO_idx]=min(arrayfun(@(x) x.parameter(17), PO_t_extra(1).branch.point));

[ml_dde_foldfuncs,fold_PO_tk(k).branch]=SetupPOfold(ml_dde_funcs,PO_t_extra(1).branch,fold_PO_idx,...
    'contpar',[ind_Kappa, ind_Tau],...
    'dir',ind_Tau,'step',-0.0001,...
    'min_bound',[ind_Kappa,  0;   ind_Tau,0],...
    'max_bound',[ind_Kappa, 71;  ind_Tau,55],...
    'max_step', [ind_Kappa,0.5; ind_Tau,0.5]);

% we increase accuracy of the computations
fold_PO_tk(k).branch.point(1)=p_remesh(fold_PO_tk(k).branch.point(1),4,64);
fold_PO_tk(k).branch.point(2)=p_remesh(fold_PO_tk(k).branch.point(2),4,64);

fold_PO_tk(k).branch.method.continuation.plot=0;
fold_PO_tk(k).branch.method.continuation.plot_progress=0;
fold_PO_tk(k).branch.method.point.print_residual_info=0;
fold_PO_tk(k).branch.method.point.adapt_mesh_before_correct=1;
fold_PO_tk(k).branch.method.point.adapt_mesh_after_correct=1;
fold_PO_tk(kk).branch.method.point.minimal_accuracy=1e-10;
fold_PO_tk(kk).branch.method.point.halting_accuracy=1e-12;

tic,
fold_PO_tk(k).branch=br_contn(ml_dde_foldfuncs,fold_PO_tk(k).branch,2500);
fold_PO_tk(k).branch=br_rvers(fold_PO_tk(k).branch);
toc,

tic,
fold_PO_tk(k).branch=br_contn(ml_dde_foldfuncs,fold_PO_tk(k).branch,2500);
toc,


%% here we compute the period-doubling bifurcation
tic
nunst=GetStability(PO_t_extra(3).branch,'funcs',ml_dde_funcs,'exclude_trivial',1);
PD_PO_idx=find(nunst>0,1,'first'),

% the continuation is tricky to initiate so we give it a bit more time to
% converge

PO_t_extra(3).branch.method.point.newton_max_iterations=20;
PO_t_extra(3).branch.method.point.newton_nmon_iterations=10;

[PDfuncs,branchPD]=SetupPeriodDoubling(ml_dde_funcs,PO_t_extra(3).branch,PD_PO_idx,'contpar',[14,17],...
    'dir',14,'step',1e-4,'min_bound',[14,0; 17,0],'max_bound',[14, 71; 17,55],...
    'max_step',[14,0.1; 17,0.1]);

branchPD.method.point.newton_max_iterations=20;
branchPD.method.point.newton_nmon_iterations=10;
branchPD.branch.method.point.minimal_accuracy=1e-10;
branchPD.branch.method.point.halting_accuracy=1e-12;

branchPD.method.point.print_residual_info=0;

tic,
branchPD=br_contn(PDfuncs,branchPD,100);
branchPD=br_rvers(branchPD);
branchPD=br_contn(PDfuncs,branchPD,100);
toc,

%% plot the first fold branches, the period doubling branch and the solution branches
for k=1:5 
    br_splot3(PO_t_extra(k).branch,x_tau,x_kappa,p_period)
end
for k=1:7 
    br_splot3(PO_k(k).branch,x_tau,x_kappa,p_period)
end
br_plot3(branchPD,x_tau,x_kappa,p_period,'k')
br_plot3(fold_PO_tk(1).branch,x_tau,x_kappa,p_period)

%% plot all the fold branches and the period doubling branch
for k=1:9
    br_plot(fold_PO_tk(k).branch,x_tau, x_kappa)
end
br_plot(branchPD,x_kappa,p_period,'k')