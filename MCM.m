clear 
clc
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MCM_initialise
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define simulation parameters

viz_iterations  = true;
save_iterations = true;

eco_pars.Nsupply    = 3;
eco_pars.kappa      = 0.01;

nyears = 10000;
years_per_iteration = 1;

eco_pars.seasonalcycle = false;

iseed = (find(eco_pars.unq_ESD<2)-1).*eco_pars.ntroph+1; % Strict Auto + ESD<2
% iseed = find(eco_pars.ESD<2); % All trophic strategies if ESD<2
% iseed = []; % Everything is Everywhere

Output_fdir = ['Run_N' num2str(eco_pars.Nsupply) '_Kappa_' num2str(eco_pars.kappa) '_' datestr(now,'yyyy_mm_dd')];
if exist(['Output/' Output_fdir])==0
    mkdir(['Output/' Output_fdir]);
    mkdir(['Output/' Output_fdir '/Figures']);
end
%% Initial conditions
N_0 = eco_pars.Nsupply;
P_0 = zeros(eco_pars.jpmax,1);

ndays_per_iteration=365*years_per_iteration; % needs to be an exact multiple of 365
n_ode_iterations=nyears./ndays_per_iteration*365;
tf=ndays_per_iteration;
eco_pars.t_res=7;
t0=0;
ndata_yr=numel([t0:eco_pars.t_res:tf]);
data=zeros(ndata_yr.*n_ode_iterations,eco_pars.jpmax+1);
data=sparse(data);


if isempty(iseed)
    N_0    = N_0./(eco_pars.jpmax); 
    P_0(:) = N_0; % seed with everything
    eco_pars.Pmut=speye(size(eco_pars.Pmut));
    minphy=eco_pars.extnct;
    eco_pars.minphy = minphy; % set minimum 'everything is everywhere' threshold
    eco_pars.extnct = eco_pars.extnct ./ 1000; % decrease 'extinction' threshold
else
    P_0(iseed) = eco_pars.seed_val;
    minphy=0;
    eco_pars.minphy = minphy;
end
v0=[N_0;P_0];

dead=[];

%% Initialise figure
if viz_iterations
    fh=figure(1);
    fh.Position=[103 1309 846 902];
    [ rgb, rgb3 ] = cmap_2d(eco_pars);
    set(gcf,'defaultAxesColorOrder',[0 0 0;rgb],...
            'Position',[50 231 1074 1114],...
            'Color','w')
    clf
end

%% Solve!
for k=1:n_ode_iterations
    tic
    % call ecosystem function
    [tout,yout] = ode45(@(t,y)  ecosystem(t,y,eco_pars,dead,minphy) ,[t0:eco_pars.t_res:tf],v0);
    
    if viz_iterations
        fh=figure(1);
        fh.Position=[103 1309 846 902];
        disp_indx=numel(tout);
        plot_output(tout,yout,eco_pars,disp_indx)
        if save_iterations
            fname=['Output/' Output_fdir '/Figures/FitLand_' num2str(k,'%05i') '.png'];
            saveas(gcf,fname);
        end
    end
    
    data((k-1).*ndata_yr+1:k.*ndata_yr,:)=yout;
    
    t0=tout(end);
    tf=t0+ndays_per_iteration;
    v0=yout(end,:);
    
    % Save Output
    save(['Output/' Output_fdir '/Workspace_dump.mat'],'data','eco_pars','nyears','years_per_iteration','iseed','dead','minphy')
    disp(['Year = ' num2str(k.*ndays_per_iteration./365) '; ' num2str(toc) ' seconds.'])
end





