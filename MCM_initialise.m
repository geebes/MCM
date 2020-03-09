%%  
addpath Diagnostics

%%  Environmental parameters

eco_pars.Nsupply    = 5;
eco_pars.kappa      = 0.01;

% growth limitation term
eco_pars.gamma = 0.1;

eco_pars.gamma_cycle = (cos(linspace(0,2.*pi,365+1))+1)./2;
eco_pars.gamma_cycle = 0.1+eco_pars.gamma_cycle/10;

%% ECOSYSTEM & EVOLUTION PARAMETERS

% Trait dimensions
eco_pars.nsize      = 51;                                                  % total number of plankton size classes
eco_pars.ntroph     = 51;                                                  % total number of plankton trophic classes
eco_pars.tau        = 1.0;                                                 % trophic trade-off parameter

% quota
eco_pars.muinf_a    = +4.70;   % /d                                        % Ward et al. AmNat 2017
eco_pars.muinf_b    = -0.26;                                               % Ward et al. AmNat 2017

eco_pars.Vmax_a     = +0.024;  % molN/cell/d                               % Ward et al. AmNat 2017
eco_pars.Vmax_b     = +1.10;         

eco_pars.Qmin_a     = +0.032;  % molN/cell                                 % Ward et al. AmNat 2017
eco_pars.Qmin_b     = 0.76;                 

eco_pars.alpha_a    = +794;    % m3/mmolN/d
eco_pars.alpha_b    = -0.63;   % 

eco_pars.linearmort = +0.025;  % /d                                        % Maranon et al EcoLett 2013

% grazing
eco_pars.g_a        = +6.48; % m3 / mmol N / d                             % coefficient of attack rate
eco_pars.g_b        = -0.16;                                               % exponent of attack rate
eco_pars.refuge     = -100;                                                % grazing refuge parameter (negative)
                                                                           
eco_pars.pp_opt     = +1000;                                               % optimum pred:prey ratio
eco_pars.pp_sig     = +2;                                                  % width grazing kernel
eco_pars.lambda     = +0.7;                                                % grazing efficiency

% initialisation
eco_pars.PHY_init   = 0.0;                                                 % initial biomass *everywhere*
eco_pars.seed_loc   = [];                                                  % spatial indices for seeding
eco_pars.seed_PHY   = [];                                                  % plankton indices for seeding
eco_pars.seed_val   = 1e-12;                                               % biomass for seeding

% evolution (trophic strategy, size)
eco_pars.mutrat     = [1e-15 1e-15];                                       % mutation rate relative to size of trait space (will be normalised)
eco_pars.extnct     = 1e-15;                                               % extinction threshold

% remineralisation
eco_pars.Dremin = 0.01;

% plankton advection scaling
eco_pars.trscale=1;
eco_pars.ngenes=50;



eco_pars.jpmax=eco_pars.nsize.*eco_pars.ntroph;


%% 


% Trait dimensions
% size
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Need to make sure size classes 'snap' to optimal pp interval
n_OoM=min(eco_pars.nsize-1,4); % Minimum number of orders of magnitude to span in length dimension
if eco_pars.nsize==1
    minsize=6; % only size class
elseif eco_pars.nsize<=5
    minsize=0.6; % size of smallest class
else
    minsize=0.06; % size of smallest class
    n_OoM=5;
end
l_ratio=log10(nthroot(eco_pars.pp_opt,3)); % find log10 of optimal length ratio (as cube root of volume ratio)
proceed=1; % n size classes in optimal pp interval
d=1;
while proceed % If smaller spacing can still achieve required length span...
    f=l_ratio/d; % calculate log spacing
    if (l_ratio/(d+1)) >= (n_OoM/(eco_pars.nsize-1))
        d=d+1; % add one more size class in optimal pp interval
    else
        proceed=0;
    end % stop iteration if spacing too close
end
ESD = minsize.*10.^([0:(eco_pars.nsize-1)].*f)'; % calculate ESDs of all plankton classes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     ESD  = logspace(log10(0.6),log10(6000),eco_pars.nsize)'; % equivalent spherical diameter
ESR  = 0.5.*ESD;     % equivalent spherical radius
eco_pars.V  = 4/3*pi.*ESR.^3; % cell volume
eco_pars.unq_ESD = ESD;

% trophic strategy
eco_pars.trophic = linspace(1,0,eco_pars.ntroph);

[~,ESD]=ndgrid(eco_pars.trophic,ESD);
eco_pars.ESD=reshape(ESD,numel(ESD),1);

% grid traits and reshape to one dimension
[trophic Vol]=ndgrid(eco_pars.trophic,eco_pars.V);
eco_pars.trophic=reshape(trophic,numel(trophic),1);
eco_pars.V=reshape(Vol,numel(Vol),1);



% Trait diffusion/mutation matrix
[xn yn]=ndgrid(1:eco_pars.ntroph,1:eco_pars.nsize);
xl=reshape(xn,numel(xn),1);
yl=reshape(yn,numel(yn),1);
tvect=[xl yl];


eco_pars.phymin = 0;

nn=eco_pars.jpmax;
tm=zeros(nn);
for i=1:nn
    for j=1:nn
        adjvect=abs(tvect(i,:)-tvect(j,:));
        % diffusion in first dimension (trophic strategy)
        if adjvect(1)==1 & adjvect(2)==0
            % normalise flux to range of trait space
            tm(i,j)=range(eco_pars.trophic)^2 .* eco_pars.mutrat(1)./(eco_pars.trophic(i)-eco_pars.trophic(j)).^2;
            %                 tm(i,j)=eco_pars.mutrat(1)./(1./(eco_pars.ntroph-1)).^2;
        end
        % diffusion in second dimension (size)
        if adjvect(2)==1 & adjvect(1)==0
            % normalise flux to range of trait space
            tm(i,j)=range(log10(eco_pars.V))^2 .* eco_pars.mutrat(2)./(log10(eco_pars.V(i))-log10(eco_pars.V(j))).^2;
            %                 tm(i,j)=eco_pars.mutrat(1)./(1./(eco_pars.nsize-1)).^2;
        end
    end
end

% adjust threshold by number of bins (more bins will have less biomass in each)
eco_pars.extnct = eco_pars.extnct ./ eco_pars.jpmax;

tm(eye(nn)==1)=1-sum(tm,2);

if true % mutation
    eco_pars.Pmut=sparse(tm);
else
    eco_pars.Pmut=speye(size(tm));
end

clear trophic Vol xl yl xn yn ESD ESR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


eco_pars.muinf   = eco_pars.muinf_a .* eco_pars.V .^ eco_pars.muinf_b;
eco_pars.Qmin    = eco_pars.Qmin_a  .* eco_pars.V .^ eco_pars.Qmin_b;
eco_pars.Vmax    = eco_pars.Vmax_a  .* eco_pars.V .^ eco_pars.Vmax_b;
% Edwards et al. L&O 2012
eco_pars.alpha   = eco_pars.alpha_a .* eco_pars.V .^ eco_pars.alpha_b;
eco_pars.g       = eco_pars.g_a     .* eco_pars.V .^ eco_pars.g_b;

eco_pars.allometric.muinf=[eco_pars.muinf_a eco_pars.muinf_b];
eco_pars.allometric.Qmin =[eco_pars.Qmin_a  eco_pars.Qmin_b ];
eco_pars.allometric.Vmax =[eco_pars.Vmax_a  eco_pars.Vmax_b ];
eco_pars.allometric.alpha=[eco_pars.alpha_a eco_pars.alpha_b];
eco_pars.allometric.g    =[eco_pars.g_a     eco_pars.g_b    ];
eco_pars=rmfield(eco_pars,{'muinf_a','muinf_b','Qmin_a','Qmin_b','Vmax_a','Vmax_b','alpha_a','alpha_b','g_a','g_b'});

eco_pars.mumax=eco_pars.muinf .* eco_pars.Vmax ...
    ./ (eco_pars.muinf .* eco_pars.Qmin + eco_pars.Vmax ); % d^{-1}

% Grazing matrix
herbmat = zeros(eco_pars.jpmax);
[Vpry Vprd]= meshgrid(eco_pars.V,eco_pars.V);
prd_pry=Vprd./Vpry;
herbmat=exp(-(log(prd_pry./eco_pars.pp_opt)).^2 ./ (2.*eco_pars.pp_sig.^2));
herbmat(herbmat<1e-3)=0;
eco_pars.prdpry=sparse(herbmat);
eco_pars.pryprd=eco_pars.prdpry';

clear herbmat Vpry Vprd prd_pry

if eco_pars.ntroph>1
    % apply trophic trade-offs
    eco_pars.mumax = eco_pars.mumax .*    eco_pars.trophic .^eco_pars.tau;
    eco_pars.g     = eco_pars.g     .* (1-eco_pars.trophic).^eco_pars.tau;
end



%% % NONNEGATIVE OPTION SEEMS TO VIOLATE MASS CONSERVATION

% set solver options: positive definite, tolerances, min/max tsteps
AbsTol      = 1e-3;
RelTol      = 1e-6;
MaxStep     = inf;
InitialStep = 1;
Refine      = 1; 

eco_pars.odeoptions=odeset('AbsTol',AbsTol,...
                           'RelTol',RelTol,...
                           'MaxStep',MaxStep,...
                           'InitialStep',InitialStep,...
                           'Refine',Refine);











