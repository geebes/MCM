function [DV invfit]=ecosystem(t,v0,p,dead,minphy);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % time dependent limitation 
    if p.seasonalcycle
        limitation = interp1(0:365, p.gamma_cycle,rem(t,365));
    else
        limitation = p.gamma;
    end
    
    % Extract state variables
    N = max(v0(1    ,:),0); % [1 x ns]
    P = max(v0(2:end,:),0); % [jp x ns]
    
    dead    = find(P < p.extnct);
    P(dead) = minphy;
    assignin('caller','dead',dead);
    evalin('caller','y(1+dead) = minphy;');
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Allee effect - brings death to populations below threshold biomass
    AlleeDeath  =  10.* exp(-P./p.extnct); % Allee =  1e3.*exp(-P./p.extnct);
    AlleeGrowth =  spfun(@(x) 1-exp(-x./p.extnct),P); % Allee =  1-exp(-P./p.extnct);
    % N.B. spfun selectively applies a function to only the nonzero elements
    % of a sparse matrix S, preserving the sparsity pattern of the original matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % NUTRIENT UPTAKE (i.e. AUTOTROOPHIC GROWTH)
    mumax    = p.mumax;
    alphaN   = p.alpha .* N;   % Implicit Expansion of row and column - expands input vectors for matrix of all possible combinations [jp x ns]
    Nuptk    = limitation .* mumax .* alphaN ./ (mumax + alphaN);
    Nuptk(mumax + alphaN<=0) = 0; % Check for possible divide by zero
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % GRAZING
    g = p.g; % [jpred x ns] % GRAZING ATTACK RATE % Implicit Expansion of row and column - expands input vectors for matrix of all possible combinations
        
    % PREDATOR GAINS
        % TOTAL PREY AVAILABLE TO EACH (Potential) PREDATOR
        pry_P    = p.prdpry * P; % [jpred x jprey] *[jprey x ns] = [jpred x ns]
        refg     = spfun(@(x) 1-exp(p.refuge.*x),pry_P); % [jpred x ns] = [jpred x ns] % PREDATOR EFFORT (PREY REFUGE)
        g_refg   = refg .* g;
        % PREDATOR ATTACK RATES ON WHOLE COMMUNITY (pre-assimilation)
        GGain  = g_refg .* pry_P;                          % [jpred x ns   ].*[jpred x ns] = [jpred x ns]  
    % PREY LOSSES  
        % PREDATOR POPULATIONS' FEEDING EFFORT   
        g_att  = AlleeGrowth .* P.*g_refg;                      % [jpred x ns].*[jpred x ns] = [jpred x ns]
        % COMMUNITY ATTACKS FELT BY EACH PREY   
        GLoss  = p.pryprd*g_att;          % [jprey x jpred]*[jpred x ns] = [jprey x ns]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % MORTALITY
    BaseMort = p.linearmort .* (1 + P) + p.kappa;
    AlleeMort = BaseMort + AlleeDeath;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % GROSS GROWTH RATE
    ggr = Nuptk + GGain.*p.lambda;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % REPRODUCTIVE FLUX (including 'mutations')
    mutflux = p.Pmut * (AlleeGrowth .* ggr .* P); % GGRs subject to mutation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ODEs
    dNdt = - sum(Nuptk.*P,1) ...
           + p.kappa .* (p.Nsupply - N);
          
    dPdt = + mutflux   ...           % ggr including mutational flux
         - ( AlleeMort + GLoss).*P;    % linear mortality + grazing losses

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    invfit.invfit = (+ ggr   ...              % gross growth rate
                  - BaseMort ...           % linear mortality
                  - GLoss);                % Grazing losses

    invfit.autotrophic  = Nuptk;
    invfit.grazing_gain = GGain.*p.lambda;
    invfit.grazing_loss = GLoss;
    invfit.linear_mort  = BaseMort;
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Output
    DV = [dNdt;dPdt];
end