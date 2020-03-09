clear

fdir = 'Run_N3_Kappa_0.01_2020_03_09';
load(['../Output/' fdir '/Workspace_dump.mat'])
ESD=nthroot(6.*eco_pars.V./pi,3);


data1=data(53:53:end,:);

addpath('..')
[ rgb, rgb3 ] = cmap_2d(eco_pars);
%%

tr_index=zeros(size(data1,1),size(data1,2)-1);

for t=1:size(data1,1)
    y=data1(t,:)';

    [DV, invfit] = ecosystem(t,y,eco_pars,dead,minphy);
    
    z = invfit.autotrophic./(invfit.autotrophic+invfit.grazing_gain);
    z(isnan(z))=0;
    
    tr_index(t,:)=z;
end
% Autotroph   --> tr_index=1
% Heterotroph --> tr_index=0

%%


Nit=data1(:,1);
Bio=data1(:,2:end);
Bio(Bio<0)=0;

figure(1)
clf
c=0;
for i=unique(round(logspace(0,4,200)))
    
    N=Nit(i);
    X=Bio(i,:)';
    X(X<1e-9)=0;
    
    clf
    
    
    M=(X.*eco_pars.g)'.*eco_pars.pryprd.*X;
    P=X .* eco_pars.mumax.*eco_pars.alpha.*N ./ (eco_pars.mumax + eco_pars.alpha.*N);
    
    M=[M;P'];
    M(:,end+1)=0;
    M(M<1e-6)=0;
    
    G=digraph(M,'OmitSelfLoops');
    
    if size(G.Edges,1)>0
        LWidths = 5.*(G.Edges.Weight/max(G.Edges.Weight)).^2;

        p=plot(G);
        
        p.XData=[1-eco_pars.trophic;0.25];
        p.YData=[ESD;0.06];
        p.MarkerSize=eps;
        p.LineWidth=LWidths;
        p.EdgeAlpha=0.1;
        p.EdgeColor='k';
    end
    
    hold on
    scatter(1-eco_pars.trophic,ESD,300.*sqrt(X)+eps,rgb,'filled')
    scatter(1-eco_pars.trophic,ESD,300.*sqrt(X)+eps,'k')
    
    scatter(0.25,0.06,300.*sqrt(N)+eps,[0.25 0.25 0.25],'filled')
    scatter(0.25,0.06,300.*sqrt(N)+eps,'k')
    
    
    set(gca,'YScale','log','YTick',[0.2 2 20 200 2000]);
    xlim([0 1])
    ylim([0.06 6000])
    
    xlabel('Trophic Strategy')
    ylabel(['Diameter (' char(181) 'm)'])
    
    title(i)
    drawnow
    set(gcf,'Color','w')
    box on
    axis square
    set(gca,'fontsize', 18)
    
    c=c+1;
    sname=['../Figures/FoodWeb_' num2str(c,'%05i') '.png'];
    export_fig(sname,'-r300')
end



