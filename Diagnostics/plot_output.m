
function plot_output(tout,yout,eco_pars,disp_indx)

    % discrete trait values
    t_troph = linspace(0,1,eco_pars.ntroph);
    t_size  = eco_pars.unq_ESD;
    
    % yticklabels
    ytickL=repmat({''},1,eco_pars.ntroph);
    ytickL{1}='Phytoplankton';
    ytickL{end}='Zooplankton';
    ytickL{ceil(eco_pars.ntroph/2)}='Mixotroph';
    for i=1:eco_pars.nsize
        if ismember(round(eco_pars.unq_ESD(i).*1e5),logspace(3,8,6).*6)
            xtickL{i}=num2str(t_size(i),'%4.2f');
        else
            xtickL{i}='';
        end
    end
    ESD=nthroot(6.*eco_pars.V./pi,3);

    az = -120;
    el = 27;

    for i=disp_indx
        % diagnose invasion fitness
        dead=[];
        minphy=eco_pars.minphy;
        [~,invfit]=ecosystem(1,yout(i,:)',eco_pars,dead,minphy);
        invfit=reshape(invfit.invfit,eco_pars.ntroph,eco_pars.nsize);
        plot_threshold = 0.25;
        isinZ=abs(invfit)>plot_threshold;
        invfit(isinZ) = sign(invfit(isinZ)) .* plot_threshold;

        clf
        axes
        view(az,el);
        ax1=gca;
        ax1.Colormap=winter;
        ax1.Color = 'none';      
        ax1.XTick=1:eco_pars.ntroph;
        ax1.XTickLabel=xtickL;
        ax1.YDir = 'reverse';  
        ax1.YTick=1:eco_pars.nsize;
        ax1.YTickLabel=ytickL;
        ax1.ZAxis.Color='none';
        hold on
        
        matrix=reshape(yout(i,2:end),eco_pars.ntroph,eco_pars.nsize);
        matrix(matrix<=eco_pars.extnct)=eco_pars.minphy;
        matrix(matrix<=0)=NaN;
        logmatrix=log10(matrix);
        
        hbar = bar3(matrix,1);
                
        hold on;
        axis tight
        for ii = 1:eco_pars.ntroph
            c     = get(hbar(ii), 'CData');
            color = repelem(repmat(logmatrix(:, ii), 1, 4), 6, 1);
            set(hbar(ii), 'CData', color);
        end
        caxis(log10([eco_pars.extnct 1]))
        zlim([0 1])
        drawnow
        
        xlabel(['Size (' char(181) 'm)'])
        ylabel('Trophic strategy')
        zlabel('Invasion fitness (d^{-1})')
        
        hc = colorbar('location','WestOutside');
        hc.Position = [0.1 0.1 0.01 0.15];
        hc.TickLabels = 10.^hc.Ticks;
        ht = text(-0.03,-0.05,{'Biomass';'(mmol m^{-3})'},...
                  'Units','normalized',...
                  'HorizontalAlignment','center');
        axes
        view(az,el);
        ax2=gca;
        ax2.Position=ax1.Position;
        ax2.XAxis.Color='none';
        ax2.XScale = 'Log';
        ax2.XAxis.Color='none';
        ax2.YAxis.Color='none';
        ax2.ZAxis.FirstCrossoverValue  = eco_pars.unq_ESD(1);
        ax2.ZAxis.SecondCrossoverValue  = 0;
        ax2.ZTick  = linspace(-plot_threshold,plot_threshold,5);
        ax2.Color = 'none';
        ax2.Colormap=flipud(redblue);
        hold on
        
        h=surf(t_size,1-t_troph,invfit);
        axis tight
        zlim([-2 1].*plot_threshold)
        caxis([-1 1].*plot_threshold)
        scatter3(ones(1,1000).*t_size(1),...
                 ones(1,1000).*t_troph(1),...
                 linspace(-plot_threshold,plot_threshold,1000),...
                 5,...
                 linspace(-plot_threshold,plot_threshold,1000),...
                 'filled')
        
        isextant=find(yout(i,2:end)>eco_pars.extnct);
        scatter3(ESD(isextant),eco_pars.trophic(isextant),invfit(isextant),2e2.*sqrt(yout(i,1+isextant)),'k','filled');
                
        ax=axis;
        ph=patch([ax(1) ax(1) ax(2) ax(2)],[ax(3) ax(4) ax(4) ax(3)],'k');
        ph.FaceAlpha=0.25;
        
        ty=floor(tout(i)./365);
        if ty==1
            ht = title(['Time = ' num2str(ty) ' year']);
        else
            ht = title(['Time = ' num2str(ty) ' years']);
        end
        ht.Units = 'normalized';
        ht.Position = [0 0.75 0];
        ht.HorizontalAlignment = 'left';
        
        axes
        view(az,el);
        ax3=gca;
        ax3.XLim=ax2.XLim;
        ax3.YLim=ax2.YLim;
        ax3.ZLim=ax2.ZLim;
        hold on
        line(ones(1,2).*t_size(1),...
             ones(1,2).*t_troph(1),...
             linspace(-2.*plot_threshold,-plot_threshold,2),'color','w','LineW',1)
        axis off
        
        drawnow
        

    end
end