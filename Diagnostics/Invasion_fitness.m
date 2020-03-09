clear

load('/Users/baw103/GitHub/MCM/output/run_N5_constant_newseed/run_N5_constant_newseed.mat')
ESD=nthroot(6.*eco_pars.V./pi,3);


data1=data(53:53:end,:);

[ rgb, rgb3 ] = cmap_2d(eco_pars);
%%

tr_index = zeros(size(data1,1),size(data1,2)-1);
aut      = zeros(size(data1,1),size(data1,2)-1);
het      = zeros(size(data1,1),size(data1,2)-1);

for t=1:size(data1,1)
    y=data1(t,:)';

    [DV, invfit] = ecosystem(t,y,eco_pars,dead,minphy);
    
    z = invfit.autotrophic./(invfit.autotrophic+invfit.grazing_gain);
    z(isnan(z))=0;
        
    tr_index(t,:)=z;
    
    aut(t,:) = invfit.autotrophic;
    het(t,:) = invfit.grazing_gain;
end
% Autotroph   --> tr_index=1
% Heterotroph --> tr_index=0

%%
[ rgb, rgb3 ] = cmap_2d(eco_pars);

isnotnano=find(ESD<2 | ESD>20);

Bio=data1(:,2:end)./eco_pars.Qmin'; % convert to abundance
Bio(Bio<0)=0;
BioN=Bio;
BioN(:,isnotnano)=0;

[i,j] = find(BioN>1e-3);
j=unique(j);
numel(j)

Bio1=BioN(:,j); 
tr_index1=tr_index(:,j);
sz_index1=eco_pars.V(j);
[~,~,sz_index1] = unique(sz_index1);

trophic=eco_pars.trophic(j);
trophic(trophic<1)=0;

hetfrac = sum(het(:,j).*Bio1,2) ./ (sum(aut(:,j).*Bio1,2) + sum(het(:,j).*Bio1,2));

return

%%
figure(1)
clf

sh1=subplot(211);
set(gca,'XScale','lin')
hold on
Biosum1=cumsum(Bio1,2)./sum(Bio1,2);

xx = [1:size(Biosum1,1) size(Biosum1,1):-1:1];
for i=fliplr(2:size(Biosum1,2)) % loop through populations
    z=trophic(i);
    a=1-sz_index1(i)./eco_pars.nsize;
    
    cdata=[1-z z 1-z].*a;
    if i>2
        yy = [Biosum1(:,i)' fliplr(Biosum1(:,i-1)')];
    else % add zero row for last size lass
        yy = [Biosum1(:,i)' zeros(1,size(Biosum1,1))];
    end
    if sz_index1(i)~=sz_index1(i-1) % block black lines between size classes
        h1 = fill(xx,yy,cdata,'EdgeColor','k','LineW',1);
    else
        h2 = fill(xx,yy,cdata,'EdgeColor','none');
        uistack(h2,'bottom')
    end
    drawnow
end
plot(zeros(size(Biosum1,1),1),'Color','k','LineWidth',1);
plot(ones(size(Biosum1,1),1),'Color','k','LineWidth',1);
xlim([0 1000])
    

ch=colorbar;
colormap(sh1,[0 1 0;1 0 1]);
ch.Ticks = [0 1];
ch.TickLabels={'Aflagellated','Flagellated'}
ylabel('Biomass (mmol N m^{-3})')



sh2=subplot(212);
set(gca,'XScale','lin')
hold on
Biosum1=cumsum(Bio1,2)./sum(Bio1,2);

for i=fliplr(1:size(Biosum1,2)) % loop through populations
    
    z=ceil(tr_index1(:,i).*255);
%     z(z<255)=0;
    cdata=uint8([255-z z 255-z z.*0+1])';
    
	h = plot(Biosum1(:,i)','k','LineW',1);
    drawnow
    set(h.Edge, 'ColorBinding','interpolated', 'ColorData',cdata)
end
ch=colorbar;
colormap(sh2,greenmag(64));
ch.Ticks = [0 1];
ch.TickLabels={'Autotrophy','Heterotrophy'}
xlabel('Time (years)')
ylabel('Biomass (mmol N m^{-3})')
xlim([0 1000])
%%

fig2=figure(2);
clf
set(0,'defaultAxesFontSize',16)

sh1=subplot(141);
set(gca,'XScale','lin','XTick',0:0.5:1,'XTickLabel',{'0','50','100'})
hold on
Biofrac1 = Bio1./sum(Bio1,2);

xx = 1:size(Biosum1,1);
for i=1:size(Biosum1,2) % loop through populations
    z=trophic(i);
    
    if z==1
        cdata='g';
    else
        cdata='b';
    end
    plot(Biofrac1(:,i)',xx,'Color',cdata,'LineW',1);
    drawnow
end
ylim([0 500])
box on
ch=colorbar('SouthOutside');
colormap(sh1,[0 1 0;0 0 1]);
ch.Ticks = [0 1];
ch.TickLabels={'Aflagellate','Flagellate'}
xlabel('Biomass fraction (%)')
ylabel('Time (years)')


sh2=subplot(142);
set(gca,'XScale','lin','XTick',0:0.5:1,'XTickLabel',{'0','50','100'})

for i=1:size(Biosum1,2) % loop through populations
    z=ceil(tr_index1(:,i).*255);
    cdata=uint8([255-z z 255-z ones(size(z)).*255])';
    
	h = plot(Biofrac1(:,i)',xx,'k','LineW',1);
    drawnow
    set(h.Edge, 'ColorBinding','interpolated', 'ColorData',cdata)
    hold on
end
ylim([0 500])
box on
ch=colorbar('SouthOutside');
colormap(sh2,[linspace(0,1,64)' linspace(1,0,64)' linspace(0,1,64)']);
ch.Ticks = [0 1];
ch.TickLabels={'Autotrophy','Heterotrophy'}
xlabel('Biomass fraction (%)')

sh3=subplot(143)
set(gca,'XScale','lin','XTick',0:0.5:1,'XTickLabel',{'0','50','100'})
hold on
isHNP=find(trophic<1); % flagellated nanoplankton
plot(sum(Bio1(:,isHNP),2)./sum(Bio1,2),1:size(Bio1,1),'k','LineW',1)
plot(hetfrac,1:size(Bio1,1),'k--','LineW',1)
ylim([0 500])
box on
xlabel('NP ''heterotroph'' fraction (%)')
ch=colorbar('SouthOutside');

subplot(144)
plot(sum(Bio>0,2),1:size(Bio,1),'k','LineW',1)
ylim([0 500])
xlabel('''species'' richness')
ch=colorbar('SouthOutside');

sname=['timeline_abundance.eps'];
fig2.Renderer='painters';
print(gcf, sname, '-depsc', '-r0');


return
%%
fh=figure(3);
fh.Position=[47 1 821 1344];
c=0;
for t=unique(round(logspace(0,4,200)))
    y=data1(t,:)';

    [DV, invfit] = ecosystem(t,y,eco_pars,dead,minphy);
    
    zclr = invfit.autotrophic./(invfit.autotrophic+invfit.grazing_gain);
    zclr(isnan(zclr))=0;

    invfit.invfit       =reshape(invfit.invfit,eco_pars.nsize,eco_pars.ntroph)';
    invfit.autotrophic  =reshape(invfit.autotrophic,eco_pars.nsize,eco_pars.ntroph)';
    invfit.grazing_gain =reshape(invfit.grazing_gain,eco_pars.nsize,eco_pars.ntroph)';
    invfit.grazing_loss =reshape(invfit.grazing_loss,eco_pars.nsize,eco_pars.ntroph)';
    invfit.linear_mort	=reshape(invfit.linear_mort,eco_pars.nsize,eco_pars.ntroph)';
    
    zt=invfit.invfit;
    za=invfit.autotrophic - invfit.grazing_loss - invfit.linear_mort;
    zh=invfit.grazing_gain - invfit.grazing_loss - invfit.linear_mort;
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % full invasion fitness
    figure(3)
    clf

    cxlim=0.2;
    nclr=101;
    nroot=3;
    szscl=100.*sign(Bio(t,:)').*abs(nthroot(Bio(t,:)',nroot))+eps;
    ctcks=linspace(-cxlim,cxlim,nclr);
    colormap(flipud(redblue(nclr)));
    
    contourf(linspace(0,1,eco_pars.ntroph),log10(eco_pars.unq_ESD)',zt,ctcks,'LineC','none')
    set(gca,'XTick',[0 0.5 1],'XTickLabel',{'Autotroph','Mixotroph','Heterotroph'})
    set(gca,'YTick',log10(2.*logspace(-1,3,5)),'YTickLabel',{'0.2','2','20','200','2000'})
    hold on
    scatter(1-eco_pars.trophic,log10(ESD),szscl,[1-zclr zclr 1-zclr],'k')
    scatter(1-eco_pars.trophic,log10(ESD),szscl,[1-zclr zclr 1-zclr],'filled','MarkerFaceAlpha',0.5)
    axis xy
    caxis([-cxlim cxlim])
    axis square
    title('Invasion fitness')
    ylabel('ESD (microns)')
    text(0.06,3.5,['Year ' num2str(t)],'FontSize',16)
    hclr=colorbar('Location','SouthOutside');
    
    spcb=axes;
    s2=colorbar('Location','SouthOutside','Position',hclr.Position-[0 0.03 0 0]);
    colormap(s2,greenmag);
    axis off
    caxis([0 1])
    s2.Ticks=[0 0.5 1];
    s2.TickLabels={'Autotrophic','Mixotrophic','Heterotrophic'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    drawnow
    c=c+1;
    sname=['Figures/Fitness_' num2str(c,'%03i') '.png'];
    set(gcf,'Color','w')
    export_fig(sname,'-r300')
    

end
%%

fh=figure(4);
clf
fh.Position=[326 1 797 1344];

c=0;
for t=[1 2 5 10 20 50 100 200 500 1000 2000 5000]
    y=data1(t,:)';

    [DV, invfit] = ecosystem(t,y,eco_pars,dead,minphy);
    
    zclr = invfit.autotrophic./(invfit.autotrophic+invfit.grazing_gain);
    zclr(isnan(zclr))=0;

    invfit.invfit       =reshape(invfit.invfit,eco_pars.nsize,eco_pars.ntroph)';
    invfit.autotrophic  =reshape(invfit.autotrophic,eco_pars.nsize,eco_pars.ntroph)';
    invfit.grazing_gain =reshape(invfit.grazing_gain,eco_pars.nsize,eco_pars.ntroph)';
    invfit.grazing_loss =reshape(invfit.grazing_loss,eco_pars.nsize,eco_pars.ntroph)';
    invfit.linear_mort	=reshape(invfit.linear_mort,eco_pars.nsize,eco_pars.ntroph)';
    
    zt=invfit.invfit;
    za=invfit.autotrophic - invfit.grazing_loss - invfit.linear_mort;
    zh=invfit.grazing_gain - invfit.grazing_loss - invfit.linear_mort;
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % full invasion fitness
    c=c+1;
    subplot(4,3,c)

    cxlim=0.2;
    nclr=101;
    nroot=3;
    Bio(Bio<eco_pars.extnct)=0;
    szscl=50.*sign(Bio(t,:)').*abs(nthroot(Bio(t,:)',nroot))+eps;
    ctcks=linspace(-cxlim,cxlim,nclr);
    colormap(flipud(redblue(nclr)));
    
    zt(zt>+cxlim)=+cxlim;
    zt(zt<-cxlim)=-cxlim;
    contourf(linspace(0,1,eco_pars.ntroph),log10(eco_pars.unq_ESD)',zt,ctcks,'LineC','none')
    set(gca,'XTick',[])
    set(gca,'YTick',[])
 
    hold on
    scatter(1-eco_pars.trophic,log10(ESD),szscl,[1-zclr zclr 1-zclr],'k')
    scatter(1-eco_pars.trophic,log10(ESD),szscl,[1-zclr zclr 1-zclr],'filled','MarkerFaceAlpha',0.5)
    axis xy
    caxis([-cxlim cxlim])
    axis square
    
    set(gca,'XTick',[0 0.5 1])
    set(gca,'YTick',log10(2.*logspace(-1,3,5)))
    if c==1 | c==4 | c==7 | c==10
        ylabel('ESD')
         set(gca,'YTickLabel',{'0.2','2','20','200','2000'})
    else
         set(gca,'YTickLabel',{''})
    end
    if c>9
        xlabel('Trophic strategy')        
        set(gca,'XTickLabel',{'A','M','H'})
    else
        set(gca,'XTickLabel',{''})
    end
    
    text(0.06,3.4,['Year ' num2str(t)],'FontSize',16)

    drawnow
end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ax=axes;
ax.Position=[0.1300 0.06 0.7750 0.8150];
caxis([-cxlim cxlim])
colormap(flipud(redblue))
hclr=colorbar('Location','SouthOutside');
axis off

ax=axes;
ax.Position=[0.1300 0.06 0.7750 0.8150];
s2=colorbar('Location','SouthOutside','Position',hclr.Position-[0 0.04 0 0]);
colormap(s2,greenmag);
axis off
caxis([0 1])
s2.Ticks=[0 0.5 1];
s2.TickLabels={'Autotrophic','Mixotrophic','Heterotrophic'};


sname=['Figures/Fitness_timeline.eps'];
set(gcf,'Color','w')
saveas(gcf,sname,'epsc')
    
