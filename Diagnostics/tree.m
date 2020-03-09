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

isnotnano=find(ESD<2 | ESD>20);

Bio=data1(:,2:end);
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

%%
fh=figure(4);
clf
set(0,'defaultAxesFontSize',16)


subplot(1,4,[1 3])
% plot first 500000 'generations' (i.e. days)
data=data1(1:round(0.5e6/365),2:end)';
data=full(data);
data=reshape(data,51,51,size(data,2));
data=padarray(data,[1 1 1],0,'both');
% data=sqrt(data);

[xx yy zz] = meshgrid(0:52,0:52,1:size(data,3));

p=patch(isosurface(xx,yy,zz,data,1e-3));
hold on
image(linspace(0.5,51.5,51),linspace(0.5,51.5,51),rgb3)
set(gca,'YDir','Reverse')
xlim([0 51])
ylim([0 51])
zlim([0 size(data,3)])
view(-60,15);
box on

% get trait indices at vertices
c1=ceil(p.Vertices(:,1));
c2=ceil(p.Vertices(:,2));
c1(c1==0)=1;
c2(c2==0)=1;
ci=(c1-1).*51+c2;
    
p.FaceVertexCData=rgb(ci,:);
    
p.FaceColor = 'interp';
p.EdgeColor = 'none';

ax = gca;
% ax.BoxStyle = 'full';

camlight HEADLIGHT
camproj('perspective')
lighting gouraud
material DULL

xb=[  11  27]; % nano size range
yb=[  0  12];  % autotrophic end of spectrum
zb=[  0 750];  % first few years
hold on
hp=patch([xb(1) xb(1) xb(1) xb(1)],[yb(1) yb(2) yb(2) yb(1)],[zb(1) zb(1) zb(2) zb(2)],'k');
hp.FaceAlpha=0.25;
hp=patch([xb(1) xb(1) xb(2) xb(2)],[yb(1) yb(2) yb(2) yb(1)],[zb(2) zb(2) zb(2) zb(2)],'k');
hp.FaceAlpha=0.25;
hp=patch([xb(1) xb(1) xb(2) xb(2)],[yb(2) yb(2) yb(2) yb(2)],[zb(1) zb(2) zb(2) zb(1)],'k');
hp.FaceAlpha=0.25;

ax.XTick=1:10:51;
ax.XTickLabel={'0.06','0.6','6','60','600','6000'};
xlabel(['Size (' char(181) 'm)'])
ax.YTick=[1 6 11 16 21 26 31 36 41 46 51];
ax.YTickLabel={'Autotrophic','','','','','Mixotrophic','','','','','Heterotrophic'};
ylabel(['Trophic strategy'])
ax.ZTick=[];
zlabel(['Time \rightarrow'])

asp=pbaspect;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inset
ax=axes;
ax.Position=[0.625 0.235 0.35 0.65];
p=patch(isosurface(xx,yy,zz,data,1e-3));
    
p.FaceVertexCData=rgb(ci,:);
    
p.FaceColor = 'interp';
p.EdgeColor = 'none';

ax = gca;
% ax.BoxStyle = 'full';

camlight HEADLIGHT
lighting gouraud
material DULL
camproj('perspective')
hold on
image(linspace(0.5,51.5,51),linspace(0.5,51.5,51),rgb3)
set(gca,'YDir','Reverse')
xlim([xb])
ylim([yb])
zlim([zb])
view(-64,15);
box on

ax.XTick=1:10:51;
ax.XTickLabel={''};
ax.YTick=[1 6 11 16 21 26 31 36 41 46 51];
ax.YTickLabel={''};
ax.ZTick=[];


set(gcf,'Color','w')

    
sname=['Figures/Tree.png'];
export_fig(sname,'-r500')
