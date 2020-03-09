function [ rgb, rgb3 ] = cmap_2d(eco_pars)


    ns=eco_pars.nsize;
    nt=eco_pars.ntroph;
    
    % define colormap
    x=[0 0 0.5 1 1]'; % x cordinates of corners and centre
    y=[0 1 0.5 0 1]'; % y cordinates of corners and centre
    %  g y w   b m
    r=[0 1 1   0 1]';   % corresponding r values
    g=[1 1 1   0 0]';   % corresponding g values
    b=[0 0 1   1 1]';   % corresponding b values

    
    xl=linspace(0,1,nt); % unique values of trait x
    yl=linspace(0,1,ns); % unique values of trait y
    
    [xx yy]=meshgrid(xl,yl); % generate trait grid
    
    Fr=scatteredInterpolant(x,y,r); % r values across grid
    Fg=scatteredInterpolant(x,y,g); % g values across grid
    Fb=scatteredInterpolant(x,y,b); % b values across grid
    
    r_out = Fr(xx,yy); r_out(r_out>1)=1;r_out(r_out<0)=0;
    g_out = Fg(xx,yy); g_out(g_out>1)=1;g_out(g_out<0)=0;
    b_out = Fb(xx,yy); b_out(b_out>1)=1;b_out(b_out<0)=0;
    
    rgb=[r_out(:) g_out(:) b_out(:)]; % rgb array

    rgb3=cat(3,r_out,g_out,b_out);
end