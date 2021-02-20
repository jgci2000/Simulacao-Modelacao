%Exact elliptical orbit plotter for structures of masses used as an input to gravity_sim.m

function elliptical_orbit_plotter

%% Fixed inputs %%%

%Number of polar values in each orbit
N = 1000;

%

%Load gravity_sim.m settings file
[filename, pathname] = uigetfile('..\settings\*.m', 'Choose gravity sim settings file');
if ~isequal(filename,0) && ~isequal(pathname,0)
    
    %Load data if file has been selected
    [Lmax,masses,view_option] = load_data( [pathname,filename] );
    
    %Create figure and axis
    figure('color',[1 1 1],'name','elliptical orbit plotter','units','normalized');
    axes('units','normalized','outerposition',[0,0,1,1]);
    grid on
    xlabel('x')
    ylabel('y')
    zlabel('z')
    
    %Generate color selection for orbits from current colormap
    map = colormap;
    dim = size(map);
    colour_i = round(linspace(1,dim(1),length(masses)));
    
    %Initialize length extent of plot
    smax = 1;
    
    %Generate and plot orbits
    for m=1:length(masses)
        
        %Possibly enlarge length extent of plot
        smax = max( [smax, masses(m).init_orbit_s] );
        
        %Generate orbital ellipse
        [ux,vx1,vx2,...
            uy,vy1,vy2,...
            uz,vz1,vz2,...
            xc,x1,x2,...
            yc,y1,y2,...
            zc,z1,z2,...
            theta,theta_dot,t,j,J,E,P] =...
            two_body_elliptical_orbit(  masses(m).init_orbit_G, masses(m).init_orbit_M, masses(m).mass, masses(m).init_orbit_s,...
            masses(m).init_orbit_ecc, masses(m).init_orbit_theta0, masses(m).init_orbit_h0,...
            masses(m).init_orbit_hdot, masses(m).init_orbit_d, masses(m).init_orbit_alpha,N );
        
        %Set legend entry from mass name
        legend_str{m} = masses(m).name;
        
        %Plot orbits
        p(m) = plot3(x2,y2,z2,'r'); hold on;
        if isempty(masses(m).marker_RGB)==1
            colour = map(colour_i(m),:);       
        else
            colour = masses(m).marker_RGB;
        end
        set(p(m),'color',colour);
        
        %Plot start and finish values, and orbit focus
        plot3(x2(1),y2(1),z2(1),'o','color',colour);
        plot3(xc(1),yc(1),zc(1),'+','color',colour);
        plot3(x2(N),y2(N),z2(N),'*','color',colour);
        plot3(xc(N),yc(N),zc(N),'+','color',colour);
    end
    xlim([-Lmax/1.5,Lmax/1.5]);
    ylim([-Lmax/1.5,Lmax/1.5]);
    zlim([-Lmax/1.5,Lmax/1.5]);
    axis vis3d;
    grid on
    view(view_option);
    if view_option==2
        set(gca,'position',[0,0,1,1])
    else
        set(gca,'position',[0.0,0.13,1,0.7])
    end    
    xlabel('x')
    ylabel('y')
    zlabel('z')
    legend(p,legend_str,'location','eastoutside');
    title( {['Elliptical orbit plotter: ',num2str(length(masses)),' masses'],...
        strrep(filename,'_','\_')} );
    print( gcf, '-dpng','-r300',['Elliptical orbit plot',...
        '. ',strrep(datestr(now),':','-'),'.png'])
end

%%

%Get settings from file
function [Lmax,masses,view_option]  = load_data( data_name )

%Add all functions in this directory
addpath(pwd);

%Note the desired outputs will be defined in the data file, which is a
%MATLAB script. All other parameters will be ignored since the file is run
%within this function.
run(data_name);

%Remove functions from path
rmpath(pwd);

%End of code