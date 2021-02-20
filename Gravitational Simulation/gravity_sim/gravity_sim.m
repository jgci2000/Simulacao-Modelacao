%gravity_sim
% Gravity simulator, where the user controls the position of a star with the
% keyboard and/or mouse. A star is surrounded by rings of effectively
% massless planets. These are placed in initially circular orbits about the
% star. Any number of stars can be defined.
%
% LAST UPDATED by Andy French. April 2013

function gravity_sim

%Add 'functions' directory to the MATLAB path. This contains supplementary
%functions necessary to define the inputs to the gravity simulator.
addpath([pwd,'\functions'],'-begin')

%Initialise structure d into which which all gravity_sim data is stored.
%This is a global variable so all sub-functions of gravity_sim.m will be
%able to see and modify it.
global d

%

%% Hard coded inputs %%

%Position change resulting from arrow key press
d.ds = 0.1;

%Mass change resulting from key press (in Solar masses)
d.dM = 0.1;

%Radius for number density
d.rho_r = 1;

%Density multiplier. Maximum colour scale for density map is set by DM*
%number of objects / area of space (Lmax^2 * 4)
d.DM = 10;

%Speed multiplier. Maximum colour scale for speed map is set by SM*
%maximum initial speed * number of objects / area of space (Lmax^2 * 4)
d.SM = 10;

%Number density map grid
d.N = 20;

%Axes positions depending on 2D or 3D view. (Normalized units)
d.view2_axes_position = [0.1,0.08,0.8,0.8];
d.view3_axes_position = [0.25,0.25,0.5,0.5];

%Marker size for planets
d.msize = 2.5;

%Marker size of star trail
d.trail_size = 2;

%Store inter-mass distances (only for masses array, not rings or clusters)
d.store_inter_mass_distances = 1;

%

%% Additional inputs from a file %%
get_file = 1;
while get_file==1
    
    %Request user for inputs script
    [filename, pathname] = uigetfile('settings\*.m', 'Select gravity sim inputs');
    
    %Select option to load defaults if no file is chosen
    if filename==0
        % Construct a questdlg with three options
        choice = questdlg('Run the default options?', ...
            'Gravity Simulator', ...  %Title of question box
            'Default options','Choose another file','Quit',...   %Three options
            'Default options');  %Default when return key is pressed
        
        % Handle response
        switch choice
            case 'Default options'
                
                %Default inputs
                [d.masses,d.rings,d.clusters,d.G,d.Lmax,d.k,d.view_option,d.map_underlay,d.dt,d.Ng,d.Nr,d.RM] = gravity_sim_defaults;
                
                %Make sure choose file while loop does not continue
                get_file = 0;
                
            case 'Choose another file'
                %get_file = 1,  so the while loop will continue
                %through another iteration
            case 'Quit'
                %Quit the gravity simulator
                rmpath([pwd,'\functions']);
                return
        end
    else
        %Run inputs script
        [d.masses,d.rings,d.clusters,d.G,d.Lmax,d.k,d.view_option,d.map_underlay,d.dt,d.Ng,d.Nr,d.RM] = load_data([pathname,'\',filename]);
        get_file = 0;
    end
end

%

%% Generate vectors of parameters for masses, rings and clusters %%

%Start simulation time
d.t = 0;
d.t_store = 0;

%Masses
[x_m,y_m,z_m,vx_m,vy_m,vz_m,R_m,M_m,d.masses] = make_masses( d.masses );

%Compite initial inter-mass distance.
d.inter_mass_distance = distance( [x_m;y_m;z_m], [x_m;y_m;z_m] );

%Rings of massless planets
[x_r,y_r,z_r,vx_r,vy_r,vz_r,R_r,M_r,d.rings] = make_rings( d.rings, d.G );


%Add globular clusters of massless planets
[x_c,y_c,z_c,vx_c,vy_c,vz_c,R_c,M_c,d.clusters] = make_clusters( d.clusters, d.G );

%Assemble masses, rings and clusters into single vectors
d.x = [x_m,x_r,x_c]; clear x_m x_r x_c
d.y = [y_m,y_r,y_c]; clear y_m y_r y_c
d.z = [z_m,z_r,z_c]; clear z_m z_r z_c
d.vx = [vx_m,vx_r,vx_c]; clear vx_m vx_r vx_c
d.vy = [vy_m,vy_r,vy_c]; clear vy_m vy_r vy_c
d.vz = [vz_m,vz_r,vz_c]; clear vz_m vz_r vz_c
d.R = [R_m,R_r,R_c]; clear R_m R_r R_c
d.M = [M_m,M_r,M_c]; clear M_m M_r M_c

%Calculate maximum initial speed
d.max_speed = max(sqrt(d.vx.^2 + d.vy.^2 + d.vz.^2));

%

%% Set up figure and axes %%

%Create figure and axes
fig = figure('name','gravity_sim','WindowButtonDownFcn',@clickmouse,...
    'Toolbar','none','menubar','none','numbertitle','off',...
    'renderer','opengl','KeyPressFcn',@keypressfunc,...
    'BackingStore','off');
ax = axes('units','normalized');

%Get screensize and make figure fill the screen
screensize = get(0,'screensize');
set(fig,'units','pixels','position',screensize);

%

%% Initialise GUI user data control flags which can be modified via user %%
%% keyboard input %%
d.run = 1;
d.restart = 0;
d.load = 0;
d.save = 0;
d.print = 0;
d.wait = 0;
d.more_mass = 0;
d.less_mass = 0;
d.increase_AU = 0;
d.decrease_AU = 0;
d.up = 0;
d.down = 0;
d.left = 0;
d.right = 0;
d.m1_x = d.x(1);
d.m1_y = d.y(1);
d.button = '';
d.view = d.view_option;
d.view2 = d.view2_axes_position;
d.view3 = d.view3_axes_position;
d.map_underlay = d.map_underlay;
d.movie_frames = 0;
d.movie_dir = [];

%Set up axis
hold on;
xlabel('x')
ylabel('y')
zlabel('z')
axis vis3d
axis equal
xlim([-d.Lmax,d.Lmax])
ylim([-d.Lmax,d.Lmax])
zlim([-d.Lmax,d.Lmax])
view(d.view_option)
grid on;
hold on;

%Initialise movie frame number
movie_frame_num = 0;

%Define new colormap, making sure very low densities are a white colour.
cbar = colorbar;
set(cbar,'visible','off');
if (d.view==2) && (strcmp(d.map_underlay,'density')==1)
    
    %Mass number density
    colormap(density_colormap);
    caxis([0,d.DM*numel(d.M)/(4*d.Lmax^2)]);
    set(cbar,'visible','on');
elseif (d.view==2) && (strcmp(d.map_underlay,'speed')==1)
    
    %Speed map
    colormap(density_colormap);
    caxis([0,d,SM*d.max_speed*numel(d.M)/(4*d.Lmax^2)]);
    set(cbar,'visible','on');
end

%Plot initial positions of masses and density map, if selected
d = plot_initial_positions(d);

%Set axes view position
if d.view_option==3
    set(ax,'position',d.view3_axes_position);
else
    set(ax,'position',d.view2_axes_position);
end


%%

%% Run simulation %%

%Store original values of parameters to enable a reset when the 'r' key is
%pressed.
x_orig = d.x;
y_orig = d.y;
z_orig = d.z;
vx_orig = d.vx;
vy_orig = d.vy;
vz_orig = d.vz;
M_orig = d.M;

%Initialize 100 iteration time diagnostic
iter = 0;
timer_zero = tic;
d.hundred_iter_times = [];

%Run simulation until user presses 'q' (this changes d.run to 0)
loop = 1;
while d.run == 1
    %Obtain figure user data and check for key presses
    if d.wait == 1
        pause(0.1);
    elseif d.wait == 0
        
        %Timer diagnostic to see if gravity_sim.m is slowing down with time
        iter = iter+1;
        if iter == 100
            d.hundred_iter_times = [d.hundred_iter_times,toc(timer_zero)];
            timer_zero = tic;
            iter = 0;
        end
        
        %Update positions, veocities and accelerations and times
        [ d.x, d.y, d.z, d.vx, d.vy, d.vz, d.acc_x, d.acc_y, d.acc_z ] =...
            gravity( d.x, d.y, d.z, d.vx, d.vy, d.vz, d.R, d.M, d.G, d.Ng, d.Nr, d.RM, d.dt, 2 );
        d.t = d.t + d.dt;
        d.t_store = [d.t_store,d.t];
        d.x = d.x(2,:);
        d.y = d.y(2,:);
        d.z = d.z(2,:);
        d.vx = d.vx(2,:);
        d.vy = d.vy(2,:);
        d.vz = d.vz(2,:);
        d.acc_x = d.acc_x(2,:);
        d.acc_y = d.acc_y(2,:);
        d.acc_z = d.acc_z(2,:);
        
        %Modify velocities of objects that have gone beyond the screen
        %boundary, which represents an elastic sheet with coefficient of
        %restitution k
        xlims = get(gca,'xlim');
        ylims = get(gca,'ylim');
        zlims = get(gca,'zlim');
        
        %Right boundary
        i = find( (d.x>xlims(2)) & (d.vx>0) );
        if ~isempty(i) && ~isempty(d.k)
            d.vx(i) = -d.k*d.vx(i);
        end
        
        %Left boundary
        i = find( (d.x<xlims(1))  & (d.vx<0) );
        if ~isempty(i) && ~isempty(d.k)
            d.vx(i) = -d.k*d.vx(i);
        end
        
        %Top boundary
        i = find( (d.y>ylims(2)) & (d.vy>0) );
        if ~isempty(i) && ~isempty(d.k)
            d.vy(i) = -d.k*d.vy(i);
        end
        
        %Bottom boundary
        i = find( (d.y<ylims(1)) & (d.vy<0) );
        if ~isempty(i) && ~isempty(d.k)
            d.vy(i) = -d.k*d.vy(i);
        end
        
        %Upper z boundary
        i = find( (d.z>zlims(2)) & (d.vz>0) );
        if ~isempty(i) && ~isempty(d.k)
            d.vz(i) = -d.k*d.vz(i);
        end
        
        %Lower z boundary
        i = find( (d.z<zlims(1)) & (d.vz<0) );
        if ~isempty(i) && ~isempty(d.k)
            d.vz(i) = -d.k*d.vz(i);
        end
        
        %Update mass density or speed map
        if strcmp(d.map_underlay,'density')==1
            [xx,yy,new_rho] = mass_density_map( d.x,d.y,d.rho_r,d.N,d.Lmax);
            d.rho = 0.5 * ( d.rho + new_rho );
            set(d.rho_map,'CData',d.rho );
            
            %Mass number density
            colormap(density_colormap);
            caxis([0,d.DM*numel(d.M)/(4*d.Lmax^2)]);
            
        elseif strcmp(d.map_underlay,'speed')==1
            [xx,yy,new_speed] = speed_map( d.x,d.y,d.vx,d.vy,d.vz,d.rho_r,d.N,d.Lmax);
            d.speed = 0.5 * ( d.speed + new_speed );
            set(d.s_map,'CData',d.speed );
            
            %Speed map
            colormap(density_colormap);
            caxis([0,d.SM*d.max_speed*numel(d.M)/(4*d.Lmax^2)]);
        end
        
        %Update positions of masses
        
        %Masses
        if ~isempty(d.masses)
            for n=1:length(d.masses)
                if isempty( d.masses(n).marker_RGB )
                    [r,g,b] = x_to_color(sqrt(d.vx(n).^2 + d.vy(n).^2 + d.vz(n).^2 ),speed_colormap,0,d.max_speed);
                else
                    r = d.masses(n).marker_RGB(1);
                    g = d.masses(n).marker_RGB(2);
                    b = d.masses(n).marker_RGB(3);
                end
                
                %Update position
                set( d.p(n), 'Xdata', d.x(n), 'Ydata', d.y(n), 'Zdata', d.z(n),'markerfacecolor',[r,g,b],'markeredgecolor',[r,g,b]);
                
                %Add to trail
                if d.masses(n).plot_trail == 1
                    plot3( d.x(n),d.y(n),d.z(n),'r.','markerfacecolor',[r,g,b],'markeredgecolor',...
                        [r,g,b],'markersize',d.trail_size);
                end
            end
        end
        
        %Update inter-mass distance store
        if d.store_inter_mass_distances==1
            if ~isempty(d.masses)
                %Add a new page to this array
                nM = length(d.masses);
                d.inter_mass_distance(:,:,loop) = distance( [d.x(1:nM);d.y(1:nM);d.z(1:nM)],[d.x(1:nM);d.y(1:nM);d.z(1:nM)] );
            end
        end
        
        %Rings
        if ~isempty(d.rings)
            n0 = length(d.masses);
            end_index = n0;
            for n=1:length(d.rings)
                start_index = end_index+1;
                end_index = start_index + d.rings(n).num_masses - 1;
                if isempty( d.rings(n).marker_RGB )
                    [r,g,b] = x_to_color( mean( sqrt(d.vx(start_index:end_index).^2 +...
                        d.vy(start_index:end_index).^2 + d.vz(start_index:end_index).^2 )),speed_colormap,0,d.max_speed);
                else
                    r = d.rings(n).marker_RGB(1);
                    g = d.rings(n).marker_RGB(2);
                    b = d.rings(n).marker_RGB(3);
                end
                set( d.p(n0+n), 'Xdata', d.x(start_index:end_index), 'Ydata', d.y(start_index:end_index), 'Zdata',...
                    d.z(start_index:end_index),...
                    'markerfacecolor',[r,g,b],'markeredgecolor',[r,g,b]);
            end
        end
        
        %Clusters
        if ~isempty(d.clusters)
            n0 = n0 + length(d.rings);
            for n=1:length(d.clusters)
                start_index = end_index+1;
                end_index = start_index + d.clusters(n).num_masses - 1;
                if isempty( d.clusters(n).marker_RGB )
                    [r,g,b] = x_to_color( mean( sqrt(d.vx(start_index:end_index).^2 +...
                        d.vy(start_index:end_index).^2 + d.vz(start_index:end_index).^2 )),speed_colormap,0,d.max_speed);
                else
                    r = d.clusters(n).marker_RGB(1);
                    g = d.clusters(n).marker_RGB(2);
                    b = d.clusters(n).marker_RGB(3);
                end
                set( d.p(n0+n), 'Xdata', d.x(start_index:end_index), 'Ydata', d.y(start_index:end_index), 'Zdata',...
                    d.z(start_index:end_index),...
                    'markerfacecolor',[r,g,b],'markeredgecolor',[r,g,b]);
            end
        end
        
        %Flush pending graphics requests
        drawnow
        
        %Modify position of first mass, unless mouse is being used to
        %position it
        if strcmp(d.button, '')==1;
            d.m1_x = d.x(1);
            d.m1_y = d.y(1);
        end
        
        %Modify position if arrow keys ae pressed
        if d.up == 1;
            d.m1_y = d.m1_y + d.ds;
            d.up = 0;
        elseif d.down == 1;
            d.m1_y = d.m1_y - d.ds;
            d.down = 0;
        elseif d.right == 1;
            d.m1_x = d.m1_x + d.ds;
            d.right = 0;
        elseif d.left == 1;
            d.m1_x = d.m1_x - d.ds;
            d.left = 0;
        end
        
        %Modify sun position based upon location defined in d.m1_x and d.m1_y,
        d.x(1) = d.m1_x;
        d.y(1) = d.m1_y;
        
        %Modify mass if button is pressed
        if d.more_mass ==1
            d.M(1) = d.M(1) + d.dM;
            d.more_mass = 0;
        elseif d.less_mass ==1
            d.M(1) = d.M(1) - d.dM;
            d.less_mass = 0;
        end
        
        %Update title. Two versions, one for the screen, on for printed
        %images.
        d.title_str_print = {['t = ',num2str(d.t)]};
        d.title_str_sim = {['Mass moveable gravity simulation. M = ',num2str(d.M(1)),'. t = ',num2str(d.t),' / 2: 2D / 3: 3D / Arrows or mouse moves sun'],...
            'm: adds mass / n: removes mass / d: density map / v: speed map / x: no map / F1: Zoom out / F12: Zoom in',...
            '[: Start movie frames / ]: stop movie frames / p: screenshot / s: save data / l: load data / w: pause / c: continue / r: restart / q: quit'};
        if d.movie_frames == 1
            set(d.title,'string', d.title_str_print);
        else
            set(d.title,'string', d.title_str_sim);
        end
        
        %Change view and turn on or off speed or density map underlays
        if d.view == 3
            view(3);
            set(ax,'position',d.view3);
            set(d.rho_map,'visible','off');
            set(d.s_map,'visible','off');
            set(cbar,'visible','off');
            d.map_underlay = 'none';
        else
            view(2);
            set(ax,'position',d.view2);
            if strcmp(d.map_underlay,'density')==1
                set(d.rho_map,'visible','on');
                set(d.s_map,'visible','off');
                set(cbar,'visible','on');
            elseif strcmp(d.map_underlay,'speed')==1
                set(d.rho_map,'visible','off');
                set(d.s_map,'visible','on');
                set(cbar,'visible','on');
            else
                set(d.rho_map,'visible','off');
                set(d.s_map,'visible','off');
                set(cbar,'visible','off');
            end
        end
        
        %Save a .mat file of the current data set
        if d.save == 1
            d = orderfields(d);
            save( ['data\mgravity_data t=',num2str(d.t),...
                '. ',strrep(datestr(now),':','-'),'.mat'], 'd' );
            d.save = 0;
        end
        
        %Print a screenshot of the current view
        if d.print == 1;
            set(d.title,'string', d.title_str_print);
            print( fig, '-dpng','-r300',['mgravity. t=',num2str(d.t),...
                '. ',strrep(datestr(now),':','-'),'.png'])
            d.print = 0;
            set(d.title,'string', d.title_str_sim);
        end
        
        %Reset parameters if 'r' is pressed
        if d.restart == 1
            d.t = 0;
            d.x = x_orig;
            d.y = y_orig;
            d.z = z_orig;
            d.vx = vx_orig;
            d.vy = vy_orig;
            d.vz = vz_orig;
            d.M = M_orig;
            
            %Clear axis and start again
            cla;
            
            %Plot initial number density map and speed map
            d = plot_initial_positions(d);
            
            %Prevent further restarting until 'r' button is pressed again
            d.restart = 0;
        end
        
        %Load parameters if 'l' is pressed
        if d.load == 1
            [filename, pathname] = uigetfile('data\*.mat', 'Select new start data');
            if filename~=0
                %Load saved data structure d
                load([pathname,filename]);
                
                %Clear axis and start again
                cla;
                
                %Plot initial number density map and speed map
                d = plot_initial_positions(d);
                
                %Update 'original' parameters required for a restart when
                %'r' is pressed
                x_orig = d.x;
                y_orig = d.y;
                z_orig = d.z;
                vx_orig = d.vx;
                vy_orig = d.vy;
                vz_orig = d.vz;
                M_orig = d.M;
            end
            d.load = 0;
            
            %For some reason this process deactivates the figure! Hence the
            %following call is necessary
            figure(fig);
            
            %Reset loop number
            loop = 1;
        end
        
        %Modify axes limits by 10% if +/- buttons pressed
        if d.increase_AU ==1
            d.Lmax = d.Lmax * 1.1;
            d.increase_AU = 0;
        elseif d.decrease_AU ==1
            d.Lmax = d.Lmax / 1.1;
            d.decrease_AU = 0;
        end
        xlim([-d.Lmax,d.Lmax])
        ylim([-d.Lmax,d.Lmax])
        zlim([-d.Lmax,d.Lmax])
        
        %Write indexed frames to current movie directory. One assumes
        %a maximum of 99999 frames (!)
        if d.movie_frames == 1
            movie_frame_num = movie_frame_num + 1;
            print( fig, '-dpng','-r300',[d.movie_dir,'\',leadingzero(movie_frame_num,5),'.png']);
        else
            movie_frame_num = 0;
        end
        
        %Update loop number
        loop = loop + 1;
    end
end
rmpath([pwd,'\functions']);
close(fig);

%%

%% FUNCTIONS %%

%Figure key press function callback
function keypressfunc( fig,evnt )
global d
if evnt.Character == 'q'
    d.run = 0;
elseif evnt.Character == 'p'
    d.print = 1;
elseif evnt.Character == 's'
    d.save = 1;
elseif evnt.Character == 'w'
    d.wait = 1;
elseif evnt.Character == 'c'
    d.wait = 0;
elseif evnt.Character == 'm'
    d.more_mass = 1;
elseif evnt.Character == 'n'
    d.less_mass = 1;
elseif evnt.Character == 'r'
    if d.wait == 0
        d.restart = 1;
    end
elseif evnt.Character == 'l'
    if d.wait == 0
        d.load = 1;
    end
elseif strcmp(get(gcf,'currentkey'),'uparrow')==1
    d.up = 1;
elseif strcmp(get(gcf,'currentkey'),'downarrow')==1
    d.down = 1;
elseif strcmp(get(gcf,'currentkey'),'leftarrow')==1
    d.left = 1;
elseif strcmp(get(gcf,'currentkey'),'rightarrow')==1
    d.right = 1;
elseif strcmp(get(gcf,'currentkey'),'f1')==1
    d.increase_AU = 1;
elseif strcmp(get(gcf,'currentkey'),'f12')==1
    d.decrease_AU = 1;
elseif strcmp(get(gcf,'currentkey'),'2')==1
    d.view = 2;
elseif strcmp(get(gcf,'currentkey'),'3')==1
    d.view = 3;
elseif strcmp(get(gcf,'currentkey'),'v')==1
    d.map_underlay = 'speed';
elseif strcmp(get(gcf,'currentkey'),'d')==1
    d.map_underlay = 'density';
elseif strcmp(get(gcf,'currentkey'),'x')==1
    d.map_underlay = 'none';
elseif strcmp(get(gcf,'currentkey'),'leftbracket')==1
    d.movie_frames = 1;
    %Make directory to store movie frames
    d.movie_dir = ['mgravity movie ',strrep(datestr(now),':','-')];
    mkdir(d.movie_dir );
elseif strcmp(get(gcf,'currentkey'),'rightbracket')==1
    d.movie_frames = 0;
end

%%

%Function which runs when the mouse is clicked inside the figure window
function clickmouse(src,evnt)
global d

%Remove mouse pointer
set(gcf,'PointerShapeCData',nan(16,16),'Pointer','custom');

%Set callbacks for GUI
thisfig = gcbf();
set(thisfig,'WindowButtonMotionFcn',@dragmouse,'WindowButtonUpFcn',@unclickmouse);

%Start clock and display time elapsed in title bar
d.t0 = tic;
set(gcf,'name','gravity_sim')

%Set flag in gui which describes the type of mouse click
button = get(gcf, 'SelectionType');

%Double left mouse button click
if strcmp(button,'open')
    d.button = 'double click';
    
    %Single left mouse click
elseif strcmp(button,'normal')
    d.button = 'single click';
    
    %Shift + click or right & left mouse click
elseif strcmp(button,'extend')
    d.button = 'shift & click';
    
    %Modify figure name
    set(gcf,'name','gravity_sim - nothing more happens with a shift + click')
    
    %Control + click or right mouse click
elseif strcmp(button,'alt')
    d.button = 'ctrl & click';
    
    %Modify figure name
    set(gcf,'name','gravity_sim - nothing more happens with a ctrl + click')
else
    d.button='nothing';
end

%Store position when mouse was clicked
pos = get(gca,'CurrentPoint');
x = pos(1,1); y = pos(1,2);
d.start_x = x; d.start_y = y;

%%

%Function which runs when the mouse is dragged in the figure window
function dragmouse(src,evnt)
global d

%Get current position of the mouse
pos = get(gca,'CurrentPoint');
d.m1_x = pos(1,1); d.m1_y = pos(1,2);

%Update time in figure title bar
if strcmp(d.button,'ctrl & click') == 1
    set(gcf,'name',['gravity_sim dragging time = ',num2str( toc(d.t0) ),'s'] );
elseif strcmp(d.button,'shift & click') == 1
    set(gcf,'name',['gravity_sim dragging time = ',num2str( toc(d.t0) ),'s'] );
else
    set(gcf,'name',['gravity_sim dragging time = ',num2str( toc(d.t0) ),'s'] );
end

%Flush the graphics buffer and update GUI user data
drawnow;
set(gcf,'UserData',d);


%%

%Function which executes when the mouse button is released
function unclickmouse(src,evnt)
global d

%Get handle to figure currently invoking a callback
thisfig = gcbf();

%Reset mouse pointer to arrow. Reset mouse dragging functions to do nothing
set(thisfig,'Pointer','arrow','WindowButtonUpFcn','','WindowButtonMotionFcn','');

%Reset clock
d.t0 = 0;
d.button = '';

%Modify figure name
set(thisfig,'name','gravity_sim')

%%

%mass_rings
% Function which creates a vectors of x,y,z coordinates dscribing the
% initial position of rings of masses orbiting in a circular fashion about
% a centra mass mc. vx,vy and vz are the corresponding x,y,z velocities.
function [x0,y0,z0,vx0,vy0,vz0,num_masses] = mass_rings( xc,yc,zc,mc,arc_separation_AU,...
    num_rings, r0, ring_radius_diff_AU, d, alpha, G)

%Planets starting their orbits in concentric rings about star with coordinates
%(xc,yc,zc)
x0 = [];
y0 = [];
z0 = [];
vx0 = [];
vy0 = [];
vz0 = [];
for n=1:num_rings
    
    %Ring radius /AU
    r = r0 + ring_radius_diff_AU*(n-1);
    
    %Ring prbital speed /AU per Year
    v = sqrt( G *mc/r );
    
    %Number of planets per ring
    num_planets_per_ring = floor( 2*pi*r/arc_separation_AU );
    
    %Compute planet positions and velocities
    if num_planets_per_ring>=1
        theta = linspace(0,2*pi,num_planets_per_ring+1);
        for k=1:num_planets_per_ring
            x0 = [x0,xc+r*cos(theta(k))];
            vx0 = [vx0,-v*sin(theta(k))];
            y0 = [y0,yc+r*sin(theta(k))];
            vy0 = [vy0,v*cos(theta(k))];
            z0 = [z0,zc];
            vz0 = [vz0,0];
        end
    end
end

%Rotate to desired orientation
[x0,y0,z0] = point(x0,y0,z0,1,0,0,d(1),d(2),d(3));
[x0,y0,z0] = vrot(x0,y0,z0,d(1),d(2),d(3),d(1),d(2),d(3),alpha);
[vx0,vy0,vz0] = point(vx0,vy0,vz0,1,0,0,d(1),d(2),d(3));
[vx0,vy0,vz0] = vrot(vx0,vy0,vz0,d(1),d(2),d(3),d(1),d(2),d(3),alpha);

%Compute total number of masses
num_masses = length(x0);

%%

%gravity
% N-body gravity simulator. Uses the Verlet algorithm + Newtonian gravity +
% repulsion to determine future x,y,z coordinates (and velocities and accelerations)
% of N masses.
%
% [ t, x, y, z, vx, vy, vz, ax, ay, az ] =...
%    gravity( x0, y0, z0, vx0, vy0, vz0, R, M, g, dt, N )
%
% t               Vector of simulation times / Earth years
% x,y,z           x,y,z coordinates in AU (distance between the Earth and
%                 the Sun
% vx,vy,vz        Cartesian velocities in AU per Earth year
% ax,ay,az        Cartesian accelerations in AU per square Earth years.
%
% x0,y0,z0        Initial object positions in AU
% vx0,vy0,vz0     Initial ibject velocities in AU per Earth year
% R               Row vector of object raddi / Sun radii
% M               Row vector of object masses / Sun masses
% g               Strength of gravity / Gravitational force constant 6.67 * 10^-11
% dt              Time step, in  Earth years
% N               Number of time steps
function [ x, y, z, vx, vy, vz, ax, ay, az ] =...
    gravity( x0, y0, z0, vx0, vy0, vz0, R, M, G, Ng,Nr,RM, dt, N )

%Initialise output arrays
x = repmat( x0, [N,1] );
y = repmat( y0, [N,1] );
z = repmat( z0, [N,1] );
vx = repmat( vx0, [N,1] );
vy = repmat( vy0, [N,1] );
vz = repmat( vz0, [N,1] );

%Work out inital acceleration
[ax0,ay0,az0] = newton( x0,y0,z0,G,M,R,Ng,Nr,RM );
ax = repmat( ax0, [N,1] );
ay = repmat( ay0, [N,1] );
az = repmat( az0, [N,1] );

%Work out subsequent dynamics of all bodies via Verlet + Newtonian Gravity
for n=2:N
    
    %Update positions
    x(n,:) = x(n-1,:) + dt * vx(n-1,:) + 0.5*dt*dt*ax(n-1,:);
    y(n,:) = y(n-1,:) + dt * vy(n-1,:) + 0.5*dt*dt*ay(n-1,:);
    z(n,:) = z(n-1,:) + dt * vz(n-1,:) + 0.5*dt*dt*az(n-1,:);
    
    %Update acceleration
    [Ax,Ay,Az] = newton( x(n,:),y(n,:),z(n,:),G,M,R,Ng,Nr,RM );
    ax(n,:) = Ax(:);
    ay(n,:) = Ay(:);
    az(n,:) = Az(:);
    
    %Update velocity of objects using acceleration
    vx(n,:) = vx(n-1,:) + 0.5*dt*( ax(n,:) + ax(n-1,:) );
    vy(n,:) = vy(n-1,:) + 0.5*dt*( ay(n,:) + ay(n-1,:) );
    vz(n,:) = vz(n-1,:) + 0.5*dt*( az(n,:) + az(n-1,:) );
end

%%

%Compute acceleration of all objects via Newtonian-ish gravity (i.e.
%include a repulsion term to model close encounters)
% Note x,y,z are in units of AU, M are in Solar masses.
function [ax,ay,az] = newton( x,y,z,G,M,R,Ng,Nr,RM )

%Number of bodies
B = length(M);

%Work out bodies which don't have zero mass
nzm = find(M~=0);
zm = find(M==0);

%Initialise accelerations
ax = zeros(1,B);
ay = zeros(1,B);
az = zeros(1,B);

%

%Compute acceleration of non-zero masses

%Compute accelerations resulting from each mass and store these in a matrix aa.
d = distance( [x(nzm);y(nzm);z(nzm)],[x(nzm);y(nzm);z(nzm)] );
m = repmat( M(nzm).',[1,length(nzm)] );

%Compute repulsion constant
rc = RM * repmat( R(nzm).',[1,length(nzm)] );

%x component of acceleration
aa = m.*( repmat( x(nzm).',[1,length(nzm)] ) - repmat( x(nzm),[length(nzm),1] ))./d.^(Ng+1) -...
    rc.*m.*( repmat( x(nzm).',[1,length(nzm)] ) - repmat( x(nzm),[length(nzm),1] ))./d.^(Nr+1);
aa(1:length(nzm)+1:end) = 0;  %Don't let self same mass displacements contribute!
ax(nzm) = G*sum(aa,1);

%y component of acceleration
aa = m.*( repmat( y(nzm).',[1,length(nzm)] ) - repmat( y(nzm),[length(nzm),1] ))./d.^(Ng+1) -...
    rc.*m.*( repmat( y(nzm).',[1,length(nzm)] ) - repmat( y(nzm),[length(nzm),1] ))./d.^(Nr+1);
aa(1:length(nzm)+1:end) = 0; %Don't let self same mass displacements contribute!
ay(nzm) = G*sum(aa,1);

%z component of acceleration
aa = m.*( repmat( z(nzm).',[1,length(nzm)] ) - repmat( z(nzm),[length(nzm),1] ))./d.^(Ng+1) -...
    rc.*m.*( repmat( z(nzm).',[1,length(nzm)] ) - repmat( z(nzm),[length(nzm),1] ))./d.^(Nr+1);
aa(1:length(nzm)+1:end) = 0; %Don't let self same mass displacements contribute!
az(nzm) = G*sum(aa,1);

%

%Compute acceleration of zero masses

%Compute accelerations resulting from each mass and store these in a matrix aa.
d = distance( [x(nzm);y(nzm);z(nzm)],[x(zm);y(zm);z(zm)] );
m = repmat( M(nzm).',[1,length(zm)] );

%Compute repulsion constant
rc = repmat( R(nzm).',[1,length(zm)] );

%x component of acceleration
aa = m.*( repmat( x(nzm).',[1,length(zm)] ) - repmat( x(zm),[length(nzm),1] ))./d.^(Ng+1) -...
    rc.*m.*( repmat( x(nzm).',[1,length(zm)] ) - repmat( x(zm),[length(nzm),1] ))./d.^(Nr+1);
ax(zm) = G*sum(aa,1);

%y component of acceleration
aa = m.*( repmat( y(nzm).',[1,length(zm)] ) - repmat( y(zm),[length(nzm),1] ))./d.^(Ng+1) -...
    rc.*m.*( repmat( y(nzm).',[1,length(zm)] ) - repmat( y(zm),[length(nzm),1] ))./d.^(Nr+1);
ay(zm) = G*sum(aa,1);

%z component of acceleration
aa = m.*( repmat( z(nzm).',[1,length(zm)] ) - repmat( z(zm),[length(nzm),1] ))./d.^(Ng+1) -...
    rc.*m.*( repmat( z(nzm).',[1,length(zm)] ) - repmat( z(zm),[length(nzm),1] ))./d.^(Nr+1);
az(zm) = G*sum(aa,1);

%%

%Create mass density map. This creates an array of mass densities within a
%specified radius r about an N*N equispaced grid of points within the gravity
%simulation grid. This can be used to plot a smoothed surface which varies
%dynamically with the distribution of masses. By mass density we actually
%mean number density i.e. the number of objects within a particular radius.
function [xx,yy,rho] = mass_density_map( x,y,r,N,Lmax)

%Create grid
xx = linspace(-Lmax,Lmax,N);
yy = linspace(-Lmax,Lmax,N);
[xx,yy] = meshgrid(xx,yy);

%Turn grid points into row vectors to enable vectorized distance
%calculation
xx = reshape(xx,[1,N*N]);
yy = reshape(yy,[1,N*N]);

%Step through grid and find the number of masses within radius r of the
%grid point
d = distance( [xx;yy],[x;y] );

%Compute density map
rho = sum( d<r, 2 ).';

%Scale by size of area
rho = rho / (pi*r^2);

%Reshape to square arrays
xx = reshape(xx,[N,N]);
yy = reshape(yy,[N,N]);
rho = reshape(rho,[N,N]);

%%

%Create speed map. This creates an array of object speeds within a
%specified radius r about an N*N equispaced grid of points within the gravity
%simulation grid. This can be used to plot a smoothed surface which varies
%dynamically with the distribution of masses. By mass density we actually
%mean number density i.e. the number of objects within a particular radius.
function [xx,yy,speed] = speed_map( x,y,vx,vy,vz,r,N,Lmax)

%Create grid
xx = linspace(-Lmax,Lmax,N);
yy = linspace(-Lmax,Lmax,N);
[xx,yy] = meshgrid(xx,yy);

%Turn grid points into row vectors to enable vectorized distance
%calculation
xx = reshape(xx,[1,N*N]);
yy = reshape(yy,[1,N*N]);

%Form a vector of object speeds and then replicate by the number of grid
%points
s = sqrt( vx.^2 + vy.^2 + vz.^2 );
s = repmat( s, [N*N,1] );

%Step through grid and find the number of masses within radius r of the
%grid point
d = distance( [xx;yy],[x;y] );

%Compute average speed map
% speed = sum( s.*(d<r) ,2 ).' ./ sum( d<r, 2 ).';
speed = sum( s.*(d<r) ,2 ).';

%Scale by size of area
speed = speed / (pi*r^2);

%Set /0 rows (i.e. where an 'object bin' contains no masses) to be 0
speed(isnan(speed)) = 0;
speed(speed==Inf) = 0;

%Reshape to square arrays
xx = reshape(xx,[N,N]);
yy = reshape(yy,[N,N]);
speed = reshape(speed,[N,N]);

%%

%distance
% Function which computes Euclidean distance matrix.
% This fully vectorized (VERY FAST!) m-file computes the
% Euclidean distance between two vectors by:
%          ||A-B|| = sqrt ( ||A||^2 + ||B||^2 - 2*A.B )
%
% Syntax: E = distance(A,B)
%    A - (DxM) matrix
%    B - (DxN) matrix
%    E - (MxN) Euclidean distances between vectors in A and B
%
% Example :
%    A = rand(400,100); B = rand(400,200);
%    d = distance(A,B);

% Author   : Roland Bunschoten
%            University of Amsterdam
%            Intelligent Autonomous Systems (IAS) group
%            Kruislaan 403  1098 SJ Amsterdam
%            tel.(+31)20-5257524
%            bunschot@wins.uva.nl
% Last Rev : Oct 29 16:35:48 MET DST 1999
% Tested   : PC Matlab v5.2 and Solaris Matlab v5.3
% Thanx    : Nikos Vlassis
function d = distance(a,b)
if (nargin ~= 2)
    error('Not enough input arguments');
end
if (size(a,1) ~= size(b,1))
    error('A and B should be of same dimensionality');
end
aa=sum(a.*a,1); bb=sum(b.*b,1); ab=a'*b;
d = sqrt(abs(repmat(aa',[1 size(bb,2)]) + repmat(bb,[size(aa,2) 1]) - 2*ab));

%%

%interp_colormap
% Function which interpolates current colourmap to yield better graduated
% shading. N is number of possible colours.
function interp_colormap( N )

%Get current colourmap
map = colormap;

%Initialise new colormap
new_map = ones(N,3);

%Get size of current colormap and initalise red,green,blue vectors
dim = size(map);
R = ones(1,dim(1));
G = ones(1,dim(1));
B = ones(1,dim(1));
RR = ones(1,N);
GG = ones(1,N);
BB = ones(1,N);

%Populate these with current colormap
R(:) = map(:,1);
G(:) = map(:,2);
B(:) = map(:,3);

%Interpolate to yield new colour map
x = linspace( 1, dim(1), N );
RR = interp1( 1:dim(1), R, x );
GG = interp1( 1:dim(1), G, x );
BB = interp1( 1:dim(1), B, x );
new_map(:,1) = RR(:);
new_map(:,2) = GG(:);
new_map(:,3) = BB(:);

%Set colormap to be new map
colormap( new_map );

%%

%Make RGB color values from a scalar and a colormap
function [r,g,b] = x_to_color(x,cmap,xmin,xmax)
red = cmap(:,1).';
green = cmap(:,2).';
blue = cmap(:,3).';
cx = linspace( 0, 1, length(red) );
x = (x - xmin)/(xmax - xmin);
x(x>1) = 0.99;
x(x<0) = 0.01;
r = interp1( cx, red, x );
g = interp1( cx, green, x );
b = interp1( cx, blue, x );
r(r>1)=1; r(r<0)=0;
g(g>1)=1; g(g<0)=0;
b(b>1)=1; b(b<0)=0;

%%

%Convert a number x into a string of N characters
function s = leadingzero(x,N)
s = [strrep(blanks(N-1 - fix(log10(x))),' ','0'),...
    num2str(x,['%',num2str(N),'.0f'])];

%%

%Colormap used to indicate speed of masses
function map = speed_colormap
N = 200;
x = linspace(0,1,N);

%Start color
R1 = 0;
G1 = 0;
B1 = 1;

%End color
R2 = 1;
G2 = 0.4;
B2 = 0;

%Create interpolated colormap
R = interp1([0,1],[R1,R2],x).';
G = interp1([0,1],[G1,G2],x).';
B = interp1([0,1],[B1,B2],x).';
map = [R,G,B];

%%

%Colormap used to indicate number density of masses
function map = density_colormap
N = 200;
map = colormap('jet');
map = [[linspace(1,map(30,1),30).',...
    linspace(1,map(30,2),30).',...
    linspace(1,map(30,3),30).'];...
    map(30:end-10,:)];
colormap(map);
interp_colormap( N );
map = colormap;

%%

%Another colormap used to indicate number density of masses
function map = density_colormap2
N = 200;
x = linspace(0,1,N);

%Start color
R1 = 1;
G1 = 1;
B1 = 1;

%Middle color
R2 = 0;
G2 = 0;
B2 = 1;

%End color
R3 = 1;
G3 = 0;
B3 = 0;

%Create interpolated colormap
R = interp1([0,0.7,1],[R1,R2,R3],x).';
G = interp1([0,0.7,1],[G1,G2,G3],x).';
B = interp1([0,0.7,1],[B1,B2,B3],x).';
map = [R,G,B];

%%

%globular_cluster
% Function which creates arrays of positions and velocities of masses that
% will form circular orbits aboyt a central mass, filling a spherical volume
% around the central masses.
function [x0,y0,z0,vx0,vy0,vz0,num_masses] = globular_cluster( xc,yc,zc,mc,shell_separation_AU,...
    num_shells, r0, num_masses_per_square_AU, G)

%Planets starting their orbits in spherical shells about star with coordinates
%(xc,yc,zc)
x0 = [];
y0 = [];
z0 = [];
vx0 = [];
vy0 = [];
vz0 = [];
for n=1:num_shells
    
    %Ring radius /AU
    r = r0 + shell_separation_AU*(n-1);
    
    %Ring prbital speed /AU per Year
    v = sqrt( G *mc/r );
    
    %Number of planets per shell
    s = floor( sqrt( num_masses_per_square_AU*4*pi*r^2 ));
    azi = linspace(0,2*pi,s+1);
    elev = linspace(-pi,pi,s+1);
    
    %Compute planet positions and velocities
    if s >=1
        for i=1:s
            for j=1:s
                %x
                x0 = [x0,xc+r*cos(elev(i))*sin(azi(j))];
                vx0 = [vx0,-v*sin(elev(i))*sin(azi(j))];
                
                %y
                y0 = [y0,yc+r*cos(elev(i))*cos(azi(j))];
                vy0 = [vy0,-v*sin(elev(i))*cos(azi(j))];
                
                %z
                z0 = [z0,zc + r*sin(elev(i))];
                vz0 = [vz0,v*cos(elev(i))];
            end
        end
    end
end

%Compute total number of masses
num_masses = length(x0);

%%

%Make vector of masses from input parameters
function [x,y,z,vx,vy,vz,R,M,masses] = make_masses( masses )
x = [];
y = [];
z = [];
vx = [];
vy = [];
vz = [];
R = [];
M = [];
if ~isempty(masses)
    for n=1:length(masses)
        x = [x,masses(n).x0];
        y = [y,masses(n).y0];
        z = [z,masses(n).z0];
        vx = [vx,masses(n).vx0];
        vy = [vy,masses(n).vy0];
        vz = [vz,masses(n).vz0];
        R = [R,masses(n).radii];
        M = [M,masses(n).mass];
    end
end

%%

%Make vector of massless rings from input parameters
function [x,y,z,vx,vy,vz,R,M,rings] = make_rings( rings, G )
x = [];
y = [];
z = [];
vx = [];
vy = [];
vz = [];
R = [];
M = [];
if ~isempty(rings)
    for n=1:length(rings)
        [xx,yy,zz,vxx,vyy,vzz,num_masses] = mass_rings( rings(n).xc,rings(n).yc,rings(n).zc,....
            rings(n).mass_at_centre, rings(n).arc_separation_AU,...
            rings(n).num_rings, rings(n).first_ring_radius_AU,...
            rings(n).ring_radius_diff_AU, rings(n).d, rings(n).alpha, G);
        x = [x,xx];
        y = [y,yy];
        z = [z,zz];
        vx = [vx,vxx + rings(n).vxc];
        vy = [vy,vyy + rings(n).vyc];
        vz = [vz,vzz + rings(n).vzc];
        R = [R,zeros(1,num_masses)];
        M = [M,zeros(1,num_masses)];
        rings(n).num_masses = num_masses;
    end
end

%%

%Make globular clusters from input parameters
function [x,y,z,vx,vy,vz,R,M,clusters] = make_clusters( clusters, G )
x = [];
y = [];
z = [];
vx = [];
vy = [];
vz = [];
R = [];
M = [];
if ~isempty(clusters)
    for n=1:length(clusters)
        [xx,yy,zz,vxx,vyy,vzz,num_masses] = globular_cluster( clusters(n).xc,clusters(n).yc,clusters(n).zc,...
            clusters(n).mass_at_centre,clusters(n).shell_separation_AU,...
            clusters(n).num_shells, clusters(n).first_shell_radius, clusters(n).num_masses_per_square_AU, G);
        x = [x,xx];
        y = [y,yy];
        z = [z,zz];
        vx = [vx,vxx + clusters(n).vxc];
        vy = [vy,vyy + clusters(n).vyc];
        vz = [vz,vzz + clusters(n).vzc];
        R = [R,zeros(1,num_masses)];
        M = [M,zeros(1,num_masses)];
        clusters(n).num_masses = num_masses;
    end
end

%%

%Plot initial positions of masses and density map, if selected
function d = plot_initial_positions(d)

%Plot initial number density map and speed map
[xx,yy,d.rho] = mass_density_map( d.x,d.y,d.rho_r,d.N,d.Lmax);
d.rho_map = surf(xx,yy,-ones(size(xx)),d.rho,'visible','off');
[xx,yy,d.speed] = speed_map( d.x,d.y,d.vx,d.vy,d.vz,d.rho_r,d.N,d.Lmax);
d.s_map = surf(xx,yy,-ones(size(xx)),d.speed,'visible','off');
shading interp;
set(d.rho_map,'facealpha',0.7)
set(d.s_map,'facealpha',0.7)
if ( d.view == 2 ) && ( strcmp(d.map_underlay,'density')==1 )
    set(d.rho_map,'visible','on');
elseif( d.view == 2 ) && ( strcmp(d.map_underlay,'speed')==1 )
    set(d.s_map,'visible','on');
end

%Plot initial positions of masses

%Masses
if ~isempty(d.masses)
    for n=1:length(d.masses)
        if isempty( d.masses(n).marker_RGB )
            [r,g,b] = x_to_color(sqrt(d.vx(n).^2 + d.vy(n).^2 + d.vz(n).^2 ),speed_colormap,0,d.max_speed);
        else
            r = d.masses(n).marker_RGB(1);
            g = d.masses(n).marker_RGB(2);
            b = d.masses(n).marker_RGB(3);
        end
        d.p(n) = plot3(d.x(n),d.y(n),d.z(n),'ro','markerfacecolor',[r,g,b],'markeredgecolor',...
            [r,g,b],'markersize',d.masses(n).marker_size);
    end
end

%Rings
if ~isempty(d.rings)
    n0 = length(d.masses);
    end_index = n0;
    for n=1:length(d.rings)
        start_index = end_index+1;
        end_index = start_index + d.rings(n).num_masses - 1;
        if isempty( d.rings(n).marker_RGB )
            [r,g,b] = x_to_color( mean( sqrt(d.vx(start_index:end_index).^2 +...
                d.vy(start_index:end_index).^2 + d.vz(start_index:end_index).^2 )),speed_colormap,0,d.max_speed);
        else
            r = d.rings(n).marker_RGB(1);
            g = d.rings(n).marker_RGB(2);
            b = d.rings(n).marker_RGB(3);
        end
        d.p(n0+n) = plot3(d.x(start_index:end_index),d.y(start_index:end_index),d.z(start_index:end_index),'r.','markersize',d.msize,...
            'markerfacecolor',[r,g,b],'markeredgecolor',[r,g,b]);
    end
end

%Clusters
if ~isempty(d.clusters)
    n0 = n0 + length(d.rings);
    for n=1:length(d.clusters)
        start_index = end_index+1;
        end_index = start_index + d.clusters(n).num_masses - 1;
        if isempty( d.clusters(n).marker_RGB )
            [r,g,b] = x_to_color( mean( sqrt(d.vx(start_index:end_index).^2 +...
                d.vy(start_index:end_index).^2 + d.vz(start_index:end_index).^2 )),speed_colormap,0,d.max_speed);
        else
            r = d.clusters(n).marker_RGB(1);
            g = d.clusters(n).marker_RGB(2);
            b = d.clusters(n).marker_RGB(3);
        end
        d.p(n0+n) = plot3(d.x(start_index:end_index),d.y(start_index:end_index),d.z(start_index:end_index),'r.','markersize',d.msize,...
            'markerfacecolor',[r,g,b],'markeredgecolor',[r,g,b]);
    end
end

%Set title
d.title_str  = {['Mass moveable gravity simulation. M = ',num2str(d.M(1)),' / 2: 2D / 3: 3D / Arrows or mouse moves sun'],...
    'm: adds mass / n: removes mass / d: density map / v: speed map / x: no map / F1: Zoom out / F12: Zoom in',...
    '[: Start movie frames / ]: stop movie frames / p: screenshot / s: save data / l: load data / w: pause / c: continue / r: restart / q: quit'};
d.title = title(d.title_str);

%%

%Default settings
function [masses,rings,clusters,G,Lmax,k,view_option,map_underlay,dt,Ng,Nr,RM] = gravity_sim_defaults

%Astronomical parameters in SI units
AU = 149.6e9;
G = 6.67e-11;
M_sun = 2e30;
Yr = 365*24*3600;

%Dimensionless number which controls the dynamics, and results from the
%scaling of mass, distance and time parameters to make them dimensionless.
G = G*Yr^2*M_sun/(AU^3);

%Strength of gravity
G_factor = 1;
G = G_factor*G;

%Timestep
dt = 0.05;

%Gravitation exponent
Ng = 2;

%Repulsion exponent
Nr = 5;

%Repulsion magnitude
RM = 0.01;

%Set dimensions of space
Lmax = 8;

%Star masses in Solar masses
M1 = 1;
M2 = 1;

%Initial star separation
s = 6;

%Initial phase (polar angle) of two-mass elliptical orbit
theta0=0;

%Elliptical orbit eccentricity
ecc = 0;

%Orientation of semi-major axis of ellipse
d = [1,0,0];

%Rotation of ellipse clockwise about d /radians
alpha = 0;

%Centre of mass initial coordinates
h0 = [0,0,0];

%Cetre of mass velocity
hdot = [0,0,0];

%Boundary coefficient of restitution (if [] then there is
%no boundary)
k = [];

%Initial view option
view_option = 3;

%For 2D views, choose density or speed map underlay ('density', 'speed',
%'none')
map_underlay = 'none';

%Binary star circular orbit conditions
[ux,vx1,vx2,...
    uy,vy1,vy2,...
    uz,vz1,vz2,...
    xc,x1,x2,...
    yc,y1,y2,...
    zc,z1,z2,...
    theta,theta_dot,t,j,J,E,P] =...
    two_body_elliptical_orbit(  G, M1, M2, s, ecc, theta0, h0, hdot, d, alpha, 1 );

%Star (the position of this one can be modified)
masses(1).name = 'Sun #1';
masses(1).mass = M1;
masses(1).radii = 0.01;
masses(1).x0 = x1;
masses(1).y0 = y1;
masses(1).z0 = z1;
masses(1).vx0 = vx1;
masses(1).vy0 = vy1;
masses(1).vz0 = vz1;
masses(1).marker_RGB = [];
masses(1).marker_size = 5;
masses(1).plot_trail = 1;
masses(1).init_orbit_M = M2;
masses(1).init_orbit_G = G;
masses(1).init_orbit_P = P;
masses(1).init_orbit_ecc = ecc;
masses(1).init_orbit_s = s;
masses(1).init_orbit_h0 = h0;
masses(1).init_orbit_hdot = hdot;
masses(1).init_orbit_d = d;
masses(1).init_orbit_alpha = alpha;
masses(1).init_orbit_theta0 = theta0;
masses(1).init_orbit_thetadot0 = theta_dot;
masses(1).init_orbit_J = J;
masses(1).init_orbit_j = j;
masses(1).init_orbit_E = E;

%Star (the position of this one can be modified)
masses(2).name = 'Sun #2';
masses(2).mass = M2;
masses(2).radii = 0.01;
masses(2).x0 = x2;
masses(2).y0 = y2;
masses(2).z0 = z2;
masses(2).vx0 = vx2;
masses(2).vy0 = vy2;
masses(2).vz0 = vz2;
masses(2).marker_RGB = [];
masses(2).marker_size = 5;
masses(2).plot_trail = 1;
masses(2).init_orbit_M = M1;
masses(2).init_orbit_G = G;
masses(2).init_orbit_P = P;
masses(2).init_orbit_ecc = ecc;
masses(2).init_orbit_s = s;
masses(2).init_orbit_h0 = h0;
masses(2).init_orbit_hdot = hdot;
masses(2).init_orbit_d = d;
masses(2).init_orbit_alpha = alpha;
masses(2).init_orbit_theta0 = theta0;
masses(2).init_orbit_thetadot0 = theta_dot;
masses(2).init_orbit_J = J;
masses(2).init_orbit_j = j;
masses(2).init_orbit_E = E;

%Rings (1)
rings(1).xc = x1;
rings(1).yc = y1;
rings(1).zc = z1;
rings(1).vxc = vx1;
rings(1).vyc = vy1;
rings(1).vzc = vz1;
rings(1).num_rings = 30;
rings(1).arc_separation_AU = 1*pi/30;
rings(1).first_ring_radius_AU = 2;
rings(1).ring_radius_diff_AU = 0.1;
rings(1).d = [1,0,0];
rings(1).alpha = 0;
rings(1).mass_at_centre = M1;
rings(1).marker_RGB = [0,0,0];
rings(1).marker_size = 1;

%Rings (2)
rings(2).xc = x2;
rings(2).yc = y2;
rings(2).zc = z2;
rings(2).vxc = vx2;
rings(2).vyc = vy2;
rings(2).vzc = vz2;
rings(2).num_rings = 30;
rings(2).arc_separation_AU = 1*pi/30;
rings(2).first_ring_radius_AU = 2;
rings(2).ring_radius_diff_AU = 0.1;
rings(2).d = [1,0,0];
rings(2).alpha = 0;
rings(2).mass_at_centre = M2;
rings(2).marker_RGB = [1,0,1];
rings(2).marker_size = 1;

%Clusters
clusters = [];

%%

%Get settings from file
function [masses,rings,clusters,G,Lmax,k,view_option,map_underlay,dt,Ng,Nr,RM]  = load_data( data_name )

%Note the desired outputs will be defined in the data file, which is a
%MATLAB script. All other parameters will be ignored since the file is run
%within this function.
run(data_name);

%End of code