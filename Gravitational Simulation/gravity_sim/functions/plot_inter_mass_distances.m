%plot_inter_mass_distances
% Function which loads data saved by gravity_sim and plots the inter-mass
% distances vs time

function plot_inter_mass_distances

%% INPUTS %%

%Choose specific pairs to plot. The text must match the data!
%{} will give all the data.
chosen_pairs = {
    'Sun & Pluto',...
    'Sun & Neptune',...
    'Sun & Uranus',...
    'Sun & Saturn',...
    'Sun & Jupiter'
    };
%chosen_pairs = {};

%%

%Get data
[filename,pathname] = uigetfile('..\data\*.mat','Select gravity_sim data');
if ~isempty(filename)
    
    %Load data
    load( [pathname,filename] );
    
    %Determine dimensions
    dim = size(d.inter_mass_distance);
    
    %Initialize inter-mass distance matrix
    IMDM = zeros(dim(3), 0.5 * ( dim(1)^2 - dim(1) ) );
    
    %Extract inter-mass-distance matrix in a readily plotable form and determine plot legend
    dim_IMDM = size(IMDM);
    legend_str = cell(1,dim_IMDM(2));
    ignore = triu( ones(length(d.masses)), 1);
    k=0;
    for r=1:length(d.masses)
        for c=1:length(d.masses)
            if ignore(r,c)~=0
                k=k+1;
                legend_str{k} = [d.masses(r).name,' & ',d.masses(c).name];
                IMDM(:,k) = d.inter_mass_distance(r,c,:);
            end
        end
    end
    
    %Work out chosen pairings
    if isempty(chosen_pairs)
        pairings = 1:dim_IMDM(2);
    else
        pairings = [];
        new_legend_str = {};
        for n=1:length(chosen_pairs);
            i = strcmp( chosen_pairs{n}, legend_str );
                pairings = [pairings,find(i==1)];
                new_legend_str = [new_legend_str,legend_str{i}];
        end
        
        %Check if there are no pairings. In this case plot everything
        if isempty(pairings)
            pairings = 1:dim_IMDM(2);
        else
            legend_str = new_legend_str;
        end       
    end
    
    %Set up figure and axes
    figure('color',[1 1 1],'name','Inter-mass distances')
    axes('Nextplot','add')
       
    %Interpolate curent colormap to determine line colours
    cmap = colormap('jet');
    dim_cmap = size(cmap);
    colour_i = round(linspace(1,dim_cmap(1),length(pairings)));
    
    %Plot inter-mass distances
    for n=1:length(pairings)
        plot( d.t_store(1:dim_IMDM(1)), IMDM( :,pairings(n) ).','color',cmap( colour_i(n), : ) );
    end
    
    legend(legend_str,'location','eastoutside');
    xlabel('time')
    ylabel('Inter-mass distance')
    title('Gravity sim inter-mass distances')
    
    %Save PNG file of output
    print( gcf, '-dpng','-r300',['Inter-mass distances',...
        '. ',strrep(datestr(now),':','-'),'.png'])
    
end

%End of code