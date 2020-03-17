%% rm_raincloud_cg - plots a raincloud of different groups and observations
% Use like: h = rm_raincloud_cg(data)
% Where 'data' is an M x N cell array of M measurements and N data series
% See below for optional inputs.
%
% Modified: Cristina Gil, 17.03.2020

function h = rm_raincloud_cg(data, varargin)
%% ---------------------------- INPUT ----------------------------
%
% data - M x N cell array of M measurements and N data series
%
% --------------------- OPTIONAL ARGUMENTS ----------------------
%
% colours               - N x 3 array defining the colour to plot each data series % color vector for rainclouds (default gray, i.e. = [.5 .5 .5])
% density_type          - choice of density algo ('ks' or 'rath'). Default = 'ks'
% bandwidth             - bandwidth of smoothing kernel (default = 1)
% plot_top_to_bottom    - logical to plot top-to-bottom (default=1)
% plot_median_lines     - logical set to 0 if you don't want median lines plotted (default = 0)
% plot_median_dots      - logical set to 0 if you don't want median dots plotted (default = 0)
% box_on                - logical to turn box plots on/off (default = 0)
% line_width            - scalar value to set global line width (default = 2)
% bxcl                  - color of box outline
% box_col_match         - logical to set it so that boxes match the colour of clouds (default = 0)
% box_dodge             - logical to turn on/off box plot dodging (default = 0)
% raindrop_size         - scalar positive value to control the size of the raindrops (default = 3)
% alpha                 - scalar positive value to increase cloud alpha (defalut = 1)
% dist_plots            - scalar that defines the distance between plots (default = 1.5)  

% ---------------------------- OUTPUT ----------------------------
% h is a cell array containing handles of the various figure parts:
% h.p{i,j}  is the handle to the density plot from data{i,j}
% h.s{i,j}  is the handle to the 'raindrops' (individual datapoints) from data{i,j}
% h.b1{i,j} [optional: only if box_on is true] is the handle for the mean line of the boxplots
% h.b2{i,j} [optional: only if box_on is true] is the handle for right whisker of the boxplots
% h.b3{i,j} [optional: only if box_on is true] is the handle for left whisker of the boxplots
% h.m(i,j)  [optional: only if plot_mean_dots is true] is the handle to the single, large dot that represents mean(data{i,j})
% h.l(i,j)  [optional: only if plot_mean_lines is true] is the handle for the line connecting h.m(i,j) and h.m(i+1,j)
%
% ------------------------ EXAMPLE USAGE -------------------------
% h = rm_raincloud_cg(raincloudData,'colours',[0 .3961 .7412; .8902 .4471 .1333],'plot_mean_lines',0,'box_on',1)
% 
%% check dimensions of data
% input parsing settings
p = inputParser;
p.CaseSensitive = true;
p.Parameters;
p.Results;
p.KeepUnmatched = true;
validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x > 0);

[n_plots_per_series, n_series] = size(data);

%% default arguments
% set the desired and optional input arguments
addRequired(p, 'data', @iscell);
addOptional(p, 'colours', [0.5 0.5 0.5; 1 1 1; 0 0 0], @isnumeric)
addOptional(p, 'density_type', 'ks', @ischar)
addOptional(p, 'bandwidth', [])
addOptional(p, 'plot_top_to_bottom', 1, @isnumeric)
addOptional(p, 'plot_median_lines', 1, @isnumeric)
addOptional(p, 'plot_median_dots', 0, @isnumeric)
addOptional(p, 'box_on', 0, @isnumeric)
addOptional(p, 'bxcl', [0 0 0], @isnumeric)
addOptional(p, 'line_width', 2, validScalarPosNum)
addOptional(p, 'box_col_match', 0, @isnumeric)
addOptional(p, 'box_dodge', 0, @isnumeric)
addOptional(p, 'raindrop_size', 3, validScalarPosNum)
addOptional(p, 'alpha', 0.5, @isnumeric)
addOptional(p, 'dist_plots', 1.5, @isnumeric)



% parse the input
parse(p,data,varargin{:});
% then set/get all the inputs out of this structure
data                = p.Results.data;
colours             = p.Results.colours;
bandwidth           = p.Results.bandwidth;
density_type        = p.Results.density_type;
plot_top_to_bottom  = p.Results.plot_top_to_bottom;
plot_median_lines   = p.Results.plot_median_lines;
plot_median_dots    = p.Results.plot_median_dots;
box_on              = p.Results.box_on;
bxcl                = p.Results.bxcl;
line_width          = p.Results.line_width;
box_col_match       = p.Results.box_col_match;
box_dodge           = p.Results.box_dodge;
raindrop_size       = p.Results.raindrop_size;
alpha               = p.Results.alpha;
dist_plots          = p.Results.dist_plots;



%% Calculate properties of density plots

% Probably okay to hard-code this as it just determines the granularity of
% the density estimate
density_granularity = 200;

n_bins = repmat(density_granularity, n_plots_per_series, n_series);

% initialize variables
ks = cell(size(data));
x = cell(size(data));
q = cell(size(data));
faces = cell(size(data));

% calculate kernel densities
for i = 1:n_plots_per_series
    for j = 1:n_series
       
        switch density_type
            
            case 'ks'
                
                % compute density using 'ksdensity'
                [ks{i, j}, x{i, j}] = ksdensity(data{i, j}, 'NumPoints', n_bins(i, j), 'bandwidth', bandwidth);
                
            case 'rash'
                
                % check for rst_RASH function (from Robust stats toolbox) in path, fail if not found 
                assert(exist('rst_RASH', 'file') == 2, 'Could not compute density using RASH method. Do you have the Robust Stats toolbox on your path?');
                
                % compute density using RASH
                [x{i, j}, ks{i, j}] = rst_RASH(data{i, j});
                
                % override default 'n_bins' as rst_RASH determines number of bins
                n_bins(i, j) = size(ks{i, j}, 2);
        end
        
        % Define the faces to connect each adjacent f(x) and the corresponding points at y = 0.
        q{i, j}     = (1:n_bins(i, j) - 1)';
        faces{i, j} = [q{i, j}, q{i, j} + 1, q{i, j} + n_bins(i, j) + 1, q{i, j} + n_bins(i, j)];
        
    end
end

% determine spacing between plots
plotting_space = mean(mean(cellfun(@max, ks)));
jit_width = plotting_space/8;
spacing = plotting_space * (n_series+dist_plots); % dist_plots to have a margin. set higher if you want more distance between plots
ks_offsets = (0:n_plots_per_series-1) .* spacing;

% flip so first plot in series is plotted on the *top*
ks_offsets  = fliplr(ks_offsets);

% calculate patch vertices from kernel density
verts = cell(size(data));
for i = 1:n_plots_per_series
    for j = 1:n_series
        verts{i, j} = [x{i, j}', ks{i, j}' + ks_offsets(i); x{i, j}', ones(n_bins(i, j), 1) * ks_offsets(i)];
    end
end

%% boxplot [CGA]
if box_on
    Y = cell(size(data));
    for i = 1:n_plots_per_series
        for j = 1:n_series
            quartiles   = quantile(data{i,j},[0.25 0.75 0.5]);
            iqr         = quartiles(2) - quartiles(1);
            Xs          = sort(data{i,j});
            whiskers(1) = min(Xs(Xs > (quartiles(1) - (1.5 * iqr))));
            whiskers(2) = max(Xs(Xs < (quartiles(2) + (1.5 * iqr))));
            Y{i,j}      = [quartiles whiskers];
        end
    end
end

%% jitter for the raindrops

drop_pos = cell(size(data));
for i = 1:n_plots_per_series
    for j = 1:n_series
        if box_dodge
            jit = rand(1, length(data{i,j}))*jit_width;
            offset = ks_offsets(i)-(4*j-1)*jit_width;
            drop_pos{i,j}=offset-jit;

        else
            jit = rand(1, length(data{i,j}))*jit_width;
            offset = ks_offsets(i)-(4*j-2.5)*jit_width;
            drop_pos{i,j}=offset-jit;

        end
    end
end

%% plot

hold on

% patches
offsets = zeros(n_plots_per_series,n_series);
for i = 1:n_plots_per_series
    for j = 1:n_series
        
        % plot patches
        h.p{i, j} = patch('Faces', faces{i, j}, 'Vertices', verts{i, j}, 'FaceVertexCData', colours(j, :), 'FaceColor', 'flat', 'EdgeColor', 'none', 'FaceAlpha', alpha);
        
        % scatter rainclouds
        h.s{i, j} = scatter(data{i, j}, drop_pos{i,j}, 'MarkerFaceColor', colours(j, :), 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.5, 'SizeData', raindrop_size);
    
        % Plot boxplots [CGA]
        if box_on
            if box_col_match
                bxcl = colours(j,:);
            end
            if box_dodge
                offsets(i,j) = ks_offsets(i)-(4*j-2.5)*jit_width;
            else
                offsets(i,j) = ks_offsets(i)-(4*j-2)*jit_width;
            end
                box_pos = [Y{i,j}(1) offsets(i,j)-jit_width*0.5 Y{i,j}(2)-Y{i,j}(1) jit_width];
                % mean line
                h.b1{i,j} = line([Y{i,j}(3) Y{i,j}(3)], [offsets(i,j)-jit_width*0.5 offsets(i,j)+jit_width*0.5], 'col', bxcl, 'LineWidth', line_width);
                % whiskers
                h.b2{i,j} = line([Y{i,j}(2) Y{i,j}(5)], [offsets(i,j) offsets(i,j)], 'col', bxcl, 'LineWidth', line_width);
                h.b3{i,j} = line([Y{i,j}(1) Y{i,j}(4)], [offsets(i,j) offsets(i,j)], 'col', bxcl, 'LineWidth', line_width);
                
                % 'box' of the 'boxplot'
                h.b4{i,j} = rectangle('Position', box_pos);
                set(h.b4{i,j}, 'EdgeColor', bxcl)
                set(h.b4{i,j}, 'LineWidth', line_width);
        
        end
    end
end

%% means (for mean dots)
cell_medians = cellfun(@median, data);

% plot meadian lines
if plot_median_lines
    for i = 1:n_plots_per_series - 1 % We have n_plots_per_series-1 lines because lines connect pairs of points
        for j = 1:n_series
            if box_col_match
                bxcl = colours(j,:);
            end
            % Mean line reaches the boxplot
            h.l(i, j) = line(cell_medians([i i+1], j), [offsets(i,j) offsets(i+1,j)], 'LineWidth', line_width, 'Color', bxcl);
        end
    end
end

% plot median dots
if plot_median_dots
    for i = 1:n_plots_per_series
        for j = 1:n_series
            h.m(i, j) = scatter(cell_medians(i, j), ks_offsets(i), 'MarkerFaceColor', colours(j, :), 'MarkerEdgeColor', [0 0 0], 'MarkerFaceAlpha', 1, 'SizeData', raindrop_size * 5, 'LineWidth', 2);
        end
    end
end

%% clear up axis labels

% 'YTick', likes values that *increase* as you go up the Y-axis, but we plot the first
% raincloud at the top. So flip the vector around
set(gca, 'YTick', fliplr(ks_offsets));

set(gca, 'YTickLabel', n_plots_per_series:-1:1);

%% determine plot rotation
% default option is left-to-right
% plot_top_to_bottom can be set to 1 
% NOTE: Because it's easier, we actually prepare everything plotted
% top-to-bottom, then - by default - we rotate it here. That's why the
% logical is constructed the way it is.

% rotate and flip
if plot_top_to_bottom
    view([90 -90]);
    axis ij
end
