function [varargout] = j_fit(x_data, y_data, fit_type, varargin)
%
% Fits a curve of specified type using least squares, returning
% fit parameters and displaying a plot.
% 
% [a b ...] = j_fit(x_data, y_data, fit_type, iterations=1000, scale=1, shift=0)
% [a b ...] = j_fit(x_data, y_data, @my_function, iterations=1000, scale=1, shift=0)
%
% supported fit types:
% Logistic curve (sigmoid) asymptoting at 0 & 1:    'logistic1'
% Logistic curve (sigmoid) asymptoting at 0.5 & 1:  'logistic2'
% (not implemented) Gaussian curve:                 'gaussian' (alt. 'gauss')
% First derivative of a Gaussian:                   'gauss_deriv' (alt. 's_curve')
% Custom (example function signature):              function retVal = f(x, params)
% 
% copyright 2010 Jason Fischer - University of California, Berkeley & 
%                Steve Haroz - University of California, Davis

%% optional arguments
optional_arguments_count = size(varargin,2);

%iterations
if optional_arguments_count >= 1 && isnumeric(varargin{1})
    iterations = varargin{1};
else
    iterations = 1000;
end

%start_point
if optional_arguments_count >= 2 && isnumeric(varargin{2})
    start_point = varargin{2};
else
    start_point = [0.92 0.12]; %hard coded (not rounded for testing)
end

%scale
if optional_arguments_count >= 3 && isnumeric(varargin{3})
    scale = varargin{3};
else
    scale = 1;
end

%shift
if optional_arguments_count >= 4 && isnumeric(varargin{4})
    shift = varargin{4};
else
    shift = 0;
end

%% select fit_type
fit_name = fit_type;

% if user specified function
if isa(fit_type, 'function_handle')
    f = fit_type; 
    fit_name = func2str(fit_type); % try to get readable name
else % assume string
    switch lower(fit_type)
        
      case 'logistic1' % sigmoid 1
        start_point = [1 mean(x_data)];
        f = make_scaled_shifted_logistic(scale, shift); % make the function
        
      case 'logistic2' % sigmoid 2
        scale = .5;
        shift = .5;
        start_point = [1 mean(x_data)];
        f = make_scaled_shifted_logistic(scale, shift); % make the function
        
      case {'gaussian','gauss'} %not yet added
        start_point = [0 1];
        f = @gaussian2;
        
      case {'s_curve','gauss_deriv'} % derivative of a gaussian
        c = sqrt(2)/exp(-.5); %this constant makes the a parameter correspond to the peak amplitude
        start_point = [10 0]; %arbitrary hard coded -- fix!!
        f = make_gaussian_derivative(c); % make the function
        
      otherwise
        error('Unknown fit type.')
    end % switch
end % if

%% main estimation
[boot_params, y_predicted, x_data, y_data, fit_parameters] = wrapper(f, start_point, x_data, y_data, iterations);

%find bootstrapped thresholds for error bars on plot
y_thr1 = (min(y_data) + .5)/2; %these values can be changed to select where thresholds are computed
thr1 = (log((scale/(y_thr1-shift))-1)-boot_params(:,1).*boot_params(:,2))./-boot_params(:,1);
y_thr2 = .5;
thr2 = (log((scale/(y_thr2-shift))-1)-boot_params(:,1).*boot_params(:,2))./-boot_params(:,1);
y_thr3 = (max(y_data) + .5)/2;
thr3 = (log((scale/(y_thr3-shift))-1)-boot_params(:,1).*boot_params(:,2))./-boot_params(:,1);

%% construct plot
% this currently works best w/ a single screen setup
screen_size = get(0,'ScreenSize');
fig = figure('Position',[screen_size(3)/6 screen_size(4)/12 4*screen_size(3)/6 10*screen_size(4)/12]);

% main data plot
subplot(5,5,1:15)
plot_main(x_data,y_data,y_predicted,[.05 .5 .8],y_thr1,y_thr2,y_thr3,thr1,thr2,thr3,fit_name);

% bootstrapped parameter plot
subplot(5,5,[16 17 21 22])
boot_scatter(boot_params(:,1),boot_params(:,2))

% output parameters
a_boot_std = std(boot_params(:,1));
b_boot_std = std(boot_params(:,2));

param_list{1} = 'Parameters:';
for i = 1:size(fit_parameters, 2)
    line_num = (i-1)*3 + 2;
    param_list{line_num + 0} = ' ';
    param_list{line_num + 1} = ['p' num2str(i) ': ' num2str(fit_parameters(:,i))];
    param_list{line_num + 2} = ['+/- ' num2str(a_boot_std) ' se'];
end

write_params(param_list,1) %second argument specifies which parameter box to print to

thresh_list{1} = 'Thresholds:';
thresh_list{2} = ' ';
thresh_list{3} = [num2str(y_thr1) ': ' num2str(mean(thr1))];
thresh_list{4} = ['+/- ' num2str(std(thr1)) ' se'];
thresh_list{5} = ' ';
thresh_list{6} = [num2str(y_thr2) ': ' num2str(mean(thr2))];
thresh_list{7} = ['+/- ' num2str(std(thr2)) ' se'];
thresh_list{8} = ' ';
thresh_list{9} = [num2str(y_thr3) ': ' num2str(mean(thr3))];
thresh_list{10} = ['+/- ' num2str(std(thr3)) ' se'];

write_params(thresh_list,2)

    function [] = plot_main(x_data,y_data,y_predicted,color,y_thr1,y_thr2,y_thr3,thr1,thr2,thr3,fit_name)
        plot(x_data,y_data,'MarkerEdgeColor',color,'MarkerFaceColor',color,'Marker','o','LineStyle','none')
        %axis([xmin xmax ymin ymax]);
        set(gca,'FontSize',14);
        hold on;
        plot((min(x_data):(max(x_data)-min(x_data))/100:max(x_data)),y_predicted,'Color',[.3 .3 .3],'LineStyle','-','LineWidth',2);
        new_title = strrep(fit_name, '_', ' ');
        title([new_title ' fit']);
        if exist('thr1','var') %allow for some fit types to skip error bars at the moment
            %line([mean(thr1)-std(thr1) mean(thr1)+std(thr1)],[y_thr1 y_thr1],'LineWidth',1,'Color',[.6 .6 .6]);
            %line([mean(thr2)-std(thr2) mean(thr2)+std(thr2)],[y_thr2 y_thr2],'LineWidth',1,'Color',[.6 .6 .6]);
           % line([mean(thr3)-std(thr3) mean(thr3)+std(thr3)],[y_thr3 y_thr3],'LineWidth',1,'Color',[.6 .6 .6]);
        end
    end

    function [] = boot_scatter(boot_params1,boot_params2)

        a_boot_mean = mean(boot_params1);
        a_boot_std = std(boot_params1);
        b_boot_mean = mean(boot_params2);
        b_boot_std = std(boot_params2);

        plot(boot_params1,boot_params2,'Color',[.8 .4 .8],'Marker','o','LineStyle','none'); hold on;
        plot(a_boot_mean,b_boot_mean,'MarkerEdgeColor',[.05 .5 .8],'MarkerFaceColor',[.05 .5 .8],'Marker','o','LineStyle','none');
        axis([a_boot_mean-3*a_boot_std a_boot_mean+3*a_boot_std b_boot_mean-3*b_boot_std b_boot_mean+3*b_boot_std])
        line([a_boot_mean-3*a_boot_std a_boot_mean+3*a_boot_std],[b_boot_mean-b_boot_std b_boot_mean-b_boot_std],'LineWidth',1,'Color',[.7 .7 .7])
        line([a_boot_mean-3*a_boot_std a_boot_mean+3*a_boot_std],[b_boot_mean-2*b_boot_std b_boot_mean-2*b_boot_std],'LineWidth',1,'Color',[.7 .7 .7])
        line([a_boot_mean-3*a_boot_std a_boot_mean+3*a_boot_std],[b_boot_mean+b_boot_std b_boot_mean+b_boot_std],'LineWidth',1,'Color',[.7 .7 .7])
        line([a_boot_mean-3*a_boot_std a_boot_mean+3*a_boot_std],[b_boot_mean+2*b_boot_std b_boot_mean+2*b_boot_std],'LineWidth',1,'Color',[.7 .7 .7])
        line([a_boot_mean-a_boot_std a_boot_mean-a_boot_std],[b_boot_mean-3*b_boot_std b_boot_mean+3*b_boot_std],'LineWidth',1,'Color',[.7 .7 .7])
        line([a_boot_mean-2*a_boot_std a_boot_mean-2*a_boot_std],[b_boot_mean-3*b_boot_std b_boot_mean+3*b_boot_std],'LineWidth',1,'Color',[.7 .7 .7])
        line([a_boot_mean+a_boot_std a_boot_mean+a_boot_std],[b_boot_mean-3*b_boot_std b_boot_mean+3*b_boot_std],'LineWidth',1,'Color',[.7 .7 .7])
        line([a_boot_mean+2*a_boot_std a_boot_mean+2*a_boot_std],[b_boot_mean-3*b_boot_std b_boot_mean+3*b_boot_std],'LineWidth',1,'Color',[.7 .7 .7])
        line([a_boot_mean a_boot_mean],[b_boot_mean-3*b_boot_std b_boot_mean+3*b_boot_std],'LineWidth',1,'Color',[.05 .5 .8])
        line([a_boot_mean-3*a_boot_std a_boot_mean+3*a_boot_std],[b_boot_mean b_boot_mean],'LineWidth',1,'Color',[.05 .5 .8])
        set(gca,'FontSize',14);
        %title('Bootstrapped parameter estimates')
        xlabel('parameter 1')
        ylabel('parameter 2')
   
    end

    function [] = write_params(param_list,position)
        uicontrol('Parent',fig,'Units','normalized','Position',[0.46 + (.24*(position-1)),0.1,0.2,0.31],'Style','edit','Max',100,'String',param_list,'Background','white','HorizontalAlignment','left')
    end


%% assign outputs
for arg = 1:length(fit_parameters)
    varargout(arg) = {fit_parameters(arg)};
end

end %j_fit main function
%----------------------------------------------------------------------

%% Core of curve fitting
% wraps the calculation process ---------------------------------------
function [boot_params, y_predicted, x_data, y_data, fit_parameters] = wrapper(f, start_point, x_data, y_data, iterations)
    fit_curve = make_fit_f_curve(f, start_point);   % create a function with a preset start value
    boot_params = bootstrp(iterations,fit_curve,x_data,y_data);
    fit_parameters = fit_curve(x_data,y_data);
    new_x = (min(x_data):(max(x_data)-min(x_data))/100:max(x_data)); % what does this do???
    y_predicted = f(new_x, fit_parameters);
    [x_data y_data] = percent_corr(x_data,y_data); %convert to % correct if necessary 
end

% creates a function estimator using 'f' ------------------------------
function [newF] = make_fit_f_curve(f, start_point)
    newF = @(x,y) fit_f_curve(f, start_point, x, y);
end

%s curve (1st derivative of a gaussian) -------------------------------
function [estimates, sum_sqr_err] = fit_f_curve(f, start_point, x, y)
% passed in: start_point = array ( size(start_point) = size(params)

    [x y] = percent_corr(x,y); %convert to % correct if necessary 

    estimates = fminsearch(@curve, start_point);
    sum_sqr_err = curve(estimates);

    function [sse] = curve(params)
        FittedCurve = f(x, params);
        ErrorVector = FittedCurve - y;
        sse = sum(ErrorVector .^ 2);
    end
end

%% support functions
function [x_adjusted y_adjusted] = percent_corr(x, y)
%convert to proportion correct if the data consists of only 2 possible
%values (correct/incorrect).    
    if size(unique(y),1) == 2
        x_set = unique(x);
        y_set = zeros(length(x_set),1);
        for i = 1:length(x_set);
            y_set(i) = mean(y(x == x_set(i))); %compute proportion correct at each unique level of x
        end
        x_adjusted = x_set;
        y_adjusted = y_set;
    else
        x_adjusted = x;
        y_adjusted = y;
    end

end

%% Functions to fit

% gaussian function
function retVal = gaussian (x, a, b, c)
    retVal = a*exp( (x-b).^2 / -2*a*c*c );
end
% gaussian function
function retVal = gaussian2 (x, params)
    mean = params(1);
    sigma = params(2);
    a = 1 / (sigma*sqrt(2*pi));
    b = mean;
    c = sigma;
    retVal = gaussian (x, a, b, c);
end

% 1st derivative of a gaussian
function retVal = gaussian_derivative (x, a, b, c)
    retVal = a*b*c*x.*exp(-((b*x).^2));
end
% make function of 1st derivative of a gaussian (with fixed c value)
function f_out = make_gaussian_derivative(c)
    f_out = @(x,params) gaussian_derivative(x, params(1), params(2), c);
end

% sigmoid curve (psychometric function)
function retVal = logistic (x, scale, shift, a, b)
    retVal = scale./(1+exp(-a*(x-b)))+shift;
end
% make function of sigmoid curve with fixed scale and shift
function f_out = make_scaled_shifted_logistic(scale, shift)
    f_out = @(x,params) logistic(x, scale, shift, params(1), params(2));
end
%----------------------------------------------------------------------

