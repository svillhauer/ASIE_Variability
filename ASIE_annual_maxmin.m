%% Load ASIE data 
clear all 
cd   '/Users/sarahvillhauer/Desktop/New ASIE/New Data'
addpath '/Users/sarahvillhauer/Desktop/New ASIE/Paper scripts'

totsieload = readtable('totalASIE.csv');
regsieload = readtable('regionalASIE.csv');


% Missing data is represented as NA, replace with NaNs
for varName = totsieload.Properties.VariableNames
    col = totsieload.(varName{1});
    
    if isnumeric(col) || islogical(col)
        % Replace NaN in numeric or logical arrays (if needed)
        col(isnan(col)) = NaN; % Already NaN, just for completeness
    elseif iscellstr(col)
        % Replace 'NA' with NaN in cell array of strings
        col(strcmp(col, 'NA')) = {NaN};
        % Convert cell array back to the table variable
        totsieload.(varName{1}) = col;
    end
end


totsie = table2cell(totsieload(:,2));  

[numRows, numCols] = size(totsie);
for i = 1:numRows
    for j = 1:numCols
        % Remove single quotes from each cell
        totsie{i,j} = strrep(totsie{i,j}, '''', '');  % Replace single quotes with nothing
    end
end

% Step 3: Convert cell array to numeric array (assuming all elements are numeric)
totsie = str2double(totsie);  % Convert cell array of strings to numeric array




% Missing data is represented as NA, replace with NaNs
for varName = regsieload.Properties.VariableNames
    col = regsieload.(varName{1});
    
    if isnumeric(col) || islogical(col)
        % Replace NaN in numeric or logical arrays (if needed)
        col(isnan(col)) = NaN; % Already NaN, just for completeness
    elseif iscellstr(col)
        % Replace 'NA' with NaN in cell array of strings
        col(strcmp(col, 'NA')) = {NaN};
        % Convert cell array back to the table variable
        regsieload.(varName{1}) = col;
    end
end

regsie = table2cell(regsieload(:,2:6));  

[numRows, numCols] = size(regsie);
for i = 1:numRows
    for j = 1:numCols
        % Remove single quotes from each cell
        regsie{i,j} = strrep(regsie{i,j}, '''', '');  % Replace single quotes with nothing
    end
end

% Step 3: Convert cell array to numeric array (assuming all elements are numeric)
regsie = str2double(regsie);  % Convert cell array of strings to numeric array



% Time
% 300:16928
time = regsieload(300:16928,1);
time = table2array(time);

%% Interpolate SIE Data 
% Column 1 = King Hakon, Column 2 = Ross, Column 3 = East Antarctica, Column 4 = Weddell, Column 5 = Bellingshausen  
% Filling in NaNs
khsie = regsie(300:16928,1); 
%khsiei = fillmissing(khsie, 'spline');

rosssie = regsie(300:16928,2);
%rosssiei = fillmissing(rosssie, 'spline');

easie = regsie(300:16928,3);
%easiei = fillmissing(easie, 'spline');

wedsie = regsie(300:16928,4);
%wedsiei = fillmissing(wedsie, 'spline');

bellsie = regsie(300:16928,5);
%bellsiei = fillmissing(bellsie, 'spline');

totsie = totsie(300:16928,1);

% Reorganize data, where rows are day-of-year s and columns are year (e.g.
% column 41 = 2019)
% Since 1978 and 2024 don't have complete data, define those years
% separately, then concatenate 
% 1979 starts at index 67 
% 2024 starts at index 16492
ross_1978 = nan(1, 365);
kh_1978 = nan(1, 365);
ea_1978 = nan(1, 365);
wed_1978  = nan(1, 365);
bell_1978  = nan(1, 365);
tot_1978 = nan(1, 365);
time_1978 = nan(1, 365);


ross_1978(365-65:365) = rosssie(1:66); 
kh_1978(365-65:365) = khsie(1:66); 
ea_1978(365-65:365) = easie(1:66); 
wed_1978(365-65:365) = wedsie(1:66); 
bell_1978(365-65:365) = bellsie(1:66); 
tot_1978(365-65:365) = totsie(1:66); 
time_1978(365-65:365) = time(1:66);

ross_2024 = nan(1, 365);
kh_2024 = nan(1, 365);
ea_2024 = nan(1, 365);
wed_2024  = nan(1, 365);
bell_2024  = nan(1, 365);
tot_2024 = nan(1, 365);
time_2024 = nan(1, 365);


ross_2024(1:138) = rosssie(16492:16629); 
kh_2024(1:138) = khsie(16492:16629); 
ea_2024(1:138) = easie(16492:16629); 
wed_2024(1:138) = wedsie(16492:16629); 
bell_2024(1:138) = bellsie(16492:16629); 
tot_2024(1:138) = totsie(16492:16629); 
time_2024(1:138) = time(16492:16629); 





timereshaped = reshape(time(67:16491), 365, 45);
rossreshaped = reshape(rosssie(67:16491), 365, 45);
khreshaped = reshape(khsie(67:16491), 365, 45);
eareshaped = reshape(easie(67:16491), 365, 45);
wedreshaped = reshape(wedsie(67:16491), 365, 45);
bellreshaped = reshape(bellsie(67:16491), 365, 45);
totreshaped = reshape(totsie(67:16491), 365,45); 


% Concatenate 1978 and 2024 to reshaped data 
rossreshaped = [ross_1978' rossreshaped ross_2024']; 
khreshaped = [kh_1978' khreshaped kh_2024']; 
eareshaped = [ea_1978' eareshaped ea_2024']; 
wedreshaped = [wed_1978' wedreshaped wed_2024']; 
bellreshaped = [bell_1978' bellreshaped bell_2024']; 
totreshaped = [tot_1978' totreshaped tot_2024']; 
timereshaped = [time_1978' timereshaped time_2024']; 




%% Extract the maximum and minimum 
% Ross
max_ross = max(rossreshaped(:,2:46)); 
min_ross = min(rossreshaped(:,2:47)); 

% King Hakon
max_kh = max(khreshaped(:,2:46)); 
min_kh = min(khreshaped(:,2:47)); 

% East Antarctica
max_ea = max(eareshaped(:,2:46)); 
min_ea = min(eareshaped(:,2:47)); 

% Weddell 
max_wed = max(wedreshaped(:,2:46)); 
min_wed = min(wedreshaped(:,2:47)); 

% Bellingshausen 
max_bell = max(bellreshaped(:,2:46)); 
min_bell = min(bellreshaped(:,2:47)); 

% Total
max_tot = max(totreshaped(:,2:46)); 
min_tot = min(totreshaped(:,2:47)); 


years_max = [1979:2023]; 
years_min = [1979:2024]; 
%% Calculate the linear trends of the annual maxima and minima 
% Maxima 
% Ross
p_maxross = polyfit(1:length(years_max), max_ross, 1);
a_post = p_maxross(1); % Slope (trend)
b_post = p_maxross(2); % Intercept
[b_post, bint_post, r_post, rint_post, stats_post] = regress(max_ross', [ones(length(years_max),1), (1:length(years_max))']);
p_value_maxross = stats_post(3);
SSR_post = sum(r_post.^2);
SST_post = sum((max_ross - mean(max_ross)).^2);
R_squared_maxross = 1 - (SSR_post / SST_post);


a_post_maxross = b_post(2); % Slope (trend)
slope_confidence_interval_maxross = bint_post(2,:); % Confidence interval for the slope
% Calculate the margin of error for the slope
slope_margin_of_error_maxross = (slope_confidence_interval_maxross(2) - slope_confidence_interval_maxross(1)) / 2;

% King Hakon 
p_maxkh = polyfit(1:length(years_max), max_kh, 1);
a_post = p_maxkh(1); % Slope (trend)
b_post = p_maxkh(2); % Intercept
[b_post, bint_post, r_post, rint_post, stats_post] = regress(max_kh', [ones(length(years_max),1), (1:length(years_max))']);
p_value_maxkh = stats_post(3);
SSR_post = sum(r_post.^2);
SST_post = sum((max_kh - mean(max_kh)).^2);
R_squared_maxkh = 1 - (SSR_post / SST_post);


a_post_maxkh = b_post(2); % Slope (trend)
slope_confidence_interval_maxkh = bint_post(2,:); % Confidence interval for the slope
% Calculate the margin of error for the slope
slope_margin_of_error_maxkh = (slope_confidence_interval_maxkh(2) - slope_confidence_interval_maxkh(1)) / 2;

% East Antarctica
p_maxea = polyfit(1:length(years_max), max_ea, 1);
a_post = p_maxea(1); % Slope (trend)
b_post = p_maxea(2); % Intercept
[b_post, bint_post, r_post, rint_post, stats_post] = regress(max_ea', [ones(length(years_max),1), (1:length(years_max))']);
p_value_maxea = stats_post(3);
SSR_post = sum(r_post.^2);
SST_post = sum((max_ea - mean(max_ea)).^2);
R_squared_maxea = 1 - (SSR_post / SST_post);

a_post_maxea = b_post(2); % Slope (trend)
slope_confidence_interval_maxea = bint_post(2,:); % Confidence interval for the slope
% Calculate the margin of error for the slope
slope_margin_of_error_maxea = (slope_confidence_interval_maxea(2) - slope_confidence_interval_maxea(1)) / 2;

% Weddell 
p_maxwed = polyfit(1:length(years_max), max_wed, 1);
a_post = p_maxwed(1); % Slope (trend)
b_post = p_maxwed(2); % Intercept
[b_post, bint_post, r_post, rint_post, stats_post] = regress(max_wed', [ones(length(years_max),1), (1:length(years_max))']);
p_value_maxwed = stats_post(3);
SSR_post = sum(r_post.^2);
SST_post = sum((max_wed - mean(max_wed)).^2);
R_squared_maxwed = 1 - (SSR_post / SST_post);


a_post_maxwed = b_post(2); % Slope (trend)
slope_confidence_interval_maxwed = bint_post(2,:); % Confidence interval for the slope
% Calculate the margin of error for the slope
slope_margin_of_error_maxwed = (slope_confidence_interval_maxwed(2) - slope_confidence_interval_maxwed(1)) / 2;

% Bellingshausen
p_maxbell = polyfit(1:length(years_max), max_bell, 1);
a_post = p_maxbell(1); % Slope (trend)
b_post = p_maxbell(2); % Intercept
[b_post, bint_post, r_post, rint_post, stats_post] = regress(max_bell', [ones(length(years_max),1), (1:length(years_max))']);
p_value_maxbell = stats_post(3);
SSR_post = sum(r_post.^2);
SST_post = sum((max_bell - mean(max_bell)).^2);
R_squared_maxbell = 1 - (SSR_post / SST_post);


a_post_maxbell = b_post(2); % Slope (trend)
slope_confidence_interval_maxbell = bint_post(2,:); % Confidence interval for the slope
% Calculate the margin of error for the slope
slope_margin_of_error_maxbell = (slope_confidence_interval_maxbell(2) - slope_confidence_interval_maxbell(1)) / 2;

% Total
p_maxtot = polyfit(1:length(years_max), max_tot, 1);
a_post = p_maxtot(1); % Slope (trend)
b_post = p_maxtot(2); % Intercept
[b_post, bint_post, r_post, rint_post, stats_post] = regress(max_tot', [ones(length(years_max),1), (1:length(years_max))']);
p_value_maxtot = stats_post(3);
SSR_post = sum(r_post.^2);
SST_post = sum((max_tot - mean(max_tot)).^2);
R_squared_maxtot = 1 - (SSR_post / SST_post);


a_post_maxtot = b_post(2); % Slope (trend)
slope_confidence_interval_maxtot = bint_post(2,:); % Confidence interval for the slope
% Calculate the margin of error for the slope
slope_margin_of_error_maxtot = (slope_confidence_interval_maxtot(2) - slope_confidence_interval_maxtot(1)) / 2;

% Minima
% Ross
p_minross = polyfit(1:length(years_min), min_ross, 1);
a_post = p_minross(1); % Slope (trend)
b_post = p_minross(2); % Intercept
[b_post, bint_post, r_post, rint_post, stats_post] = regress(min_ross', [ones(length(years_min),1), (1:length(years_min))']);
p_value_minross = stats_post(3);
SSR_post = sum(r_post.^2);
SST_post = sum((min_ross - mean(min_ross)).^2);
R_squared_minross = 1 - (SSR_post / SST_post);


a_post_minross = b_post(2); % Slope (trend)
slope_confidence_interval_minross = bint_post(2,:); % Confidence interval for the slope
% Calculate the margin of error for the slope
slope_margin_of_error_minross = (slope_confidence_interval_minross(2) - slope_confidence_interval_minross(1)) / 2;

% King Hakon 
p_minkh = polyfit(1:length(years_min), min_kh, 1);
a_post = p_minkh(1); % Slope (trend)
b_post = p_minkh(2); % Intercept
[b_post, bint_post, r_post, rint_post, stats_post] = regress(min_kh', [ones(length(years_min),1), (1:length(years_min))']);
p_value_minkh = stats_post(3);
SSR_post = sum(r_post.^2);
SST_post = sum((min_kh - mean(min_kh)).^2);
R_squared_minkh = 1 - (SSR_post / SST_post);


a_post_minkh = b_post(2); % Slope (trend)
slope_confidence_interval_minkh = bint_post(2,:); % Confidence interval for the slope
% Calculate the margin of error for the slope
slope_margin_of_error_minkh = (slope_confidence_interval_minkh(2) - slope_confidence_interval_minkh(1)) / 2;


% East Antarctica
p_minea = polyfit(1:length(years_min), min_ea, 1);
a_post = p_minea(1); % Slope (trend)
b_post = p_minea(2); % Intercept
[b_post, bint_post, r_post, rint_post, stats_post] = regress(min_ea', [ones(length(years_min),1), (1:length(years_min))']);
p_value_minea = stats_post(3);
SSR_post = sum(r_post.^2);
SST_post = sum((min_ea - mean(min_ea)).^2);
R_squared_minea = 1 - (SSR_post / SST_post);


a_post_minea = b_post(2); % Slope (trend)
slope_confidence_interval_minea = bint_post(2,:); % Confidence interval for the slope
% Calculate the margin of error for the slope
slope_margin_of_error_minea = (slope_confidence_interval_minea(2) - slope_confidence_interval_minea(1)) / 2;


% Weddell 
p_minwed = polyfit(1:length(years_min), min_wed, 1);
a_post = p_minwed(1); % Slope (trend)
b_post = p_minwed(2); % Intercept
[b_post, bint_post, r_post, rint_post, stats_post] = regress(min_wed', [ones(length(years_min),1), (1:length(years_min))']);
p_value_minwed = stats_post(3);
SSR_post = sum(r_post.^2);
SST_post = sum((min_wed - mean(min_wed)).^2);
R_squared_minwed = 1 - (SSR_post / SST_post);


a_post_minwed = b_post(2); % Slope (trend)
slope_confidence_interval_minwed = bint_post(2,:); % Confidence interval for the slope
% Calculate the margin of error for the slope
slope_margin_of_error_minwed = (slope_confidence_interval_minwed(2) - slope_confidence_interval_minwed(1)) / 2;


% Bellingshausen
p_minbell = polyfit(1:length(years_min), min_bell, 1);
a_post = p_minbell(1); % Slope (trend)
b_post = p_minbell(2); % Intercept
[b_post, bint_post, r_post, rint_post, stats_post] = regress(min_bell', [ones(length(years_min),1), (1:length(years_min))']);
p_value_minbell = stats_post(3);
SSR_post = sum(r_post.^2);
SST_post = sum((min_bell - mean(min_bell)).^2);
R_squared_minbell = 1 - (SSR_post / SST_post);


a_post_minbell = b_post(2); % Slope (trend)
slope_confidence_interval_minbell = bint_post(2,:); % Confidence interval for the slope
% Calculate the margin of error for the slope
slope_margin_of_error_minbell = (slope_confidence_interval_minbell(2) - slope_confidence_interval_minbell(1)) / 2;


% Total
p_mintot = polyfit(1:length(years_min), min_tot, 1);
a_post = p_mintot(1); % Slope (trend)
b_post = p_mintot(2); % Intercept
[b_post, bint_post, r_post, rint_post, stats_post] = regress(min_tot', [ones(length(years_min),1), (1:length(years_min))']);
p_value_mintot = stats_post(3);
SSR_post = sum(r_post.^2);
SST_post = sum((min_tot - mean(min_tot)).^2);
R_squared_mintot = 1 - (SSR_post / SST_post);


a_post_mintot = b_post(2); % Slope (trend)
slope_confidence_interval_mintot = bint_post(2,:); % Confidence interval for the slope
% Calculate the margin of error for the slope
slope_margin_of_error_mintot = (slope_confidence_interval_mintot(2) - slope_confidence_interval_mintot(1)) / 2;




x_values = 1:length(years_max);
y_values_fit_maxross = polyval(p_maxross, x_values);
y_values_fit_maxkh = polyval(p_maxkh, x_values);
y_values_fit_maxea = polyval(p_maxea, x_values);
y_values_fit_maxwed = polyval(p_maxwed, x_values);
y_values_fit_maxbell = polyval(p_maxbell, x_values);
y_values_fit_maxtot = polyval(p_maxtot, x_values);

x_values = 1:length(years_min);
y_values_fit_minross = polyval(p_minross, x_values);
y_values_fit_minkh = polyval(p_minkh, x_values);
y_values_fit_minea = polyval(p_minea, x_values);
y_values_fit_minwed = polyval(p_minwed, x_values);
y_values_fit_minbell = polyval(p_minbell, x_values);
y_values_fit_mintot = polyval(p_mintot, x_values);


%% Correlate across regions 
% Maximum correlations 
[r_rosskh_max, p_rosskh_max] = corrcoef(max_ross,max_kh);
[r_rossea_max, p_rossea_max] = corrcoef(max_ross,max_ea);
[r_rosswed_max, p_rosswed_max] = corrcoef(max_ross,max_wed);
[r_rossbell_max, p_rossbell_max] = corrcoef(max_ross,max_bell);

[r_khea_max, p_khea_max] = corrcoef(max_kh,max_ea);
[r_khwed_max, p_khwed_max] = corrcoef(max_kh,max_wed);
[r_khbell_max, p_khbell_max] = corrcoef(max_kh,max_bell);

[r_eawed_max, p_eawed_max] = corrcoef(max_ea,max_wed);
[r_eabell_max, p_eabell_max] = corrcoef(max_ea,max_bell);

[r_wedbell_max, p_wedbell_max] = corrcoef(max_wed,max_bell);


[r_totross_max, p_totross_max] = corrcoef(max_tot,max_ross);
[r_totkh_max, p_totkh_max] = corrcoef(max_tot,max_kh);
[r_totea_max, p_totea_max] = corrcoef(max_tot,max_ea);
[r_totwed_max, p_totwed_max] = corrcoef(max_tot,max_wed);
[r_totbell_max, p_totbell_max] = corrcoef(max_tot,max_bell);


% Minimum correlations
[r_rosskh_min, p_rosskh_min] = corrcoef(min_ross,min_kh);
[r_rossea_min, p_rossea_min] = corrcoef(min_ross,min_ea);
[r_rosswed_min, p_rosswed_min] = corrcoef(min_ross,min_wed);
[r_rossbell_min, p_rossbell_min] = corrcoef(min_ross,min_bell);

[r_khea_min, p_khea_min] = corrcoef(min_kh,min_ea);
[r_khwed_min, p_khwed_min] = corrcoef(min_kh,min_wed);
[r_khbell_min, p_khbell_min] = corrcoef(min_kh,min_bell);

[r_eawed_min, p_eawed_min] = corrcoef(min_ea,min_wed);
[r_eabell_min, p_eabell_min] = corrcoef(min_ea,min_bell);

[r_wedbell_min, p_wedbell_min] = corrcoef(min_wed,min_bell);


[r_totross_min, p_totross_min] = corrcoef(min_tot,min_ross);
[r_totkh_min, p_totkh_min] = corrcoef(min_tot,min_kh);
[r_totea_min, p_totea_min] = corrcoef(min_tot,min_ea);
[r_totwed_min, p_totwed_min] = corrcoef(min_tot,min_wed);
[r_totbell_min, p_totbell_min] = corrcoef(min_tot,min_bell);


%% Plot the maximum



plot(1) % Subplot of annual maxima for [1] five regions and [2] entire continent
subplot(2,1,1)
plot(years_max,max_ross,'linewidth',4,'color','red', 'LineStyle', '-', 'marker','none')
hold on
plot(years_max,max_bell,'linewidth',4,'color','#7E2F8E', 'LineStyle', '-', 'marker','none')
hold on
plot(years_max,max_wed,'linewidth',4,'color','magenta', 'LineStyle', '-', 'marker','none')
hold on
plot(years_max,max_kh,'linewidth',4,'color','#77AC30', 'LineStyle', '-', 'marker','none')
hold on
plot(years_max,max_ea,'linewidth',4,'color', 'blue', 'LineStyle', '-', 'marker','none')

ylabel('Sea Ice Extent [$10^6$ km$^2$]', 'Fontsize', 14, 'interpreter','latex')
%xlabel('Year', 'Fontsize', 14, 'interpreter', 'latex')
set(gca, 'FontSize', 20)
set(gca,'XTickLabel',[])
set(gcf,'color','white')
ylim([0 7])
xlim([1979 2023])
grid on

title('Maximum Antarctic Sea Ice Extent', 'interpreter', 'latex')
legend('Ross','Amundsen-Bellingshausen','Weddell','King Hakon','East Antarctica', 'fontsize',12,'interpreter','latex')


subplot(2,1,2)
plot(years_max,max_tot,'linewidth',4,'color','black', 'LineStyle', '-', 'marker','none')

ylabel('Sea Ice Extent [$10^6$ km$^2$]', 'Fontsize', 14, 'interpreter','latex')
xlabel('Year', 'Fontsize', 14, 'interpreter', 'latex')
set(gca, 'FontSize', 20)
%set(gca,'YTickLabel',[])
set(gcf,'color','white')
%ylim([14.5 17.5])
xlim([1979 2023])
grid on

%title('Retreat Amplitude', 'interpreter', 'latex')
legend('Total','fontsize',12,'interpreter','latex')








    plot(3) % Subplot of annual minima for [1] five regions and [2] entire continent
    subplot(2,1,1)
    plot(years_min,min_ross,'linewidth',4,'color','red', 'LineStyle', '-', 'marker','none')
    hold on
    plot(years_min,min_bell,'linewidth',4,'color','#7E2F8E', 'LineStyle', '-', 'marker','none')
    hold on
    plot(years_min,min_wed,'linewidth',4,'color','magenta', 'LineStyle', '-', 'marker','none')
    hold on
    plot(years_min,min_kh,'linewidth',4,'color','#77AC30', 'LineStyle', '-', 'marker','none')
    hold on
    plot(years_min,min_ea,'linewidth',4,'color', 'blue', 'LineStyle', '-', 'marker','none')
    
   


    ylabel('Sea Ice Extent [$10^6$ km$^2$]', 'Fontsize', 14, 'interpreter','latex')
    %xlabel('Year', 'Fontsize', 14, 'interpreter', 'latex')
    set(gca, 'FontSize', 20)
    set(gca,'XTickLabel',[])
    set(gcf,'color','white')
    ylim([0 2])
    xlim([1979 2024])
    grid on

    title('Minimum Antarctic Sea Ice Extent', 'interpreter', 'latex')
    legend('Ross','Amundsen-Bellingshausen','Weddell','King Hakon','East Antarctica', 'fontsize',12,'interpreter','latex')


    subplot(2,1,2)
    plot(years_min,min_tot,'linewidth',4,'color','black', 'LineStyle', '-', 'marker','none')

    ylabel('Sea Ice Extent [$10^6$ km$^2$]', 'Fontsize', 14, 'interpreter','latex')
    xlabel('Year', 'Fontsize', 14, 'interpreter', 'latex')
    set(gca, 'FontSize', 20)
    %set(gca,'YTickLabel',[])
    set(gcf,'color','white')
    %ylim([14.5 17.5])
    xlim([1979 2024])
    grid on

    %title('Retreat Amplitude', 'interpreter', 'latex')
    legend('Total','fontsize',12,'interpreter','latex')





    plot(7)
    subplot(2,3,1)
    plot(years_min,min_ross,'linewidth',4,'color','red', 'LineStyle', '-', 'marker','none')
    hold on
    plot(years_min,y_values_fit_minross,'linewidth',4,'color','red', 'LineStyle', '--', 'marker','none' )
    title('Ross Minimum', 'interpreter', 'latex')
    ylim([0 2])
    set(gca, 'xticklabel', {})
    ylabel('Regional Sea Ice Extent [$10^6$ km$^2$]', 'Fontsize', 14,'interpreter','latex')
    set(gca, 'FontSize', 20)
    set(gcf,'color','white')
    xlim([1979 2024])
    grid on

    subplot(2,3,2)
    plot(years_min,min_bell,'linewidth',4,'color','#7E2F8E', 'LineStyle', '-', 'marker','none')
    hold on
    plot(years_min,y_values_fit_minbell,'linewidth',4,'color','#7E2F8E', 'LineStyle', '--', 'marker','none' )


    title('Amundsen-Bellingshausen Minimum', 'interpreter', 'latex')
    %legend('minimum','Minimum', 'fontsize', 12,'interpreter','latex')
    ylim([0 2])
    set(gca, 'xticklabel', {})
    set(gca, 'yticklabel', {})
    %ylabel('Regional Sea Ice Extent ', 'Fontsize', 14,'interpreter','latex')
    set(gca, 'FontSize', 20)
    set(gcf,'color','white')
    xlim([1979 2024])
    grid on

    subplot(2,3,3)
    plot(years_min,min_wed,'linewidth',4,'color','magenta', 'LineStyle', '-', 'marker','none')
    hold on
    plot(years_min,y_values_fit_minwed,'linewidth',4,'color','magenta', 'LineStyle', '--', 'marker','none' )


    title('Weddell Minimum', 'interpreter', 'latex')
    %legend('minimum','Minimum', 'fontsize', 12,'interpreter','latex')
    ylim([0 2])
    set(gca, 'xticklabel', {})
    set(gca, 'yticklabel', {})

    set(gca, 'FontSize', 20)
    set(gcf,'color','white')
    xlim([1979 2024])
    grid on




    subplot(2,3,4)
    plot(years_min,min_kh,'linewidth',4,'color','#77AC30', 'LineStyle', '-', 'marker','none')
    hold on
    plot(years_min,y_values_fit_minkh,'linewidth',4,'color','#77AC30', 'LineStyle', '--', 'marker','none' )

    title('King Hakon Minimum', 'interpreter', 'latex')
    %legend('minimum','Minimum', 'fontsize', 12,'interpreter','latex')
    ylim([0 2])
    %ylabel('Regional Sea Ice Extent ', 'Fontsize', 14,'interpreter','latex')
    xlabel('Year', 'Fontsize', 14, 'interpreter', 'latex')
    ylabel('Regional Sea Ice Extent [$10^6$ km$^2$] ', 'Fontsize', 14,'interpreter','latex')
    set(gca, 'FontSize', 20)
    set(gcf,'color','white')
    xlim([1979 2024])
    grid on

    subplot(2,3,5)
    plot(years_min,min_ea,'linewidth',4,'color', 'blue', 'LineStyle', '-', 'marker','none')
    hold on
    plot(years_min,y_values_fit_minea,'linewidth',4,'color','blue', 'LineStyle', '--', 'marker','none' )


    title('East Antarctica Minimum', 'interpreter', 'latex')
    %legend('minimum','Minimum', 'fontsize', 12,'interpreter','latex')
    ylim([0 2])
    xlabel('Year', 'Fontsize', 14, 'interpreter', 'latex')
    set(gca, 'yticklabel', {})
    %ylabel('Sea Ice Extent ', 'Fontsize', 14,'interpreter','latex')
    set(gca, 'FontSize', 20)
    set(gcf,'color','white')
    xlim([1979 2024])
    grid on


    subplot(2,3,6)
    yyaxis left



    plot(years_min,min_tot,'linewidth',4,'color','black', 'LineStyle', '-', 'marker','none')
    hold on
    plot(years_min,y_values_fit_mintot,'linewidth',4,'color','black', 'LineStyle', '--', 'marker','none' )
    ylabel('Total Sea Ice Extent [$10^6$ km$^2$]', 'Fontsize', 14,'interpreter','latex')

    yyaxis right
    set(gca, 'yticklabel', {})

    ax = gca;
    ax.YAxis(1).Color = 'black'; % Left y-axis
    ax.YAxis(2).Color = 'black'; % Right y-axis

    title('Total Minimum', 'interpreter', 'latex')
    %legend('minimum','Minimum','fontsize', 12,'interpreter','latex')
    %ylim([14.5 17.5])

    xlabel('Year', 'Fontsize', 14, 'interpreter', 'latex')
    set(gca, 'FontSize', 20)
    set(gcf,'color','white')
    xlim([1979 2024])
    grid on






    plot(8)
    subplot(2,3,1)
    plot(years_max,max_ross,'linewidth',4,'color','red', 'LineStyle', '-', 'marker','none')
    hold on
    plot(years_max,y_values_fit_maxross,'linewidth',4,'color','red', 'LineStyle', '--', 'marker','none' )
    title('Ross Maximum', 'interpreter', 'latex')
    ylim([0 7])
    set(gca, 'xticklabel', {})
    ylabel('Regional Sea Ice Extent [$10^6$ km$^2$]', 'Fontsize', 14,'interpreter','latex')
    set(gca, 'FontSize', 20)
    set(gcf,'color','white')
    xlim([1979 2023])
    grid on


    subplot(2,3,2)
    plot(years_max,max_bell,'linewidth',4,'color','#7E2F8E', 'LineStyle', '-', 'marker','none')
    hold on
    plot(years_max,y_values_fit_maxbell,'linewidth',4,'color','#7E2F8E', 'LineStyle', '--', 'marker','none' )
    title('Amundsen-Bellingshausen Maximum', 'interpreter', 'latex')
    %legend('Maximum','Minimum', 'fontsize', 12,'interpreter','latex')
    ylim([0 7])
    set(gca, 'xticklabel', {})
    set(gca, 'yticklabel', {})
    %ylabel('Regional Sea Ice Extent ', 'Fontsize', 14,'interpreter','latex')
    set(gca, 'FontSize', 20)
    set(gcf,'color','white')
   xlim([1979 2023])
    grid on

    subplot(2,3,3)
    plot(years_max,max_wed,'linewidth',4,'color','magenta', 'LineStyle', '-', 'marker','none')
    hold on
    plot(years_max,y_values_fit_maxwed,'linewidth',4,'color','magenta', 'LineStyle', '--', 'marker','none' )
    title('Weddell Maximum', 'interpreter', 'latex')
    %legend('Maximum','Minimum', 'fontsize', 12,'interpreter','latex')
    set(gca, 'xticklabel', {})
    set(gca, 'yticklabel', {})
    ylim([0 7])
    set(gca, 'FontSize', 20)
    set(gcf,'color','white')
   xlim([1979 2023])
    grid on


    subplot(2,3,4)
    plot(years_max,max_kh,'linewidth',4,'color','#77AC30', 'LineStyle', '-', 'marker','none')
    hold on
    plot(years_max,y_values_fit_maxkh,'linewidth',4,'color','#77AC30', 'LineStyle', '--', 'marker','none' )
    title('King Hakon Maximum', 'interpreter', 'latex')
    %legend('Maximum','Minimum', 'fontsize', 12,'interpreter','latex')
    ylim([0 7])
    xlabel('Year', 'Fontsize', 14, 'interpreter', 'latex')
    ylabel('Regional Sea Ice Extent [$10^6$ km$^2$] ', 'Fontsize', 14,'interpreter','latex')
    %ylabel('Regional Sea Ice Extent ', 'Fontsize', 14,'interpreter','latex')
    set(gca, 'FontSize', 20)
    set(gcf,'color','white')
    xlim([1979 2023])
    grid on

    subplot(2,3,5)
    plot(years_max,max_ea,'linewidth',4,'color', 'blue', 'LineStyle', '-', 'marker','none')
    hold on
    plot(years_max,y_values_fit_maxea,'linewidth',4,'color','blue', 'LineStyle', '--', 'marker','none' )
    title('East Antarctica Maximum', 'interpreter', 'latex')
    %legend('Maximum','Minimum', 'fontsize', 12,'interpreter','latex')
    ylim([0 7])
    xlabel('Year', 'Fontsize', 14, 'interpreter', 'latex')
    %set(gca, 'xticklabel', {})
    set(gca, 'yticklabel', {})
    %ylabel('Sea Ice Extent ', 'Fontsize', 14,'interpreter','latex')
    set(gca, 'FontSize', 20)
    set(gcf,'color','white')
    xlim([1979 2023])
    grid on



    subplot(2,3,6)
    yyaxis left
    % set(gca, 'yticklabel', {})

    %yyaxis right
    plot(years_max,max_tot,'linewidth',4,'color','black', 'LineStyle', '-', 'marker','none')
    hold on
    plot(years_max,y_values_fit_maxtot,'linewidth',4,'color','black', 'LineStyle', '--', 'marker','none' )
    ylabel('Total Sea Ice Extent [$10^6$ km$^2$]', 'Fontsize', 14,'interpreter','latex')
    yyaxis right
    set(gca, 'yticklabel', {})

    ax = gca;
    ax.YAxis(1).Color = 'black'; % Left y-axis
    ax.YAxis(2).Color = 'black'; % Right y-axis

    title('Total Maximum', 'interpreter', 'latex')
    %legend('Maximum','Minimum','fontsize', 12,'interpreter','latex')
    %ylim([14.5 17.5])

    xlabel('Year', 'Fontsize', 14, 'interpreter', 'latex')
    set(gca, 'FontSize', 20)
    set(gcf,'color','white')
    xlim([1979 2023])
    grid on





    %% Write data to csv

% Combine year and time of advance into a matrix

max_data = [years_max', max_ross', max_bell',  max_wed',  max_kh', max_ea', max_tot'];
min_data = [years_min', min_ross', min_bell', min_wed', min_kh', min_ea', min_tot'];



% Write the data to a CSV file
csvwrite('annualmax.csv', max_data);
csvwrite('annualmin.csv', min_data);


% Create a README file
readme_text = [
    "This dataset contains the annual maximum (greatest magnitude of each year) " + ...
    "and minimum (smallest magnitude of each year) of Antarctic sea ice extent [10^6 km^2]. From 1979-2024 for minimum. " + ...
    "From 1979-2023 for maximum."
    ""
    "Files included:"
    " 1. annualmax.csv: Contains the following columns:"
    "    Column 1 - Years"
    "    Column 2 - Ross Maximum"
    "    Column 3 - Amundsen-Bellingshausen Maximum"
    "    Column 4 - Weddell Maximum"
    "    Column 5 - King Hakon Maximum"
    "    Column 6 - East Antarctica Maximum"
    "    Column 7 - Total (Continental) Maximum"
    ""
    " 2. annualmin.csv: Contains the following columns:"
    "    Column 1 - Years"
    "    Column 2 - Ross Day Minimum"
    "    Column 3 - Amundsen-Bellingshausen Minimum"
    "    Column 4 - Weddell Day Minimum"
    "    Column 5 - King Hakon Minimum"
    "    Column 6 - East Antarctica Minimum"
    "    Column 7 - Total (Continental) Minimum"
];

% Write the README file
fid = fopen('README.txt', 'wt');
fprintf(fid, '%s\n', readme_text{:});
fclose(fid);

% Optionally, package the files into a ZIP file
zip('ASIE_annualmax_annualmin_csv.zip', {'annualmax.csv', 'annualmin.csv', 'README.txt'});










