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


%% Finding index of maximum and minimum 
% Ross 
% Find the index of the maximum value in each column
[~, maxIndex_ross] = max(rossreshaped);

% Find the index of the minimum value in each column
[~, minIndex_ross] = min(rossreshaped);

% King Hakon 
% Find the index of the maximum value in each column
[~, maxIndex_kh] = max(khreshaped);

% Find the index of the minimum value in each column
[~, minIndex_kh] = min(khreshaped);

% East Antarctica 
% Find the index of the maximum value in each column
[~, maxIndex_ea] = max(eareshaped);

% Find the index of the minimum value in each column
[~, minIndex_ea] = min(eareshaped);

% Weddell
% Find the index of the maximum value in each column
[~, maxIndex_wed] = max(wedreshaped);

% Find the index of the minimum value in each column
[~, minIndex_wed] = min(wedreshaped);

% Bellingshausen 
% Find the index of the maximum value in each column
[~, maxIndex_bell] = max(bellreshaped);

% Find the index of the minimum value in each column
[~, minIndex_bell] = min(bellreshaped);

% Total
% Find the index of the maximum value in each column
[~, maxIndex_tot] = max(totreshaped);

% Find the index of the minimum value in each column
[~, minIndex_tot] = min(totreshaped);

%% Calculate advance amplitudes
for j = 2:size(rossreshaped,2)-1
    advamp_ross(j-1) = rossreshaped(maxIndex_ross(j),j) - rossreshaped(minIndex_ross(j),j);
    advamp_kh(j-1) = khreshaped(maxIndex_kh(j),j) - khreshaped(minIndex_kh(j),j);
    advamp_ea(j-1) = eareshaped(maxIndex_ea(j),j) -eareshaped(minIndex_ea(j),j);
    advamp_wed(j-1) = wedreshaped(maxIndex_wed(j),j) - wedreshaped(minIndex_wed(j),j);
    advamp_bell(j-1) = bellreshaped(maxIndex_bell(j),j) - bellreshaped(minIndex_bell(j),j);
    advamp_tot(j-1) = totreshaped(maxIndex_tot(j),j) - totreshaped(minIndex_tot(j),j);
end


for j = 2:size(rossreshaped,2)-1
    retamp_ross(j-1) = rossreshaped(maxIndex_ross(j),j) - rossreshaped(minIndex_ross(j+1),j);
    retamp_kh(j-1) = khreshaped(maxIndex_kh(j),j) - khreshaped(minIndex_kh(j+1),j);
    retamp_ea(j-1) = eareshaped(maxIndex_ea(j),j) -eareshaped(minIndex_ea(j+1),j);
    retamp_wed(j-1) = wedreshaped(maxIndex_wed(j),j) - wedreshaped(minIndex_wed(j+1),j);
    retamp_bell(j-1) = bellreshaped(maxIndex_bell(j),j) - bellreshaped(minIndex_bell(j+1),j);
    retamp_tot(j-1) = totreshaped(maxIndex_tot(j),j) - totreshaped(minIndex_tot(j+1),j);
end


years_adv = [1979:2023]; 
years_ret = [1979:2023];

% Advance 
p_advross = polyfit(1:length(years_adv), advamp_ross, 1);
a_post = p_advross(1); % Slope (trend)
b_post = p_advross(2); % Intercept
[b_post, bint_post, r_post, rint_post, stats_post] = regress(advamp_ross', [ones(length(years_adv),1), (1:length(years_adv))']);
p_value_advross = stats_post(3);
SSR_post = sum(r_post.^2);
SST_post = sum((advamp_ross - mean(advamp_ross)).^2);
R_squared_advross = 1 - (SSR_post / SST_post);



a_post_advross = b_post(2); % Slope (trend)
slope_confidence_interval_advross = bint_post(2,:); % Confidence interval for the slope
% Calculate the margin of error for the slope
slope_margin_of_error_advross = (slope_confidence_interval_advross(2) - slope_confidence_interval_advross(1)) / 2;





p_advkh = polyfit(1:length(years_adv), advamp_kh, 1);
a_post = p_advkh(1); % Slope (trend)
b_post = p_advkh(2); % Intercept
[b_post, bint_post, r_post, rint_post, stats_post] = regress(advamp_kh', [ones(length(years_adv),1), (1:length(years_adv))']);
p_value_advkh = stats_post(3);
SSR_post = sum(r_post.^2);
SST_post = sum((advamp_kh - mean(advamp_kh)).^2);
R_squared_advkh = 1 - (SSR_post / SST_post);



a_post_advkh = b_post(2); % Slope (trend)
slope_confidence_interval_advkh = bint_post(2,:); % Confidence interval for the slope
% Calculate the margin of error for the slope
slope_margin_of_error_advkh = (slope_confidence_interval_advkh(2) - slope_confidence_interval_advkh(1)) / 2;





p_advea = polyfit(1:length(years_adv), advamp_ea, 1);
a_post = p_advea(1); % Slope (trend)
b_post = p_advea(2); % Intercept
[b_post, bint_post, r_post, rint_post, stats_post] = regress(advamp_ea', [ones(length(years_adv),1), (1:length(years_adv))']);
p_value_advea = stats_post(3);
SSR_post = sum(r_post.^2);
SST_post = sum((advamp_ea - mean(advamp_ea)).^2);
R_squared_advea = 1 - (SSR_post / SST_post);



a_post_advea = b_post(2); % Slope (trend)
slope_confidence_interval_advea = bint_post(2,:); % Confidence interval for the slope
% Calculate the margin of error for the slope
slope_margin_of_error_advea = (slope_confidence_interval_advea(2) - slope_confidence_interval_advea(1)) / 2;





p_advwed = polyfit(1:length(years_adv), advamp_wed, 1);
a_post = p_advwed(1); % Slope (trend)
b_post = p_advwed(2); % Intercept
[b_post, bint_post, r_post, rint_post, stats_post] = regress(advamp_wed', [ones(length(years_adv),1), (1:length(years_adv))']);
p_value_advwed = stats_post(3);
SSR_post = sum(r_post.^2);
SST_post = sum((advamp_wed - mean(advamp_wed)).^2);
R_squared_advwed = 1 - (SSR_post / SST_post);



a_post_advwed = b_post(2); % Slope (trend)
slope_confidence_interval_advwed = bint_post(2,:); % Confidence interval for the slope
% Calculate the margin of error for the slope
slope_margin_of_error_advwed = (slope_confidence_interval_advwed(2) - slope_confidence_interval_advwed(1)) / 2;





p_advbell = polyfit(1:length(years_adv), advamp_bell, 1);
a_post = p_advbell(1); % Slope (trend)
b_post = p_advbell(2); % Intercept
[b_post, bint_post, r_post, rint_post, stats_post] = regress(advamp_bell', [ones(length(years_adv),1), (1:length(years_adv))']);
p_value_advbell = stats_post(3);
SSR_post = sum(r_post.^2);
SST_post = sum((advamp_bell - mean(advamp_bell)).^2);
R_squared_advbell = 1 - (SSR_post / SST_post);



a_post_advbell = b_post(2); % Slope (trend)
slope_confidence_interval_advbell = bint_post(2,:); % Confidence interval for the slope
% Calculate the margin of error for the slope
slope_margin_of_error_advbell = (slope_confidence_interval_advbell(2) - slope_confidence_interval_advbell(1)) / 2;





p_advtot = polyfit(1:length(years_adv), advamp_tot, 1);
a_post = p_advtot(1); % Slope (trend)
b_post = p_advtot(2); % Intercept
[b_post, bint_post, r_post, rint_post, stats_post] = regress(advamp_tot', [ones(length(years_adv),1), (1:length(years_adv))']);
p_value_advtot = stats_post(3);
SSR_post = sum(r_post.^2);
SST_post = sum((advamp_tot - mean(advamp_tot)).^2);
R_squared_advtot = 1 - (SSR_post / SST_post);



a_post_advtot = b_post(2); % Slope (trend)
slope_confidence_interval_advtot = bint_post(2,:); % Confidence interval for the slope
% Calculate the margin of error for the slope
slope_margin_of_error_advtot = (slope_confidence_interval_advtot(2) - slope_confidence_interval_advtot(1)) / 2;




% Retreat 
p_retross = polyfit(1:length(years_ret), retamp_ross, 1);
a_post = p_retross(1); % Slope (trend)
b_post = p_retross(2); % Intercept
[b_post, bint_post, r_post, rint_post, stats_post] = regress(retamp_ross', [ones(length(years_ret),1), (1:length(years_ret))']);
p_value_retross = stats_post(3);
SSR_post = sum(r_post.^2);
SST_post = sum((retamp_ross - mean(retamp_ross)).^2);
R_squared_retross = 1 - (SSR_post / SST_post);


a_post_retross = b_post(2); % Slope (trend)
slope_confidence_interval_retross = bint_post(2,:); % Confidence interval for the slope
% Calculate the margin of error for the slope
slope_margin_of_error_retross = (slope_confidence_interval_retross(2) - slope_confidence_interval_retross(1)) / 2;





p_retkh = polyfit(1:length(years_ret), retamp_kh, 1);
a_post = p_retkh(1); % Slope (trend)
b_post = p_retkh(2); % Intercept
[b_post, bint_post, r_post, rint_post, stats_post] = regress(retamp_kh', [ones(length(years_ret),1), (1:length(years_ret))']);
p_value_retkh = stats_post(3);
SSR_post = sum(r_post.^2);
SST_post = sum((retamp_kh - mean(retamp_kh)).^2);
R_squared_retkh = 1 - (SSR_post / SST_post);


a_post_retkh = b_post(2); % Slope (trend)
slope_confidence_interval_retkh = bint_post(2,:); % Confidence interval for the slope
% Calculate the margin of error for the slope
slope_margin_of_error_retkh = (slope_confidence_interval_retkh(2) - slope_confidence_interval_retkh(1)) / 2;





p_retea = polyfit(1:length(years_ret), retamp_ea, 1);
a_post = p_retea(1); % Slope (trend)
b_post = p_retea(2); % Intercept
[b_post, bint_post, r_post, rint_post, stats_post] = regress(retamp_ea', [ones(length(years_ret),1), (1:length(years_ret))']);
p_value_retea = stats_post(3);
SSR_post = sum(r_post.^2);
SST_post = sum((retamp_ea - mean(retamp_ea)).^2);
R_squared_retea = 1 - (SSR_post / SST_post);


a_post_retea = b_post(2); % Slope (trend)
slope_confidence_interval_retea = bint_post(2,:); % Confidence interval for the slope
% Calculate the margin of error for the slope
slope_margin_of_error_retea = (slope_confidence_interval_retea(2) - slope_confidence_interval_retea(1)) / 2;





p_retwed = polyfit(1:length(years_ret), retamp_wed, 1);
a_post = p_retwed(1); % Slope (trend)
b_post = p_retwed(2); % Intercept
[b_post, bint_post, r_post, rint_post, stats_post] = regress(retamp_wed', [ones(length(years_ret),1), (1:length(years_ret))']);
p_value_retwed = stats_post(3);
SSR_post = sum(r_post.^2);
SST_post = sum((retamp_wed - mean(retamp_wed)).^2);
R_squared_retwed = 1 - (SSR_post / SST_post);


a_post_retwed = b_post(2); % Slope (trend)
slope_confidence_interval_retwed = bint_post(2,:); % Confidence interval for the slope
% Calculate the margin of error for the slope
slope_margin_of_error_retwed = (slope_confidence_interval_retwed(2) - slope_confidence_interval_retwed(1)) / 2;





p_retbell = polyfit(1:length(years_ret), retamp_bell, 1);
a_post = p_retbell(1); % Slope (trend)
b_post = p_retbell(2); % Intercept
[b_post, bint_post, r_post, rint_post, stats_post] = regress(retamp_bell', [ones(length(years_ret),1), (1:length(years_ret))']);
p_value_retbell = stats_post(3);
SSR_post = sum(r_post.^2);
SST_post = sum((retamp_bell - mean(retamp_bell)).^2);
R_squared_retbell = 1 - (SSR_post / SST_post);


a_post_retbell = b_post(2); % Slope (trend)
slope_confidence_interval_retbell = bint_post(2,:); % Confidence interval for the slope
% Calculate the margin of error for the slope
slope_margin_of_error_retbell = (slope_confidence_interval_retbell(2) - slope_confidence_interval_retbell(1)) / 2;





p_rettot = polyfit(1:length(years_ret), retamp_tot, 1);
a_post = p_rettot(1); % Slope (trend)
b_post = p_rettot(2); % Intercept
[b_post, bint_post, r_post, rint_post, stats_post] = regress(retamp_tot', [ones(length(years_ret),1), (1:length(years_ret))']);
p_value_rettot = stats_post(3);
SSR_post = sum(r_post.^2);
SST_post = sum((retamp_tot - mean(retamp_tot)).^2);
R_squared_rettot = 1 - (SSR_post / SST_post);


a_post_rettot = b_post(2); % Slope (trend)
slope_confidence_interval_rettot = bint_post(2,:); % Confidence interval for the slope
% Calculate the margin of error for the slope
slope_margin_of_error_rettot = (slope_confidence_interval_rettot(2) - slope_confidence_interval_rettot(1)) / 2;






x_values_adv = 1:length(years_adv);
y_values_fit_advross = polyval(p_advross, x_values_adv);
y_values_fit_advkh = polyval(p_advkh, x_values_adv);
y_values_fit_advea = polyval(p_advea, x_values_adv);
y_values_fit_advwed = polyval(p_advwed, x_values_adv);
y_values_fit_advbell = polyval(p_advbell, x_values_adv);
y_values_fit_advtot = polyval(p_advtot, x_values_adv);



x_values_ret = 1:length(years_ret);
y_values_fit_retross = polyval(p_retross, x_values_ret);
y_values_fit_retkh = polyval(p_retkh, x_values_ret);
y_values_fit_retea = polyval(p_retea, x_values_ret);
y_values_fit_retwed = polyval(p_retwed, x_values_ret);
y_values_fit_retbell = polyval(p_retbell, x_values_ret);
y_values_fit_rettot = polyval(p_rettot, x_values_ret);

%% Correlate across regions 
% Advance amplitude correlations 
[r_rosskh_advamp, p_rosskh_advamp] = corrcoef(advamp_ross,advamp_kh);
[r_rossea_advamp, p_rossea_advamp] = corrcoef(advamp_ross,advamp_ea);
[r_rosswed_advamp, p_rosswed_advamp] = corrcoef(advamp_ross,advamp_wed);
[r_rossbell_advamp, p_rossbell_advamp] = corrcoef(advamp_ross,advamp_bell);

[r_khea_advamp, p_khea_advamp] = corrcoef(advamp_kh,advamp_ea);
[r_khwed_advamp, p_khwed_advamp] = corrcoef(advamp_kh,advamp_wed);
[r_khbell_advamp, p_khbell_advamp] = corrcoef(advamp_kh,advamp_bell);

[r_eawed_advamp, p_eawed_advamp] = corrcoef(advamp_ea,advamp_wed);
[r_eabell_advamp, p_eabell_advamp] = corrcoef(advamp_ea,advamp_bell);

[r_wedbell_advamp, p_wedbell_advamp] = corrcoef(advamp_wed,advamp_bell);


[r_totross_advamp, p_totross_advamp] = corrcoef(advamp_tot,advamp_ross);
[r_totkh_advamp, p_totkh_advamp] = corrcoef(advamp_tot,advamp_kh);
[r_totea_advamp, p_totea_advamp] = corrcoef(advamp_tot,advamp_ea);
[r_totwed_advamp, p_totwed_advamp] = corrcoef(advamp_tot,advamp_wed);
[r_totbell_advamp, p_totbell_advamp] = corrcoef(advamp_tot,advamp_bell);


% Retreat amplitude correlations
[r_rosskh_retamp, p_rosskh_retamp] = corrcoef(retamp_ross,retamp_kh);
[r_rossea_retamp, p_rossea_retamp] = corrcoef(retamp_ross,retamp_ea);
[r_rosswed_retamp, p_rosswed_retamp] = corrcoef(retamp_ross,retamp_wed);
[r_rossbell_retamp, p_rossbell_retamp] = corrcoef(retamp_ross,retamp_bell);

[r_khea_retamp, p_khea_retamp] = corrcoef(retamp_kh,retamp_ea);
[r_khwed_retamp, p_khwed_retamp] = corrcoef(retamp_kh,retamp_wed);
[r_khbell_retamp, p_khbell_retamp] = corrcoef(retamp_kh,retamp_bell);

[r_eawed_retamp, p_eawed_retamp] = corrcoef(retamp_ea,retamp_wed);
[r_eabell_retamp, p_eabell_retamp] = corrcoef(retamp_ea,retamp_bell);

[r_wedbell_retamp, p_wedbell_retamp] = corrcoef(retamp_wed,retamp_bell);


[r_totross_retamp, p_totross_retamp] = corrcoef(retamp_tot,retamp_ross);
[r_totkh_retamp, p_totkh_retamp] = corrcoef(retamp_tot,retamp_kh);
[r_totea_retamp, p_totea_retamp] = corrcoef(retamp_tot,retamp_ea);
[r_totwed_retamp, p_totwed_retamp] = corrcoef(retamp_tot,retamp_wed);
[r_totbell_retamp, p_totbell_retamp] = corrcoef(retamp_tot,retamp_bell);


%% Plot


plot(9) % Subplot of advance amplitude for [1] five regions and [2] entire continent 
subplot(2,1,1)
plot(years_adv,advamp_ross,'linewidth',4,'color','red', 'LineStyle', '-', 'marker','none')
hold on 
plot(years_adv,advamp_bell,'linewidth',4,'color','#7E2F8E', 'LineStyle', '-', 'marker','none')
hold on 
plot(years_adv,advamp_wed,'linewidth',4,'color','magenta', 'LineStyle', '-', 'marker','none')
hold on 
plot(years_adv,advamp_kh,'linewidth',4,'color','#77AC30', 'LineStyle', '-', 'marker','none')
hold on 
plot(years_adv,advamp_ea,'linewidth',4,'color', 'blue', 'LineStyle', '-', 'marker','none')
ylabel('Sea Ice Extent [$10^6$ km$^2$]', 'Fontsize', 14, 'interpreter','latex')

%xlabel('Year', 'Fontsize', 14, 'interpreter', 'latex')
set(gca, 'FontSize', 20)
set(gca,'XTickLabel',[])
set(gcf,'color','white')
xlim([1979 2023])
ylim([0 7])
grid on

title('Advance Amplitude', 'interpreter', 'latex')
    legend('Ross','Amundsen-Bellingshausen','Weddell','King Hakon','East Antarctica', 'fontsize',12,'interpreter','latex')



subplot(2,1,2)
plot(years_adv,advamp_tot,'linewidth',4,'color','black', 'LineStyle', '-', 'marker','none')

ylabel('Sea Ice Extent [$10^6$ km$^2$]', 'Fontsize', 14, 'interpreter','latex')
xlabel('Year', 'Fontsize', 14, 'interpreter', 'latex')
set(gca, 'FontSize', 20)
%set(gca,'YTickLabel',[])
set(gcf,'color','white')
xlim([1979 2023])
ylim([14.5 17.5])
grid on

%title('Advance Amplitude', 'interpreter', 'latex')
    legend('Total', 'fontsize',12,'interpreter','latex')





 
plot(10) % Subplot of retreat amplitude for [1] five regions and [2] entire continent 
subplot(2,1,1)
plot(years_ret,retamp_ross,'linewidth',4,'color','red', 'LineStyle', '-', 'marker','none')
hold on 
plot(years_ret,retamp_bell,'linewidth',4,'color','#7E2F8E', 'LineStyle', '-', 'marker','none')
hold on 
plot(years_ret,retamp_wed,'linewidth',4,'color','magenta', 'LineStyle', '-', 'marker','none')
hold on 
plot(years_ret,retamp_kh,'linewidth',4,'color','#77AC30', 'LineStyle', '-', 'marker','none')
hold on 
plot(years_ret,retamp_ea,'linewidth',4,'color', 'blue', 'LineStyle', '-', 'marker','none')
ylabel('Sea Ice Extent [$10^6$ km$^2$]', 'Fontsize', 14, 'interpreter','latex')

%xlabel('Year', 'Fontsize', 14, 'interpreter', 'latex')
set(gca, 'FontSize', 20)
set(gca,'XTickLabel',[])
set(gcf,'color','white')
ylim([0 7])
xlim([1979 2023])
grid on

title('Retreat Amplitude', 'interpreter', 'latex')
    legend('Ross','Amundsen-Bellingshausen','Weddell','King Hakon','East Antarctica', 'fontsize',12,'interpreter','latex')



subplot(2,1,2)
plot(years_ret,retamp_tot,'linewidth',4,'color','black', 'LineStyle', '-', 'marker','none')

ylabel('Sea Ice Extent [$10^6$ km$^2$]', 'Fontsize', 14, 'interpreter','latex')
xlabel('Year', 'Fontsize', 14, 'interpreter', 'latex')
set(gca, 'FontSize', 20)
%set(gca,'YTickLabel',[])
set(gcf,'color','white')
ylim([14.5 17.5])
xlim([1979 2023])
grid on

%title('Retreat Amplitude', 'interpreter', 'latex')
    legend('Total','fontsize',12,'interpreter','latex')






plot(11) % Subplot of advance amplitude for [1] five regions and [2] entire continent 
subplot(2,1,1)
plot(years_adv,advamp_ross,'linewidth',4,'color','red', 'LineStyle', '-', 'marker','none')
hold on 
plot(years_adv,advamp_bell,'linewidth',4,'color','#7E2F8E', 'LineStyle', '-', 'marker','none')
hold on 
plot(years_adv,advamp_wed,'linewidth',4,'color','magenta', 'LineStyle', '-', 'marker','none')
hold on 
plot(years_adv,advamp_kh,'linewidth',4,'color','#77AC30', 'LineStyle', '-', 'marker','none')
hold on 
plot(years_adv,advamp_ea,'linewidth',4,'color', 'blue', 'LineStyle', '-', 'marker','none')

%{
hold on 
plot(year,y_values_fit_advross,'linewidth',4,'color','red', 'LineStyle', '--', 'marker','none' )
hold on
plot(year,y_values_fit_advkh,'linewidth',4,'color','#77AC30', 'LineStyle', '--', 'marker','none' )
hold on 
plot(year,y_values_fit_advea,'linewidth',4,'color','blue', 'LineStyle', '--', 'marker','none' )
hold on
plot(year,y_values_fit_advwed,'linewidth',4,'color','magenta', 'LineStyle', '--', 'marker','none' )
hold on
plot(year,y_values_fit_advbell,'linewidth',4,'color','#7E2F8E', 'LineStyle', '--', 'marker','none' )
%}

ylabel('Sea Ice Extent [$10^6$ km$^2$]', 'Fontsize', 14, 'interpreter','latex')
%xlabel('Year', 'Fontsize', 14, 'interpreter', 'latex')
set(gca, 'FontSize', 20)
set(gca,'XTickLabel',[])
set(gcf,'color','white')
xlim([1979 2023])
ylim([0 7])
grid on

title('Advance Amplitude', 'interpreter', 'latex')
    legend('Ross','Amundsen-Bellingshausen','Weddell','King Hakon','East Antarctica', 'fontsize',12,'interpreter','latex')



subplot(2,1,2)
plot(years_adv,advamp_tot,'linewidth',4,'color','black', 'LineStyle', '-', 'marker','none')
%{
hold on
plot(year,y_values_fit_advtot,'linewidth',4,'color','black', 'LineStyle', '--', 'marker','none' )
%}

ylabel('Sea Ice Extent [$10^6$ km$^2$]', 'Fontsize', 14, 'interpreter','latex')
xlabel('Year', 'Fontsize', 14, 'interpreter', 'latex')
set(gca, 'FontSize', 20)
%set(gca,'YTickLabel',[])
set(gcf,'color','white')
xlim([1979 2023])
ylim([14.5 17.5])
grid on

%title('Advance Amplitude', 'interpreter', 'latex')
    legend('Total', 'fontsize',12,'interpreter','latex')



figure(13) 
subplot(2,3,1)
plot(years_adv,advamp_ross,'linewidth',4,'color','red', 'LineStyle', '-', 'marker','none')
hold on 
plot(years_ret,retamp_ross,'linewidth',4,'color',[1, 0.6, 0.6], 'LineStyle', '--', 'marker','none')

title('Ross Amplitudes', 'interpreter', 'latex')
legend('Advance','Retreat', 'fontsize', 12,'interpreter','latex')
ylim([0 7])
set(gca, 'xticklabel', {})
ylabel('Regional Sea Ice Extent [$10^6$ km$^2$]', 'Fontsize', 14,'interpreter','latex')
set(gca, 'FontSize', 20)
set(gcf,'color','white')
xlim([1979 2023])
grid on


subplot(2,3,2)
plot(years_adv,advamp_bell,'linewidth',4,'color','#7E2F8E', 'LineStyle', '-', 'marker','none')
hold on 
plot(years_ret,retamp_bell,'linewidth',4,'color',[0.6, 0.3, 0.7], 'LineStyle', '--', 'marker','none')
set(gca, 'xticklabel', {})
set(gca, 'yticklabel', {})
title('Amundsen-Bellingshausen Amplitudes', 'interpreter', 'latex')
legend('Advance','Retreat', 'fontsize', 12,'interpreter','latex')
ylim([0 7])
set(gca, 'FontSize', 20)
set(gcf,'color','white')
xlim([1979 2023])
grid on


subplot(2,3,3)
plot(years_adv,advamp_wed,'linewidth',4,'color','magenta', 'LineStyle', '-', 'marker','none')
hold on 
plot(years_ret,retamp_wed,'linewidth',4,'color',[1, 0.4, 1], 'LineStyle', '--', 'marker','none')
set(gca, 'xticklabel', {})
set(gca, 'yticklabel', {})
title('Weddell Amplitudes', 'interpreter', 'latex')
legend('Advance','Retreat', 'fontsize', 12,'interpreter','latex')
ylim([0 7])
set(gca, 'FontSize', 20)
set(gcf,'color','white')
xlim([1979 2023])
grid on


subplot(2,3,4)
plot(years_adv,advamp_kh,'linewidth',4,'color','#77AC30', 'LineStyle', '-', 'marker','none')
hold on 
plot(years_ret,retamp_kh,'linewidth',4,'color',[0.667, 0.875, 0.388], 'LineStyle', '--', 'marker','none')
title('King Hakon Amplitudes', 'interpreter', 'latex')
legend('Advance','Retreat', 'fontsize', 12,'interpreter','latex')
ylim([0 7])
xlabel('Year', 'Fontsize', 14, 'interpreter', 'latex')
ylabel('Regional Sea Ice Extent [$10^6$ km$^2$] ', 'Fontsize', 14,'interpreter','latex')
set(gca, 'FontSize', 20)
set(gcf,'color','white')
xlim([1979 2023])
grid on



subplot(2,3,5)
plot(years_adv,advamp_ea,'linewidth',4,'color', 'blue', 'LineStyle', '-', 'marker','none')
hold on 
plot(years_ret,retamp_ea,'linewidth',4,'color', [0.4, 0.4, 1], 'LineStyle', '--', 'marker','none')

title('East Antarctica Amplitudes', 'interpreter', 'latex')
legend('Advance','Retreat', 'fontsize', 12,'interpreter','latex')
ylim([0 7])
xlabel('Year', 'Fontsize', 14, 'interpreter', 'latex')
set(gca, 'yticklabel', {})
set(gca, 'FontSize', 20)
set(gcf,'color','white')
xlim([1979 2023])
grid on


subplot(2,3,6)
yyaxis left 

plot(years_adv,advamp_tot,'linewidth',4,'color','black', 'LineStyle', '-', 'marker','none')
hold on 
plot(years_ret,retamp_tot,'linewidth',4,'color',[0.2, 0.2, 0.2], 'LineStyle', '--', 'marker','none')
ylabel('Total Sea Ice Extent [$10^6$ km$^2$]', 'Fontsize', 14,'interpreter','latex')


yyaxis right
set(gca, 'yticklabel', {})

ax = gca;
ax.YAxis(1).Color = 'black'; % Left y-axis
ax.YAxis(2).Color = 'black'; % Right y-axis

title('Total Amplitudes', 'interpreter', 'latex')
legend('Advance','Retreat','fontsize', 12,'interpreter','latex')
ylim([14.5 17.5])

xlabel('Year', 'Fontsize', 14, 'interpreter', 'latex')
set(gca, 'FontSize', 20)
set(gcf,'color','white')
xlim([1979 2023])
grid on









figure(14) 
subplot(2,3,1)
plot(years_adv,advamp_ross,'linewidth',4,'color','red', 'LineStyle', '-', 'marker','none')
hold on 
plot(years_adv,y_values_fit_advross,'linewidth',4,'color','red', 'LineStyle', '--', 'marker','none' )

title('Ross', 'interpreter', 'latex')
%legend('Advance','Retreat', 'fontsize', 12,'interpreter','latex')
ylim([0 7])
set(gca, 'xticklabel', {})
ylabel('Regional Advance Amplitude [$10^6$ km$^2$]', 'Fontsize', 14,'interpreter','latex')
set(gca, 'FontSize', 20)
set(gcf,'color','white')
xlim([1979 2023])
grid on

subplot(2,3,2)
plot(years_adv,advamp_bell,'linewidth',4,'color','#7E2F8E', 'LineStyle', '-', 'marker','none')
hold on 
plot(years_adv,y_values_fit_advbell,'linewidth',4,'color','#7E2F8E', 'LineStyle', '--', 'marker','none' )
set(gca, 'xticklabel', {})
set(gca, 'yticklabel', {})
title('Amundsen-Bellingshausen', 'interpreter', 'latex')
%legend('Advance','Retreat', 'fontsize', 12,'interpreter','latex')
ylim([0 7])
%xlabel('Year', 'Fontsize', 14, 'interpreter', 'latex')
set(gca, 'yticklabel', {})
%ylabel('Regional Sea Ice Extent ', 'Fontsize', 14,'interpreter','latex')
set(gca, 'FontSize', 20)
set(gcf,'color','white')
xlim([1979 2023])
grid on

subplot(2,3,3)
plot(years_adv,advamp_wed,'linewidth',4,'color','magenta', 'LineStyle', '-', 'marker','none')
hold on 
plot(years_adv,y_values_fit_advwed,'linewidth',4,'color','magenta', 'LineStyle', '--', 'marker','none' )
set(gca, 'xticklabel', {})
set(gca, 'yticklabel', {})
title('Weddell', 'interpreter', 'latex')
%legend('Advance','Retreat', 'fontsize', 12,'interpreter','latex')
ylim([0 7])

set(gca, 'FontSize', 20)
set(gcf,'color','white')
xlim([1979 2023])
grid on


subplot(2,3,4)
plot(years_adv,advamp_kh,'linewidth',4,'color','#77AC30', 'LineStyle', '-', 'marker','none')
hold on 
plot(years_adv,y_values_fit_advkh,'linewidth',4,'color','#77AC30', 'LineStyle', '--', 'marker','none' )

title('King Hakon', 'interpreter', 'latex')
%legend('Advance','Retreat', 'fontsize', 12,'interpreter','latex')
ylim([0 7])
xlabel('Year', 'Fontsize', 14, 'interpreter', 'latex')
ylabel('Regional Advance Amplitude [$10^6$ km$^2$]', 'Fontsize', 14,'interpreter','latex')
%ylabel('Regional Sea Ice Extent ', 'Fontsize', 14,'interpreter','latex')
set(gca, 'FontSize', 20)
set(gcf,'color','white')
xlim([1979 2023])
grid on


subplot(2,3,5)
plot(years_adv,advamp_ea,'linewidth',4,'color', 'blue', 'LineStyle', '-', 'marker','none')
hold on 
plot(years_adv,y_values_fit_advea,'linewidth',4,'color','blue', 'LineStyle', '--', 'marker','none' )
title('East Antarctica', 'interpreter', 'latex')
%legend('Advance','Retreat', 'fontsize', 12,'interpreter','latex')
ylim([0 7])
xlabel('Year', 'Fontsize', 14, 'interpreter', 'latex')
set(gca, 'yticklabel', {})
%ylabel('Sea Ice Extent ', 'Fontsize', 14,'interpreter','latex')
set(gca, 'FontSize', 20)
set(gcf,'color','white')
xlim([1979 2023])
grid on


subplot(2,3,6)
yyaxis left 


plot(years_adv,advamp_tot,'linewidth',4,'color','black', 'LineStyle', '-', 'marker','none')
hold on
plot(years_adv,y_values_fit_advtot,'linewidth',4,'color','black', 'LineStyle', '--', 'marker','none' )
ylabel('Total Advance Amplitude [$10^6$ km$^2$]', 'Fontsize', 14,'interpreter','latex')

yyaxis right
set(gca, 'yticklabel', {})

ax = gca;
ax.YAxis(1).Color = 'black'; % Left y-axis
ax.YAxis(2).Color = 'black'; % Right y-axis

title('Total', 'interpreter', 'latex')
%legend('Advance','Retreat','fontsize', 12,'interpreter','latex')
ylim([14.5 17.5])
xlabel('Year', 'Fontsize', 14, 'interpreter', 'latex')
set(gca, 'FontSize', 20)
set(gcf,'color','white')
xlim([1979 2023])
grid on




figure(15) 
subplot(2,3,1)
plot(years_ret,retamp_ross,'linewidth',4,'color','red', 'LineStyle', '-', 'marker','none')
hold on 
plot(years_ret,y_values_fit_retross,'linewidth',4,'color','red', 'LineStyle', '--', 'marker','none' )

title('Ross', 'interpreter', 'latex')
%legend('retance','Retreat', 'fontsize', 12,'interpreter','latex')
ylim([0 7])
set(gca, 'xticklabel', {})
ylabel('Regional Retreat Amplitude [$10^6$ km$^2$]', 'Fontsize', 14,'interpreter','latex')
set(gca, 'FontSize', 20)
set(gcf,'color','white')
xlim([1979 2023])
grid on


subplot(2,3,2)
plot(years_ret,retamp_bell,'linewidth',4,'color','#7E2F8E', 'LineStyle', '-', 'marker','none')
hold on 
plot(years_ret,y_values_fit_retbell,'linewidth',4,'color','#7E2F8E', 'LineStyle', '--', 'marker','none' )

title('Amundsen-Bellingshausen', 'interpreter', 'latex')
%legend('retance','Retreat', 'fontsize', 12,'interpreter','latex')
ylim([0 7])
set(gca, 'xticklabel', {})
set(gca, 'yticklabel', {})
set(gca, 'yticklabel', {})
%ylabel('Regional Sea Ice Extent ', 'Fontsize', 14,'interpreter','latex')
set(gca, 'FontSize', 20)
set(gcf,'color','white')
xlim([1979 2023])
grid on


subplot(2,3,3)
plot(years_ret,retamp_wed,'linewidth',4,'color','magenta', 'LineStyle', '-', 'marker','none')
hold on 
plot(years_ret,y_values_fit_retwed,'linewidth',4,'color','magenta', 'LineStyle', '--', 'marker','none' )

title('Weddell', 'interpreter', 'latex')
%legend('retance','Retreat', 'fontsize', 12,'interpreter','latex')
ylim([0 7])
set(gca, 'xticklabel', {})
set(gca, 'yticklabel', {})
set(gca, 'FontSize', 20)
set(gcf,'color','white')
xlim([1979 2023])
grid on

subplot(2,3,4)
plot(years_ret,retamp_kh,'linewidth',4,'color','#77AC30', 'LineStyle', '-', 'marker','none')
hold on 
plot(years_ret,y_values_fit_retkh,'linewidth',4,'color','#77AC30', 'LineStyle', '--', 'marker','none' )

title('King Hakon', 'interpreter', 'latex')
%legend('retance','Retreat', 'fontsize', 12,'interpreter','latex')
ylim([0 7])
xlabel('Year', 'Fontsize', 14, 'interpreter', 'latex')
ylabel('Regional Retreat Amplitude [$10^6$ km$^2$]', 'Fontsize', 14,'interpreter','latex')
set(gca, 'FontSize', 20)
set(gcf,'color','white')
xlim([1979 2023])
grid on

subplot(2,3,5)
plot(years_ret,retamp_ea,'linewidth',4,'color', 'blue', 'LineStyle', '-', 'marker','none')
hold on 
plot(years_ret,y_values_fit_retea,'linewidth',4,'color','blue', 'LineStyle', '--', 'marker','none' )
xlabel('Year', 'Fontsize', 14, 'interpreter', 'latex')
set(gca, 'yticklabel', {})
title('East Antarctica', 'interpreter', 'latex')
%legend('retance','Retreat', 'fontsize', 12,'interpreter','latex')
ylim([0 7])

%ylabel('Sea Ice Extent ', 'Fontsize', 14,'interpreter','latex')
set(gca, 'FontSize', 20)
set(gcf,'color','white')
xlim([1979 2023])
grid on



subplot(2,3,6)
yyaxis left 


plot(years_ret,retamp_tot,'linewidth',4,'color','black', 'LineStyle', '-', 'marker','none')
hold on
plot(years_ret,y_values_fit_rettot,'linewidth',4,'color','black', 'LineStyle', '--', 'marker','none' )
ylabel('Total Retreat Amplitude [$10^6$ km$^2$]', 'Fontsize', 14,'interpreter','latex')

yyaxis right
set(gca, 'yticklabel', {})

ax = gca;
ax.YAxis(1).Color = 'black'; % Left y-axis
ax.YAxis(2).Color = 'black'; % Right y-axis

title('Total', 'interpreter', 'latex')
%legend('Advance','Retreat','fontsize', 12,'interpreter','latex')
ylim([14.5 17.5])
xlabel('Year', 'Fontsize', 14, 'interpreter', 'latex')
set(gca, 'FontSize', 20)
set(gcf,'color','white')
xlim([1979 2023])
grid on





    %% Write data to csv

% Combine year and time of advance into a matrix

advamplitude = [years_adv', advamp_ross', advamp_bell',  advamp_wed',  advamp_kh', advamp_ea', advamp_tot'];
retamplitude = [years_ret', retamp_ross', retamp_bell', retamp_wed', retamp_kh', retamp_ea', retamp_tot'];



% Write the data to a CSV file
csvwrite('advamplitude.csv', advamplitude);
csvwrite('retamplitude.csv', retamplitude);


% Create a README file
readme_text = [
    "This dataset contains the annual advance amplitude (max_extent(year)-min_extent(year)) " + ...
    "and annual retreat amplitude (max_extent(year)-min_extent(year+1)) of Antarctic sea ice extent [10^6 km^2]. " + ...
    "From 1979-2023 for both advance and retreat amplitudes."
    ""
    "Files included:"
    " 1. advamplitude.csv: Contains the following columns:"
    "    Column 1 - Years"
    "    Column 2 - Ross Advance Amplitude"
    "    Column 3 - Amundsen-Bellingshausen Advance Amplitude"
    "    Column 4 - Weddell Advance Amplitude"
    "    Column 5 - King Hakon Advance Amplitude"
    "    Column 6 - East Antarctica Advance Amplitude"
    "    Column 7 - Total (Continental) Advance Amplitude"
    ""
    " 2. retamplitude.csv: Contains the following columns:"
    "    Column 1 - Years"
    "    Column 2 - Ross Day Retreat Amplitude"
    "    Column 3 - Amundsen-Bellingshausen Retreat Amplitude"
    "    Column 4 - Weddell Day Retreat Amplitude"
    "    Column 5 - King Hakon Retreat Amplitude"
    "    Column 6 - East Antarctica Retreat Amplitude"
    "    Column 7 - Total (Continental) Retreat Amplitude"
];

% Write the README file
fid = fopen('README.txt', 'wt');
fprintf(fid, '%s\n', readme_text{:});
fclose(fid);

% Optionally, package the files into a ZIP file
zip('ASIE_advamplitude_retamplitude_csv.zip', {'advamplitude.csv', 'retamplitude.csv', 'README.txt'});











