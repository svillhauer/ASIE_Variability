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



%% Maximum 
rossadvindex=[]; 
[rossadv, idxrossadv] = min(rossreshaped);
rossadvindex = [rossadvindex idxrossadv];


khadvindex=[]; 
[khadv, idxkhadv] = min(khreshaped);
khadvindex = [khadvindex idxkhadv];


eaadvindex=[]; 
[eaadv, idxeaadv] = min(eareshaped);
eaadvindex = [eaadvindex idxeaadv];


wedadvindex=[]; 
[wedadv, idxwedadv] = min(wedreshaped);
wedadvindex = [wedadvindex idxwedadv];


belladvindex=[]; 
[belladv, idxbelladv] = min(bellreshaped);
belladvindex = [belladvindex idxbelladv];




totadvindex=[]; 
[totadv, idxtotadv] = min(totreshaped);
totadvindex = [totadvindex idxtotadv];

years_adv = [1979:2024];
for j = 2:47
    ross_adv(j-1) = timereshaped(rossadvindex(j),j);
end

ross_adv_decyears = ross_adv; 

for j = 1:length(years_adv)
    ross_adv(j) = 365*(ross_adv(j)-years_adv(j));
end


for j = 2:47
    kh_adv(j-1) = timereshaped(khadvindex(j),j);
end

kh_adv_decyears = kh_adv; 


for j = 1:length(years_adv)
    kh_adv(j) = 365*(kh_adv(j)-years_adv(j));
end

for j = 2:47
    ea_adv(j-1) = timereshaped(eaadvindex(j),j);
end

ea_adv_decyears = ea_adv; 


for j = 1:length(years_adv)
    ea_adv(j) = 365*(ea_adv(j)-years_adv(j));
end

for j = 2:47
    wed_adv(j-1) = timereshaped(wedadvindex(j),j);
end

wed_adv_decyears = wed_adv; 


for j = 1:length(years_adv)
    wed_adv(j) = 365*(wed_adv(j)-years_adv(j));
end

for j = 2:47
    bell_adv(j-1) = timereshaped(belladvindex(j),j);
end

bell_adv_decyears = bell_adv; 


for j = 1:length(years_adv)
    bell_adv(j) = 365*(bell_adv(j)-years_adv(j));
end

for j = 2:47
    tot_adv(j-1) = timereshaped(totadvindex(j),j);
end

tot_adv_decyears = tot_adv; 


for j = 1:length(years_adv)
    tot_adv(j) = 365*(tot_adv(j)-years_adv(j));
end




%% Minimum 
rossretindex=[]; 
[rossret, idxrossret] = max(rossreshaped);
rossretindex = [rossretindex idxrossret];


khretindex=[]; 
[khret, idxkhret] = max(khreshaped);
khretindex = [khretindex idxkhret];


earetindex=[]; 
[earet, idxearet] = max(eareshaped);
earetindex = [earetindex idxearet];


wedretindex=[]; 
[wedret, idxwedret] = max(wedreshaped);
wedretindex = [wedretindex idxwedret];


bellretindex=[]; 
[bellret, idxbellret] = max(bellreshaped);
bellretindex = [bellretindex idxbellret];


totretindex=[]; 
[totret, idxtotret] = max(totreshaped);
totretindex = [totretindex idxtotret];



years_ret = [1979:2023];
for j = 2:46
    ross_ret(j-1) = timereshaped(rossretindex(j),j);
end

ross_ret_decyears = ross_ret; 

for j = 1:length(years_ret)
    ross_ret(j) = 365*(ross_ret(j)-years_ret(j));
end


for j = 2:46
    kh_ret(j-1) = timereshaped(khretindex(j),j);
end

kh_ret_decyears = kh_ret; 

for j = 1:length(years_ret)
    kh_ret(j) = 365*(kh_ret(j)-years_ret(j));
end

for j = 2:46
    ea_ret(j-1) = timereshaped(earetindex(j),j);
end

ea_ret_decyears = ea_ret; 

for j = 1:length(years_ret)
    ea_ret(j) = 365*(ea_ret(j)-years_ret(j));
end

for j = 2:46
    wed_ret(j-1) = timereshaped(wedretindex(j),j);
end

wed_ret_decyears = wed_ret; 

for j = 1:length(years_ret)
    wed_ret(j) = 365*(wed_ret(j)-years_ret(j));
end

for j = 2:46
    bell_ret(j-1) = timereshaped(bellretindex(j),j);
end

bell_ret_decyears = bell_ret; 

for j = 1:length(years_ret)
    bell_ret(j) = 365*(bell_ret(j)-years_ret(j));
end

for j = 2:46
    tot_ret(j-1) = timereshaped(totretindex(j),j);
end

tot_ret_decyears = tot_ret; 

for j = 1:length(years_ret)
    tot_ret(j) = 365*(tot_ret(j)-years_ret(j));
end

%% Calculate linear trends of time of advance and retreat 
% Advance 

% Advance 
p_advross = polyfit(1:length(years_adv), ross_adv, 1);
a_post = p_advross(1); % Slope (trend)
b_post = p_advross(2); % Intercept
[b_post, bint_post, r_post, rint_post, stats_post] = regress(ross_adv', [ones(length(years_adv),1), (1:length(years_adv))']);
p_value_advross = stats_post(3);
SSR_post = sum(r_post.^2);
SST_post = sum((ross_adv - mean(ross_adv)).^2);
R_squared_advross = 1 - (SSR_post / SST_post);



a_post_advross = b_post(2); % Slope (trend)
slope_confidence_interval_advross = bint_post(2,:); % Confidence interval for the slope
% Calculate the margin of error for the slope
slope_margin_of_error_advross = (slope_confidence_interval_advross(2) - slope_confidence_interval_advross(1)) / 2;





p_advkh = polyfit(1:length(years_adv), kh_adv, 1);
a_post = p_advkh(1); % Slope (trend)
b_post = p_advkh(2); % Intercept
[b_post, bint_post, r_post, rint_post, stats_post] = regress(kh_adv', [ones(length(years_adv),1), (1:length(years_adv))']);
p_value_advkh = stats_post(3);
SSR_post = sum(r_post.^2);
SST_post = sum((kh_adv - mean(kh_adv)).^2);
R_squared_advkh = 1 - (SSR_post / SST_post);



a_post_advkh = b_post(2); % Slope (trend)
slope_confidence_interval_advkh = bint_post(2,:); % Confidence interval for the slope
% Calculate the margin of error for the slope
slope_margin_of_error_advkh = (slope_confidence_interval_advkh(2) - slope_confidence_interval_advkh(1)) / 2;





p_advea = polyfit(1:length(years_adv), ea_adv, 1);
a_post = p_advea(1); % Slope (trend)
b_post = p_advea(2); % Intercept
[b_post, bint_post, r_post, rint_post, stats_post] = regress(ea_adv', [ones(length(years_adv),1), (1:length(years_adv))']);
p_value_advea = stats_post(3);
SSR_post = sum(r_post.^2);
SST_post = sum((ea_adv - mean(ea_adv)).^2);
R_squared_advea = 1 - (SSR_post / SST_post);



a_post_advea = b_post(2); % Slope (trend)
slope_confidence_interval_advea = bint_post(2,:); % Confidence interval for the slope
% Calculate the margin of error for the slope
slope_margin_of_error_advea = (slope_confidence_interval_advea(2) - slope_confidence_interval_advea(1)) / 2;





p_advwed = polyfit(1:length(years_adv), wed_adv, 1);
a_post = p_advwed(1); % Slope (trend)
b_post = p_advwed(2); % Intercept
[b_post, bint_post, r_post, rint_post, stats_post] = regress(wed_adv', [ones(length(years_adv),1), (1:length(years_adv))']);
p_value_advwed = stats_post(3);
SSR_post = sum(r_post.^2);
SST_post = sum((wed_adv - mean(wed_adv)).^2);
R_squared_advwed = 1 - (SSR_post / SST_post);



a_post_advwed = b_post(2); % Slope (trend)
slope_confidence_interval_advwed = bint_post(2,:); % Confidence interval for the slope
% Calculate the margin of error for the slope
slope_margin_of_error_advwed = (slope_confidence_interval_advwed(2) - slope_confidence_interval_advwed(1)) / 2;





p_advbell = polyfit(1:length(years_adv), bell_adv, 1);
a_post = p_advbell(1); % Slope (trend)
b_post = p_advbell(2); % Intercept
[b_post, bint_post, r_post, rint_post, stats_post] = regress(bell_adv', [ones(length(years_adv),1), (1:length(years_adv))']);
p_value_advbell = stats_post(3);
SSR_post = sum(r_post.^2);
SST_post = sum((bell_adv - mean(bell_adv)).^2);
R_squared_advbell = 1 - (SSR_post / SST_post);



a_post_advbell = b_post(2); % Slope (trend)
slope_confidence_interval_advbell = bint_post(2,:); % Confidence interval for the slope
% Calculate the margin of error for the slope
slope_margin_of_error_advbell = (slope_confidence_interval_advbell(2) - slope_confidence_interval_advbell(1)) / 2;





p_advtot = polyfit(1:length(years_adv), tot_adv, 1);
a_post = p_advtot(1); % Slope (trend)
b_post = p_advtot(2); % Intercept
[b_post, bint_post, r_post, rint_post, stats_post] = regress(tot_adv', [ones(length(years_adv),1), (1:length(years_adv))']);
p_value_advtot = stats_post(3);
SSR_post = sum(r_post.^2);
SST_post = sum((tot_adv - mean(tot_adv)).^2);
R_squared_advtot = 1 - (SSR_post / SST_post);



a_post_advtot = b_post(2); % Slope (trend)
slope_confidence_interval_advtot = bint_post(2,:); % Confidence interval for the slope
% Calculate the margin of error for the slope
slope_margin_of_error_advtot = (slope_confidence_interval_advtot(2) - slope_confidence_interval_advtot(1)) / 2;




% Retreat 
p_retross = polyfit(1:length(years_ret), ross_ret, 1);
a_post = p_retross(1); % Slope (trend)
b_post = p_retross(2); % Intercept
[b_post, bint_post, r_post, rint_post, stats_post] = regress(ross_ret', [ones(length(years_ret),1), (1:length(years_ret))']);
p_value_retross = stats_post(3);
SSR_post = sum(r_post.^2);
SST_post = sum((ross_ret - mean(ross_ret)).^2);
R_squared_retross = 1 - (SSR_post / SST_post);



a_post_retross = b_post(2); % Slope (trend)
slope_confidence_interval_retross = bint_post(2,:); % Confidence interval for the slope
% Calculate the margin of error for the slope
slope_margin_of_error_retross = (slope_confidence_interval_retross(2) - slope_confidence_interval_retross(1)) / 2;





p_retkh = polyfit(1:length(years_ret), kh_ret, 1);
a_post = p_retkh(1); % Slope (trend)
b_post = p_retkh(2); % Intercept
[b_post, bint_post, r_post, rint_post, stats_post] = regress(kh_ret', [ones(length(years_ret),1), (1:length(years_ret))']);
p_value_retkh = stats_post(3);
SSR_post = sum(r_post.^2);
SST_post = sum((kh_ret - mean(kh_ret)).^2);
R_squared_retkh = 1 - (SSR_post / SST_post);



a_post_retkh = b_post(2); % Slope (trend)
slope_confidence_interval_retkh = bint_post(2,:); % Confidence interval for the slope
% Calculate the margin of error for the slope
slope_margin_of_error_retkh = (slope_confidence_interval_retkh(2) - slope_confidence_interval_retkh(1)) / 2;





p_retea = polyfit(1:length(years_ret), ea_ret, 1);
a_post = p_retea(1); % Slope (trend)
b_post = p_retea(2); % Intercept
[b_post, bint_post, r_post, rint_post, stats_post] = regress(ea_ret', [ones(length(years_ret),1), (1:length(years_ret))']);
p_value_retea = stats_post(3);
SSR_post = sum(r_post.^2);
SST_post = sum((ea_ret - mean(ea_ret)).^2);
R_squared_retea = 1 - (SSR_post / SST_post);



a_post_retea = b_post(2); % Slope (trend)
slope_confidence_interval_retea = bint_post(2,:); % Confidence interval for the slope
% Calculate the margin of error for the slope
slope_margin_of_error_retea = (slope_confidence_interval_retea(2) - slope_confidence_interval_retea(1)) / 2;




p_retwed = polyfit(1:length(years_ret), wed_ret, 1);
a_post = p_retwed(1); % Slope (trend)
b_post = p_retwed(2); % Intercept
[b_post, bint_post, r_post, rint_post, stats_post] = regress(wed_ret', [ones(length(years_ret),1), (1:length(years_ret))']);
p_value_retwed = stats_post(3);
SSR_post = sum(r_post.^2);
SST_post = sum((wed_ret- mean(wed_ret)).^2);
R_squared_retwed = 1 - (SSR_post / SST_post);



a_post_retwed = b_post(2); % Slope (trend)
slope_confidence_interval_retwed = bint_post(2,:); % Confidence interval for the slope
% Calculate the margin of error for the slope
slope_margin_of_error_retwed = (slope_confidence_interval_retwed(2) - slope_confidence_interval_retwed(1)) / 2;





p_retbell = polyfit(1:length(years_ret), bell_ret, 1);
a_post = p_retbell(1); % Slope (trend)
b_post = p_retbell(2); % Intercept
[b_post, bint_post, r_post, rint_post, stats_post] = regress(bell_ret', [ones(length(years_ret),1), (1:length(years_ret))']);
p_value_retbell = stats_post(3);
SSR_post = sum(r_post.^2);
SST_post = sum((bell_ret - mean(bell_ret)).^2);
R_squared_retbell = 1 - (SSR_post / SST_post);



a_post_retbell = b_post(2); % Slope (trend)
slope_confidence_interval_retbell = bint_post(2,:); % Confidence interval for the slope
% Calculate the margin of error for the slope
slope_margin_of_error_retbell = (slope_confidence_interval_retbell(2) - slope_confidence_interval_retbell(1)) / 2;





p_rettot = polyfit(1:length(years_ret), tot_ret, 1);
a_post = p_rettot(1); % Slope (trend)
b_post = p_rettot(2); % Intercept
[b_post, bint_post, r_post, rint_post, stats_post] = regress(tot_ret', [ones(length(years_ret),1), (1:length(years_ret))']);
p_value_rettot = stats_post(3);
SSR_post = sum(r_post.^2);
SST_post = sum((tot_ret - mean(tot_ret)).^2);
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
% Time of advance correlations 
[r_rosskh_adv, p_rosskh_adv] = corrcoef(ross_adv,kh_adv);
[r_rossea_adv, p_rossea_adv] = corrcoef(ross_adv,ea_adv);
[r_rosswed_adv, p_rosswed_adv] = corrcoef(ross_adv,wed_adv);
[r_rossbell_adv, p_rossbell_adv] = corrcoef(ross_adv,bell_adv);

[r_khea_adv, p_khea_adv] = corrcoef(kh_adv,ea_adv);
[r_khwed_adv, p_khwed_adv] = corrcoef(kh_adv,wed_adv);
[r_khbell_adv, p_khbell_adv] = corrcoef(kh_adv,bell_adv);

[r_eawed_adv, p_eawed_adv] = corrcoef(ea_adv,wed_adv);
[r_eabell_adv, p_eabell_adv] = corrcoef(ea_adv,bell_adv);

[r_wedbell_adv, p_wedbell_adv] = corrcoef(wed_adv,bell_adv);


[r_totross_adv, p_totross_adv] = corrcoef(tot_adv,ross_adv);
[r_totkh_adv, p_totkh_adv] = corrcoef(tot_adv,kh_adv);
[r_totea_adv, p_totea_adv] = corrcoef(tot_adv,ea_adv);
[r_totwed_adv, p_totwed_adv] = corrcoef(tot_adv,wed_adv);
[r_totbell_adv, p_totbell_adv] = corrcoef(tot_adv,bell_adv);


% Time of retreat correlations 
[r_rosskh_ret, p_rosskh_ret] = corrcoef(ross_ret,kh_ret);
[r_rossea_ret, p_rossea_ret] = corrcoef(ross_ret,ea_ret);
[r_rosswed_ret, p_rosswed_ret] = corrcoef(ross_ret,wed_ret);
[r_rossbell_ret, p_rossbell_ret] = corrcoef(ross_ret,bell_ret);

[r_khea_ret, p_khea_ret] = corrcoef(kh_ret,ea_ret);
[r_khwed_ret, p_khwed_ret] = corrcoef(kh_ret,wed_ret);
[r_khbell_ret, p_khbell_ret] = corrcoef(kh_ret,bell_ret);

[r_eawed_ret, p_eawed_ret] = corrcoef(ea_ret,wed_ret);
[r_eabell_ret, p_eabell_ret] = corrcoef(ea_ret,bell_ret);

[r_wedbell_ret, p_wedbell_ret] = corrcoef(wed_ret,bell_ret);


[r_totross_ret, p_totross_ret] = corrcoef(tot_ret,ross_ret);
[r_totkh_ret, p_totkh_ret] = corrcoef(tot_ret,kh_ret);
[r_totea_ret, p_totea_ret] = corrcoef(tot_ret,ea_ret);
[r_totwed_ret, p_totwed_ret] = corrcoef(tot_ret,wed_ret);
[r_totbell_ret, p_totbell_ret] = corrcoef(tot_ret,bell_ret);






  %% Subplots 
figure(1)  
subplot(2,3,1)
plot(years_adv,ross_adv,'linewidth',4,'color','red', 'LineStyle', '-', 'marker','none')
hold on 
plot(years_adv,y_values_fit_advross,'linewidth',4,'color','red', 'LineStyle', '--', 'marker','none' )
title('Ross', 'interpreter', 'latex')
ylim([30 120])
set(gca, 'xticklabel', {})
ylabel('Julian Day of Advance', 'Fontsize', 14,'interpreter','latex')
set(gca, 'FontSize', 20)
set(gcf,'color','white')
xlim([1979 2023])
grid on

subplot(2,3,2)
plot(years_adv,bell_adv,'linewidth',4,'color','#7E2F8E', 'LineStyle', '-', 'marker','none')
hold on 
plot(years_adv,y_values_fit_advbell,'linewidth',4,'color','#7E2F8E', 'LineStyle', '--', 'marker','none' )
set(gca, 'xticklabel', {})
set(gca, 'yticklabel', {})
title('Amundsen-Bellingshausen', 'interpreter', 'latex')
ylim([30 120])
set(gca, 'yticklabel', {})
set(gca, 'FontSize', 20)
set(gcf,'color','white')
xlim([1979 2023])
grid on

subplot(2,3,3)
plot(years_adv,wed_adv,'linewidth',4,'color','magenta', 'LineStyle', '-', 'marker','none')
hold on 
plot(years_adv,y_values_fit_advwed,'linewidth',4,'color','magenta', 'LineStyle', '--', 'marker','none' )
set(gca, 'xticklabel', {})
set(gca, 'yticklabel', {})
title('Weddell', 'interpreter', 'latex')
ylim([30 120])
set(gca, 'FontSize', 20)
set(gcf,'color','white')
xlim([1979 2023])
grid on


subplot(2,3,4)
plot(years_adv,kh_adv,'linewidth',4,'color','#77AC30', 'LineStyle', '-', 'marker','none')
hold on 
plot(years_adv,y_values_fit_advkh,'linewidth',4,'color','#77AC30', 'LineStyle', '--', 'marker','none' )
title('King Hakon', 'interpreter', 'latex')
ylim([30 120])
xlabel('Year', 'Fontsize', 14, 'interpreter', 'latex')
ylabel('Julian Day of Advance', 'Fontsize', 14,'interpreter','latex')
set(gca, 'FontSize', 20)
set(gcf,'color','white')
xlim([1979 2023])
grid on


subplot(2,3,5)
plot(years_adv,ea_adv,'linewidth',4,'color', 'blue', 'LineStyle', '-', 'marker','none')
hold on 
plot(years_adv,y_values_fit_advea,'linewidth',4,'color','blue', 'LineStyle', '--', 'marker','none' )
title('East Antarctica', 'interpreter', 'latex')
ylim([30 120])
xlabel('Year', 'Fontsize', 14, 'interpreter', 'latex')
set(gca, 'yticklabel', {})
set(gca, 'FontSize', 20)
set(gcf,'color','white')
xlim([1979 2023])
grid on


subplot(2,3,6)
plot(years_adv,tot_adv,'linewidth',4,'color','black', 'LineStyle', '-', 'marker','none')
hold on
plot(years_adv,y_values_fit_advtot,'linewidth',4,'color','black', 'LineStyle', '--', 'marker','none' )
set(gca, 'yticklabel', {})
title('Total', 'interpreter', 'latex')
ylim([30 120])
xlabel('Year', 'Fontsize', 14, 'interpreter', 'latex')
set(gca, 'FontSize', 20)
set(gcf,'color','white')
xlim([1979 2023])
grid on



figure(2)  
subplot(2,3,1)
plot(years_ret,ross_ret,'linewidth',4,'color','red', 'LineStyle', '-', 'marker','none')
hold on 
plot(years_ret,y_values_fit_retross,'linewidth',4,'color','red', 'LineStyle', '--', 'marker','none' )
title('Ross', 'interpreter', 'latex')
ylim([190 310])
set(gca, 'xticklabel', {})
ylabel('Julian Day of Retreat', 'Fontsize', 14,'interpreter','latex')
set(gca, 'FontSize', 20)
set(gcf,'color','white')
xlim([1979 2023])
grid on

subplot(2,3,2)
plot(years_ret,bell_ret,'linewidth',4,'color','#7E2F8E', 'LineStyle', '-', 'marker','none')
hold on 
plot(years_ret,y_values_fit_retbell,'linewidth',4,'color','#7E2F8E', 'LineStyle', '--', 'marker','none' )
set(gca, 'xticklabel', {})
set(gca, 'yticklabel', {})
title('Amundsen-Bellingshausen', 'interpreter', 'latex')
ylim([190 310])
set(gca, 'yticklabel', {})
set(gca, 'FontSize', 20)
set(gcf,'color','white')
xlim([1979 2023])
grid on

subplot(2,3,3)
plot(years_ret,wed_ret,'linewidth',4,'color','magenta', 'LineStyle', '-', 'marker','none')
hold on 
plot(years_ret,y_values_fit_retwed,'linewidth',4,'color','magenta', 'LineStyle', '--', 'marker','none' )
set(gca, 'xticklabel', {})
set(gca, 'yticklabel', {})
title('Weddell', 'interpreter', 'latex')
ylim([190 310])
set(gca, 'FontSize', 20)
set(gcf,'color','white')
xlim([1979 2023])
grid on


subplot(2,3,4)
plot(years_ret,kh_ret,'linewidth',4,'color','#77AC30', 'LineStyle', '-', 'marker','none')
hold on 
plot(years_ret,y_values_fit_retkh,'linewidth',4,'color','#77AC30', 'LineStyle', '--', 'marker','none' )
title('King Hakon', 'interpreter', 'latex')
ylim([190 310])
xlabel('Year', 'Fontsize', 14, 'interpreter', 'latex')
ylabel('Julian Day of Retreat', 'Fontsize', 14,'interpreter','latex')
set(gca, 'FontSize', 20)
set(gcf,'color','white')
xlim([1979 2023])
grid on


subplot(2,3,5)
plot(years_ret,ea_ret,'linewidth',4,'color', 'blue', 'LineStyle', '-', 'marker','none')
hold on 
plot(years_ret,y_values_fit_retea,'linewidth',4,'color','blue', 'LineStyle', '--', 'marker','none' )
title('East Antarctica', 'interpreter', 'latex')
ylim([190 310])
xlabel('Year', 'Fontsize', 14, 'interpreter', 'latex')
set(gca, 'yticklabel', {})
set(gca, 'FontSize', 20)
set(gcf,'color','white')
xlim([1979 2023])
grid on


subplot(2,3,6)
plot(years_ret,tot_ret,'linewidth',4,'color','black', 'LineStyle', '-', 'marker','none')
hold on
plot(years_ret,y_values_fit_rettot,'linewidth',4,'color','black', 'LineStyle', '--', 'marker','none' )
ylabel('Julian Day of Retreat', 'Fontsize', 14,'interpreter','latex')
set(gca, 'yticklabel', {})
title('Total', 'interpreter', 'latex')
ylim([190 310])
xlabel('Year', 'Fontsize', 14, 'interpreter', 'latex')
set(gca, 'FontSize', 20)
set(gcf,'color','white')
xlim([1979 2023])
grid on






plot(3) % Subplot of advance amplitude for [1] five regions and [2] entire continent 
plot(years_adv,ross_adv,'linewidth',4,'color','red', 'LineStyle', '-', 'marker','none')
hold on 
plot(years_adv,bell_adv,'linewidth',4,'color','#7E2F8E', 'LineStyle', '-', 'marker','none')
hold on 
plot(years_adv,wed_adv,'linewidth',4,'color','magenta', 'LineStyle', '-', 'marker','none')
hold on
plot(years_adv,kh_adv,'linewidth',4,'color','#77AC30', 'LineStyle', '-', 'marker','none')
hold on 
plot(years_adv,ea_adv,'linewidth',4,'color', 'blue', 'LineStyle', '-', 'marker','none')
hold on 
plot(years_adv,tot_adv,'linewidth',4,'color','black', 'LineStyle', '-', 'marker','none')

ylabel('Julian Day', 'Fontsize', 14, 'interpreter','latex')
xlabel('Year', 'Fontsize', 14, 'interpreter','latex')
set(gca, 'FontSize', 20)
set(gcf,'color','white')
xlim([1979 2024])
ylim([30 120])
grid on

title('Time of Advance', 'interpreter', 'latex')
    legend('Ross','Amundsen-Bellingshausen','Weddell','King Hakon','East Antarctica','Total', 'fontsize',12,'interpreter','latex')




plot(4) % Subplot of advance amplitude for [1] five regions and [2] entire continent 
plot(years_ret,ross_ret,'linewidth',4,'color','red', 'LineStyle', '-', 'marker','none')
hold on 
plot(years_ret,bell_ret,'linewidth',4,'color','#7E2F8E', 'LineStyle', '-', 'marker','none')
hold on 
plot(years_ret,wed_ret,'linewidth',4,'color','magenta', 'LineStyle', '-', 'marker','none')
hold on 
plot(years_ret,kh_ret,'linewidth',4,'color','#77AC30', 'LineStyle', '-', 'marker','none')
hold on 
plot(years_ret,ea_ret,'linewidth',4,'color', 'blue', 'LineStyle', '-', 'marker','none')
hold on
plot(years_ret,tot_ret,'linewidth',4,'color','black', 'LineStyle', '-', 'marker','none')

ylabel('Julian Day', 'Fontsize', 14, 'interpreter','latex')
xlabel('Year', 'Fontsize', 14, 'interpreter','latex')
set(gca, 'FontSize', 20)
set(gcf,'color','white')
xlim([1979 2023])
ylim([190 310])
grid on

title('Time of Retreat', 'interpreter', 'latex')
    legend('Ross','Amundsen-Bellingshausen','Weddell','King Hakon','East Antarctica','Total', 'fontsize',12,'interpreter','latex')




%% Write data to csv

% Combine year and time of advance into a matrix

toa_data = [years_adv', ross_adv', bell_adv',  wed_adv',  kh_adv', ea_adv', tot_adv'];
tor_data = [years_ret', ross_ret', bell_ret', wed_ret', kh_ret', ea_ret', tot_ret'];


toa_data_decyears = [years_adv', ross_adv_decyears', bell_adv_decyears',  wed_adv_decyears',  kh_adv_decyears', ea_adv_decyears', tot_adv_decyears'];
tor_data_decyears = [years_ret', ross_ret_decyears', bell_ret_decyears', wed_ret_decyears', kh_ret_decyears', ea_ret_decyears', tot_ret_decyears'];



% Write the data to a CSV file
% Decimal Julian days
csvwrite('advance_decimaljuliandays.csv', toa_data);
csvwrite('retreat_decimaljuliandays.csv', tor_data);

% Decimal years 
csvwrite('advance_decimalyears.csv', toa_data_decyears);
csvwrite('retreat_decimalyears.csv', tor_data_decyears);



% Create a README file
readme_text = [
    "This dataset contains the time of Antarctic sea ice " + ...
    "extent advance and retreat in decimal Julian days and decimal years. From 1979-2024 for day of advance. " + ...
    "From 1979-2023 for day of retreat."
    ""
    "Files included:"
    " 1. advance_decimaljuliandays.csv: Contains the following columns:"
    "    Column 1 - Years"
    "    Column 2 - Ross Time of Advance in Decimal Julian Days"
    "    Column 3 - Amundsen-Bellingshausen JTime of Advance in Decimal Julian Days"
    "    Column 4 - Weddell Time of Advance in Decimal Julian Days"
    "    Column 5 - King Hakon Time of Advance in Decimal Julian Days"
    "    Column 6 - East Antarctica Time of Advance in Decimal Julian Days"
    "    Column 7 - Total (Continental) Time of Advance in Decimal Julian Days"
    ""
    " 2. advance_decimalyears.csv: Contains the following columns:"
    "    Column 1 - Years"
    "    Column 2 - Ross Time of Advance in Decimal Years"
    "    Column 3 - Amundsen-Bellingshausen Time of Advance in Decimal Years"
    "    Column 4 - Weddell Time of Advance in Decimal Years"
    "    Column 5 - King Hakon Time of Advance in Decimal Years"
    "    Column 6 - East Antarctica Time of Advance in Decimal Years"
    "    Column 7 - Total (Continental) Time of Advance in Decimal Years"
    ""
    " 3. retreat_decimaljuliandays.csv: Contains the following columns:"
    "    Column 1 - Years"
    "    Column 2 - Ross Time of Retreat in Decimal Julian Days"
    "    Column 3 - Amundsen-Bellingshausen Time of Retreat in Decimal Julian Days"
    "    Column 4 - Weddell Time of Retreat in Decimal Julian Days"
    "    Column 5 - King Hakon Time of Retreat in Decimal Julian Days"
    "    Column 6 - East Antarctica Time of Retreat in Decimal Julian Days"
    "    Column 7 - Total (Continental) Time of Retreat in Decimal Julian Days"
    ""
    " 4. retreat_decimalyears.csv: Contains the following columns:"
    "    Column 1 - Years"
    "    Column 2 - Ross Time of Retreat in Decimal Years"
    "    Column 3 - Amundsen-Bellingshausen Time of Retreat in Decimal Years"
    "    Column 4 - Weddell Time of Retreat in Decimal Years"
    "    Column 5 - King Hakon Time of Retreat in Decimal Years"
    "    Column 6 - East Antarctica Time of Retreat in Decimal Years"
    "    Column 7 - Total (Continental) Time of Retreat in Decimal Years"

];

% Write the README file
fid = fopen('README.txt', 'wt');
fprintf(fid, '%s\n', readme_text{:});
fclose(fid);

% Optionally, package the files into a ZIP file
zip('ASIE_timeofadvance_timeofretreat_csv.zip', {'advance_decimaljuliandays.csv','advance_decimalyears.csv' , 'retreat_decimaljuliandays.csv', 'retreat_decimalyears.csv' , 'README.txt'});



%{
%% Write data to netCDF 
% Define the dimensions
years_dim = length(years_adv);

% Create NetCDF file for time of advance
toa_filename = 'dayofadvance.nc';
nccreate(toa_filename, 'years', 'Dimensions', {'years', years_dim});
nccreate(toa_filename, 'ross', 'Dimensions', {'ross', years_dim});
nccreate(toa_filename, 'bell', 'Dimensions', {'bell', years_dim});
nccreate(toa_filename, 'wed', 'Dimensions', {'wed', years_dim});
nccreate(toa_filename, 'kh', 'Dimensions', {'kh', years_dim});
nccreate(toa_filename, 'ea', 'Dimensions', {'ea', years_dim});
nccreate(toa_filename, 'tot', 'Dimensions', {'tot', years_dim});

% Write data to the variables
ncwrite(toa_filename, 'years', years_adv');
ncwrite(toa_filename, 'ross', ross_adv');
ncwrite(toa_filename, 'bell', bell_adv');
ncwrite(toa_filename, 'wed', wed_adv');
ncwrite(toa_filename, 'kh', kh_adv');
ncwrite(toa_filename, 'ea', ea_adv');
ncwrite(toa_filename, 'tot', tot_adv');


% Define the dimensions
years_dim = length(years_ret);

% Create NetCDF file for time of retreat
tor_filename = 'dayofretreat.nc';
nccreate(tor_filename, 'years', 'Dimensions', {'years', years_dim});
nccreate(tor_filename, 'ross', 'Dimensions', {'ross', years_dim});
nccreate(tor_filename, 'bell', 'Dimensions', {'bell', years_dim});
nccreate(tor_filename, 'wed', 'Dimensions', {'wed', years_dim});
nccreate(tor_filename, 'kh', 'Dimensions', {'kh', years_dim});
nccreate(tor_filename, 'ea', 'Dimensions', {'ea', years_dim});
nccreate(tor_filename, 'tot', 'Dimensions', {'tot', years_dim});

% Write data to the variables
ncwrite(tor_filename, 'years', years_ret');
ncwrite(tor_filename, 'ross', ross_ret');
ncwrite(tor_filename, 'bell', bell_ret');
ncwrite(tor_filename, 'wed', wed_ret');
ncwrite(tor_filename, 'kh', kh_ret');
ncwrite(tor_filename, 'ea', ea_ret');
ncwrite(tor_filename, 'tot', tot_ret');

% Create a README file
readme_text = [
    "This dataset contains the Julian day of Antarctic sea ice " + ...
    "extent advance and retreat. From 1979-2023 for day of advance. " + ...
    "From 1978-2023 for day of retreat."
    ""
    "Files included:"
    " 1. dayofadvance.nc: Contains the following columns:"
    "    - Years"
    "    - Ross Day of Advance"
    "    - Amundsen-Bellingshausen Day of Advance"
    "    - Weddell Day ofAdvance"
    "    - King Hakon Day of Advance"
    "    - East Antarctica Day of Advance"
    "    - Total (Continental) Day of Advance"
    ""
    " 2. dayofretreat.nc: Contains the following columns:"
    "    - Years"
    "    - Ross Day of Retreat"
    "    - Amundsen-Bellingshausen Day of Retreat"
    "    - Weddell Day of Retreat"
    "    - King Hakon Day of Retreat"
    "    - East Antarctica Day of Retreat"
    "    - Total (Continental) Day of Retreat"
];

% Write the README file
fid = fopen('README.txt', 'wt');
fprintf(fid, '%s\n', readme_text{:});
fclose(fid);

% Optionally, package the files into a ZIP file
zip('ASIE_dayofadvance_dayofretreat_nc.zip', {'timeofadvance.csv', 'timeofretreat.csv', 'README.txt'});

%}

