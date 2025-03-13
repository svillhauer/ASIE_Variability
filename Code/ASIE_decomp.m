clear all


% Initialize arrays for data storage

amplitude_adjusted_cycle_ross_all = [];
amplitude_adjusted_cycle_kh_all = [];
amplitude_adjusted_cycle_ea_all = []; 
amplitude_adjusted_cycle_wed_all = []; 
amplitude_adjusted_cycle_bell_all = [];
amplitude_adjusted_cycle_tot_all = [];

phase_adjusted_cycle_ross_all = []; 
phase_adjusted_cycle_kh_all = [];
phase_adjusted_cycle_ea_all = []; 
phase_adjusted_cycle_wed_all = [];
phase_adjusted_cycle_bell_all = [];
phase_adjusted_cycle_tot_all = []; 


xx_phase_adjusted_cycle_ross_all = []; 
xx_phase_adjusted_cycle_kh_all = [];
xx_phase_adjusted_cycle_ea_all = []; 
xx_phase_adjusted_cycle_wed_all = [];
xx_phase_adjusted_cycle_bell_all = [];
xx_phase_adjusted_cycle_tot_all = []; 


phase_years_ross_all = [];
phase_years_kh_all = [];
phase_years_ea_all = [];
phase_years_wed_all = [];
phase_years_bell_all = [];
phase_years_tot_all = [];


APAC_adjusted_cycle_ross_all = []; 
APAC_adjusted_cycle_kh_all = []; 
APAC_adjusted_cycle_ea_all = [];
APAC_adjusted_cycle_wed_all = [];
APAC_adjusted_cycle_bell_all = [];
APAC_adjusted_cycle_tot_all = [];

xx_APAC_adjusted_cycle_ross_all = []; 
xx_APAC_adjusted_cycle_kh_all = []; 
xx_APAC_adjusted_cycle_ea_all = [];
xx_APAC_adjusted_cycle_wed_all = [];
xx_APAC_adjusted_cycle_bell_all = [];
xx_APAC_adjusted_cycle_tot_all = [];

APAC_years_ross_all = [];
APAC_years_kh_all = [];
APAC_years_ea_all = [];
APAC_years_wed_all = [];
APAC_years_bell_all = [];
APAC_years_tot_all = [];


% Years in which the summer minimum was outside (below) the inter-decile range: 1984, 1993, 1997, 2017, 2018, 2022, 2023, 2024
% Year in which the winter maximum was below the inter-decile range: 1982, 1986, 1992, 2002, 2008, 2017, 2018, 2022, 2023
years = [1979:2022];
saveIterations = [4, 6, 8, 14, 15, 19, 24, 30, 39, 40, 44];

amplitude_adjusted_cycle = []; 
phase_adjusted_cycle =[];
APAC_cycle = []; 

% year_index = 40; % Where index 1 = 1978 and index 47 = 2024


%% Load ASIE data 
cd   '/Users/sarahvillhauer/Desktop/New ASIE/New Data'
addpath '/Users/sarahvillhauer/Desktop/New ASIE/Raw Data'

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


% Concatenate 2024 to reshaped data 
rossreshaped_2 = [rossreshaped ross_2024']; 
khreshaped_2 = [khreshaped kh_2024']; 
eareshaped_2 = [eareshaped ea_2024']; 
wedreshaped_2 = [wedreshaped wed_2024']; 
bellreshaped_2 = [bellreshaped bell_2024']; 
totreshaped_2 = [totreshaped tot_2024']; 
timereshaped_2 = [timereshaped time_2024']; 


%% Find the value of the minimum for each year (first day of the cycle)
[~, minIndex_ross] = min(rossreshaped);
[~, minIndex_kh] = min(khreshaped);
[~, minIndex_ea] = min(eareshaped);
[~, minIndex_wed] = min(wedreshaped);
[~, minIndex_bell] = min(bellreshaped);
[~, minIndex_tot] = min(totreshaped);



% Appending values of the proceeding year to the preceeding year 
%timenew = []; 
rossnew = []; 
khnew = []; 
eanew = []; 
wednew = []; 
bellnew = []; 
totnew = []; 


for i = 1:size(totreshaped_2,2)-1
   % timenew(:,i) = [timereshaped_2(:, i ) ; timereshaped_2(:,i+1)];
    rossnew(:,i) = [rossreshaped_2(:, i ) ; rossreshaped_2(:,i+1)];
    khnew(:,i) = [khreshaped_2(:, i ) ; khreshaped_2(:,i+1)];
    eanew(:,i) = [eareshaped_2(:, i ) ; eareshaped_2(:,i+1)];
    wednew(:,i) = [wedreshaped_2(:, i ) ; wedreshaped_2(:,i+1)];
    bellnew(:,i) = [bellreshaped_2(:, i ) ; bellreshaped_2(:,i+1)];
    totnew(:,i) = [totreshaped_2(:, i ) ; totreshaped_2(:,i+1)]; 
     
end

%% Define the SIE cycle for each year 
ross_SIE = []; 
kh_SIE = [];
ea_SIE = []; 
wed_SIE = []; 
bell_SIE = []; 
tot_SIE = []; 

%{
time_ross = [];
time_kh = [];
time_ea = [];
time_wed = []; 
time_bell = []; 
time_tot = []; 
%}


for i = 1:size(rossnew,2)
    ross_SIE(:,i) = rossnew(minIndex_ross(i):minIndex_ross(i)+364,i);
    kh_SIE(:,i) = khnew(minIndex_kh(i):minIndex_kh(i)+364,i);
    ea_SIE(:,i) = eanew(minIndex_ea(i):minIndex_ea(i)+364,i);
    wed_SIE(:,i) = wednew(minIndex_wed(i):minIndex_wed(i)+364,i);
    bell_SIE(:,i) = bellnew(minIndex_bell(i):minIndex_bell(i)+364,i);
    tot_SIE(:,i) = totnew(minIndex_tot(i):minIndex_tot(i)+364,i);

    %{
    time_ross(:,i) = timenew(minIndex_ross(i):minIndex_ross(i)+364,i);
    time_kh(:,i) = timenew(minIndex_kh(i):minIndex_kh(i)+364,i);
    time_ea(:,i) = timenew(minIndex_ea(i):minIndex_ea(i)+364,i);
    time_wed(:,i) = timenew(minIndex_wed(i):minIndex_wed(i)+364,i);
    time_bell(:,i) = timenew(minIndex_bell(i):minIndex_bell(i)+364,i);
    time_tot(:,i) = timenew(minIndex_tot(i):minIndex_tot(i)+364,i);
    %}
   
end

%% Traditional cycle (annual mean of Antarctic sea ice extent)

% Initialize a mask for NaN values
nanMask_ross = isnan(rossreshaped);
nanMask_kh = isnan(khreshaped);
nanMask_ea = isnan(eareshaped);
nanMask_wed = isnan(wedreshaped);
nanMask_bell = isnan(bellreshaped);
nanMask_tot = isnan(totreshaped);

% Compute the sum along rows ignoring NaN values
rowSums_ross = sum(rossreshaped, 2, 'omitnan');
rowSums_kh = sum(khreshaped, 2, 'omitnan');
rowSums_ea = sum(eareshaped, 2, 'omitnan');
rowSums_wed = sum(wedreshaped, 2, 'omitnan');
rowSums_bell = sum(bellreshaped, 2, 'omitnan');
rowSums_tot = sum(totreshaped, 2, 'omitnan');

% Count non-NaN elements along rows
counts_ross = sum(~nanMask_ross, 2);
counts_kh = sum(~nanMask_kh, 2);
counts_ea = sum(~nanMask_ea, 2);
counts_wed = sum(~nanMask_wed, 2);
counts_bell = sum(~nanMask_bell, 2);
counts_tot = sum(~nanMask_tot, 2);

% Calculate the mean row-wise
a_ross = rowSums_ross ./ counts_ross;
a_kh = rowSums_kh ./ counts_kh;
a_ea = rowSums_ea ./ counts_ea;
a_wed = rowSums_wed ./ counts_wed;
a_bell = rowSums_bell ./ counts_bell;
a_tot = rowSums_tot ./ counts_tot;

% Find the index of the minima for all regions 
[~, minIndex_ross] = min(a_ross);
[~, minIndex_kh] = min(a_kh);
[~, minIndex_ea] = min(a_ea);
[~, minIndex_wed] = min(a_wed);
[~, minIndex_bell] = min(a_bell);
[~, minIndex_tot] = min(a_tot);


% Reorganize a[s] so day-of-cycle = 0 corresponds to annual minimum 
a_ross_cycle = [a_ross(minIndex_ross:end) ; a_ross(1:minIndex_ross-1)];
a_kh_cycle = [a_kh(minIndex_kh:end) ; a_kh(1:minIndex_kh-1)];
a_ea_cycle = [a_ea(minIndex_ea:end) ; a_ea(1:minIndex_ea-1)];
a_wed_cycle = [a_wed(minIndex_wed:end) ; a_wed(1:minIndex_wed-1)];
a_bell_cycle = [a_bell(minIndex_bell:end) ; a_bell(1:minIndex_bell-1)];
a_tot_cycle = [a_tot(minIndex_tot:end) ; a_tot(1:minIndex_tot-1)];


day_of_cycle = [0:364]; 


figure;
yyaxis left
plot(day_of_cycle,a_ross_cycle,'linewidth',4,'color','red', 'LineStyle', '-', 'Marker', 'none')
hold on 
plot(day_of_cycle,a_bell_cycle,'linewidth',4, 'color', '#7E2F8E', 'LineStyle', '-', 'Marker', 'none')
hold on 
plot(day_of_cycle,a_wed_cycle,'linewidth',4,'color','magenta', 'LineStyle', '-', 'Marker', 'none')
hold on 
plot(day_of_cycle,a_kh_cycle,'linewidth',4,'color','#77AC30', 'LineStyle', '-', 'Marker', 'none')
hold on 
plot(day_of_cycle,a_ea_cycle, 'linewidth',4, 'color', 'blue', 'LineStyle', '-', 'Marker', 'none')
ylim([0 6])
ylabel('Regional Sea Ice Extent [$10^6$ km$^2$]', 'Fontsize', 14,'interpreter','latex')
set(gca, 'Color', 'white', 'XColor', 'black', 'YColor', 'black') % color for the plot area
grid on


yyaxis right
plot(day_of_cycle,a_tot_cycle,'linewidth',4,'color','black','Marker', 'none')
ylabel('Total Sea Ice Extent [$10^6$ km$^2$]', 'Fontsize', 14,'interpreter','latex')
xlabel('Day of Cycle', 'interpreter','latex')
xlim([0 364])
ylim([2 20])
title('Traditional Annual Cycle', 'Fontsize', 25 , 'interpreter','latex')
legend('Ross', 'Amundsen-Bellingshausen', 'Weddell', 'King Hakon','East Antarctica','Total','interpreter','latex')
set(gca, 'Color', 'white', 'XColor', 'black', 'YColor', 'black') % color for the plot area
set(gcf, 'Color', 'white') % Set figure background color to white
set(gca,'fontsize',20)
grid on


%% Invariant cycle (does not change from year to year) 
for year_index = 1:length(years)
% Define variables
s = day_of_cycle(:); 
extent_ross = a_ross; 
extent_kh = a_kh; 
extent_ea = a_ea; 
extent_wed = a_wed; 
extent_bell = a_bell; 
extent_tot = a_tot; 

[vI_ross,p_ross,V_ross,VAR_ross,CI_ross]  = csapsGCV(s,extent_ross); 
[vI_kh,p_kh,V_kh,VAR_kh,CI_kh]  = csapsGCV(s,extent_kh);
[vI_ea,p_ea,V_ea,VAR_ea,CI_ea]  = csapsGCV(s,extent_ea, [],[]);
[vI_wed,p_wed,V_wed,VAR_wed,CI_wed]  = csapsGCV(s,extent_wed, [],[]);
[vI_bell,p_bell,V_bell,VAR_bell,CI_bell]  = csapsGCV(s,extent_bell);
[vI_tot,p_tot,V_tot,VAR_tot,CI_tot]  = csapsGCV(s,extent_tot);


% Find the index of the minima for all regions 
[~, minIndex_ross] = min(vI_ross);
[~, minIndex_kh] = min(vI_kh);
[~, minIndex_ea] = min(vI_ea);
[~, minIndex_wed] = min(vI_wed);
[~, minIndex_bell] = min(vI_bell);
[~, minIndex_tot] = min(vI_tot);


% Reorganize a[s] so day-of-cycle = 0 corresponds to annual minimum 
vI_ross_cycle = [vI_ross(minIndex_ross:end) ; vI_ross(1:minIndex_ross-1)];
vI_kh_cycle = [vI_kh(minIndex_kh:end) ; vI_kh(1:minIndex_kh-1)];
vI_ea_cycle = [vI_ea(minIndex_ea:end) ; vI_ea(1:minIndex_ea-1)];
vI_wed_cycle = [vI_wed(minIndex_wed:end) ; vI_wed(1:minIndex_wed-1)];
vI_bell_cycle = [vI_bell(minIndex_bell:end) ; vI_bell(1:minIndex_bell-1)];
vI_tot_cycle = [vI_tot(minIndex_tot:end) ; vI_tot(1:minIndex_tot-1)];


%{
figure;
subplot(2,3,1)
plot(day_of_cycle,a_ross_cycle,'linewidth',4,'color','black', 'LineStyle', ':', 'Marker', 'none')
hold on 
plot(day_of_cycle,ross_SIE(:,year_index),'linewidth',4,'color','blue', 'LineStyle', ':', 'Marker', 'none')
hold on 
plot(day_of_cycle,vI_ross_cycle,'linewidth',2,'color','green', 'LineStyle', '-', 'Marker', 'none')
%xlabel('Day of Cycle','fontsize',25)
set(gca,'XTickLabel',[])
ylabel('Regional Sea Ice Extent [$10^{6}$ $\mathrm{km}^{2}$]','fontsize',25,'interpreter', 'latex')
title('Ross','fontsize',25,'interpreter','latex')
xlim([0 364])
set(gcf,'color','white')
ylim([0 6])
set(gca,'fontsize',20)
grid on 

subplot(2,3,2)
plot(day_of_cycle,a_bell_cycle,'linewidth',4,'color','black', 'LineStyle', ':', 'Marker', 'none')
hold on 
plot(day_of_cycle,bell_SIE(:,year_index),'linewidth',4,'color','blue', 'LineStyle', ':', 'Marker', 'none')
hold on 
plot(day_of_cycle,vI_bell_cycle,'linewidth',2,'color','green', 'LineStyle', '-', 'Marker', 'none')
%xlabel('Day of Cycle','fontsize',25)
set(gca,'XTickLabel',[])
%ylabel('Regional Sea Ice Extent [10^6 km^2]','fontsize',25)
set(gca,'YTickLabel',[])
title('Amundsen-Bellingshausen','fontsize',25, 'interpreter','latex')
xlim([0 364])
set(gcf,'color','white')
ylim([0 6])
set(gca,'fontsize',20)
grid on

subplot(2,3,3)
plot(day_of_cycle,a_wed_cycle,'linewidth',4,'color','black', 'LineStyle', ':', 'Marker', 'none')
hold on 
plot(day_of_cycle,wed_SIE(:,year_index),'linewidth',4,'color','blue', 'LineStyle', ':', 'Marker', 'none')
hold on 
plot(day_of_cycle,vI_wed_cycle,'linewidth',2,'color','green', 'LineStyle', '-', 'Marker', 'none')
%xlabel('Day of Cycle','fontsize',25)
set(gca,'XTickLabel',[])
%ylabel('Regional Sea Ice Extent [10^6 km^2]','fontsize',25)
set(gca,'YTickLabel',[])
title('Weddell','fontsize',25, 'interpreter','latex')
xlim([0 364])
set(gcf,'color','white')
ylim([0 6])
set(gca,'fontsize',20)
grid on 

subplot(2,3,4)
plot(day_of_cycle,a_kh_cycle,'linewidth',4,'color','black', 'LineStyle', ':', 'Marker', 'none')
hold on 
plot(day_of_cycle,kh_SIE(:,year_index),'linewidth',4,'color','blue', 'LineStyle', ':', 'Marker', 'none')
hold on 
plot(day_of_cycle,vI_kh_cycle,'linewidth',2,'color','green', 'LineStyle', '-', 'Marker', 'none')
xlabel('Day of Cycle','fontsize',25, 'interpreter','latex')
ylabel('Regional Sea Ice Extent [$10^{6}$ $\mathrm{km}^{2}$]','fontsize',25, 'interpreter','latex')
title('King Hakon','fontsize',25, 'interpreter','latex')
xlim([0 364])
set(gcf,'color','white')
set(gca,'fontsize',20)
ylim([0 6])
grid on 

subplot(2,3,5)
plot(day_of_cycle,a_ea_cycle,'linewidth',4,'color','black', 'LineStyle', ':', 'Marker', 'none')
hold on 
plot(day_of_cycle,ea_SIE(:,year_index),'linewidth',4,'color','blue', 'LineStyle', ':', 'Marker', 'none')
hold on 
plot(day_of_cycle,vI_ea_cycle,'linewidth',2,'color','green', 'LineStyle', '-', 'Marker', 'none')
xlabel('Day of Cycle','fontsize',25,'interpreter','latex')
%ylabel('Regional Sea Ice Extent [10^6 km^2]','fontsize',25)
set(gca,'YTickLabel',[])
title('East Antarctica','fontsize',25,'interpreter','latex')
xlim([0 364])
set(gcf,'color','white')
set(gca,'fontsize',20)
ylim([0 6])
grid on 

subplot(2,3,6)
plot(day_of_cycle,a_tot_cycle,'linewidth',4,'color','black', 'LineStyle', ':', 'Marker', 'none')
hold on 
plot(day_of_cycle,tot_SIE(:,year_index),'linewidth',4,'color','blue', 'LineStyle', ':', 'Marker', 'none')
hold on 
plot(day_of_cycle,vI_tot_cycle,'linewidth',2,'color','green', 'LineStyle', '-', 'Marker', 'none')
xlabel('Day of Cycle','fontsize',25, 'interpreter','latex')
ylabel('Total Sea Ice Extent [$10^{6}$ $\mathrm{km}^{2}$]','fontsize',25,'interpreter','latex')
title('Total','fontsize',25,'interpreter','latex')
xlim([0 364])
set(gcf,'color','white')
set(gca,'fontsize',20)
year_value = years(year_index);
year_string = num2str(year_value);
legend('Traditional', ['Recorded in ', year_string],'Invariant', 'interpreter','latex')
grid on 



%}
%{
 % Check if the current iteration is in the list of saveIterations
    if ismember(i, saveIterations)
        % Construct the filename
        filename = sprintf('invariant_iteration_%d.fig', i);
        % Save the figure
        savefig(filename);
    end
    
    % Close the figure to avoid having too many open figures
    close;

%}

%{
figure(3)
plot(day_of_cycle,a_ross_cycle,'linewidth',4,'color','blue', 'LineStyle', '-', 'Marker', 'none')
hold on 
plot(day_of_cycle,ross_SIE(:,year_index),'linewidth',4,'color','black', 'LineStyle', '-', 'Marker', 'none')
hold on 
plot(day_of_cycle,vI_ross_cycle,'linewidth',2,'color','red', 'LineStyle', '-', 'Marker', 'none')
xlabel('Day of Cycle','fontsize',25)
ylabel('Regional Sea Ice Extent [10^6 km^2]','fontsize',25)
title('Ross','fontsize',25)
xlim([0 365])
set(gcf,'color','white')
set(gca,'fontsize',20)
legend('Traditional', 'Recorded in 2016','Invariant')
grid on 


figure(4)
plot(day_of_cycle,a_kh_cycle,'linewidth',4,'color','blue', 'LineStyle', '-', 'Marker', 'none')
hold on 
plot(day_of_cycle,kh_SIE(:,year_index),'linewidth',4,'color','black', 'LineStyle', '-', 'Marker', 'none')
hold on 
plot(day_of_cycle,vI_kh_cycle,'linewidth',2,'color','red', 'LineStyle', '-', 'Marker', 'none')
xlabel('Day of Cycle','fontsize',25)
ylabel('Regional Sea Ice Extent [10^6 km^2]','fontsize',25)
title('King Hakon','fontsize',25)
xlim([0 365])
set(gcf,'color','white')
set(gca,'fontsize',20)
legend('Traditional', 'Recorded in 2016','Invariant')
grid on 

figure(5)
plot(day_of_cycle,a_ea_cycle,'linewidth',4,'color','blue', 'LineStyle', '-', 'Marker', 'none')
hold on 
plot(day_of_cycle,ea_SIE(:,year_index),'linewidth',4,'color','black', 'LineStyle', '-', 'Marker', 'none')
hold on 
plot(day_of_cycle,vI_ea_cycle,'linewidth',2,'color','red', 'LineStyle', '-', 'Marker', 'none')
xlabel('Day of Cycle','fontsize',25)
ylabel('Regional Sea Ice Extent [10^6 km^2]','fontsize',25)
title('East Antarctica','fontsize',25)
xlim([0 365])
set(gcf,'color','white')
set(gca,'fontsize',20)
legend('Traditional', 'Recorded in 2016','Invariant')
grid on 


figure(6)
plot(day_of_cycle,a_wed_cycle,'linewidth',4,'color','blue', 'LineStyle', '-', 'Marker', 'none')
hold on 
plot(day_of_cycle,wed_SIE(:,year_index),'linewidth',4,'color','black', 'LineStyle', '-', 'Marker', 'none')
hold on 
plot(day_of_cycle,vI_wed_cycle,'linewidth',2,'color','red', 'LineStyle', '-', 'Marker', 'none')
xlabel('Day of Cycle','fontsize',25)
ylabel('Regional Sea Ice Extent [10^6 km^2]','fontsize',25)
title('Weddell','fontsize',25)
xlim([0 365])
set(gcf,'color','white')
set(gca,'fontsize',20)
legend('Traditional', 'Recorded in 2016','Invariant')
grid on 


figure(7)
plot(day_of_cycle,a_bell_cycle,'linewidth',4,'color','blue', 'LineStyle', '-', 'Marker', 'none')
hold on 
plot(day_of_cycle,bell_SIE(:,year_index),'linewidth',4,'color','black', 'LineStyle', '-', 'Marker', 'none')
hold on 
plot(day_of_cycle,vI_bell_cycle,'linewidth',2,'color','red', 'LineStyle', '-', 'Marker', 'none')
xlabel('Day of Cycle','fontsize',25)
ylabel('Regional Sea Ice Extent [10^6 km^2]','fontsize',25)
title('Bellinghausen','fontsize',25)
xlim([0 365])
set(gcf,'color','white')
set(gca,'fontsize',20)
legend('Traditional', 'Recorded in 2016','Invariant')
grid on


figure(8)
plot(day_of_cycle,a_tot_cycle,'linewidth',4,'color','blue', 'LineStyle', '-', 'Marker', 'none')
hold on 
plot(day_of_cycle,tot_SIE(:,year_index),'linewidth',4,'color','black', 'LineStyle', '-', 'Marker', 'none')
hold on 
plot(day_of_cycle,vI_tot_cycle,'linewidth',2,'color','red', 'LineStyle', '-', 'Marker', 'none')
xlabel('Day of Cycle','fontsize',25)
ylabel('Total Sea Ice Extent [10^6 km^2]','fontsize',25)
xlim([0 365])
set(gcf,'color','white')
set(gca,'fontsize',20)
legend('Traditional', 'Recorded in 2016','Invariant')
grid on 
%}


%% Amplitude-adjusted cycle 

% Define variables
s = day_of_cycle(:); 

% Scale inputs by annual maxima and minimum
for j = 1:length(years) %2:length(years)-1
    input_ross(:,j) = (rossreshaped(:,j)-min(rossreshaped(:,j)))/(max(rossreshaped(:,j))-min(rossreshaped(:,j)));
    input_kh(:,j) = (khreshaped(:,j)-min(khreshaped(:,j)))/(max(khreshaped(:,j))-min(khreshaped(:,j))); 
    input_ea(:,j) = (eareshaped(:,j)-min(eareshaped(:,j)))/(max(eareshaped(:,j))-min(rossreshaped(:,j))); 
    input_wed(:,j) = (wedreshaped(:,j)-min(wedreshaped(:,j)))/(max(wedreshaped(:,j))-min(wedreshaped(:,j))); 
    input_bell(:,j) = (bellreshaped(:,j)-min(bellreshaped(:,j)))/(max(bellreshaped(:,j))-min(bellreshaped(:,j))); 
    input_tot(:,j) = (totreshaped(:,j)-min(totreshaped(:,j)))/(max(totreshaped(:,j))-min(totreshaped(:,j))); 
end

%{
input_ross = [rossreshaped(:,1) input_ross rossreshaped(:,47)]; 
input_kh = [khreshaped(:,1) input_kh khreshaped(:,47)]; 
input_ea = [eareshaped(:,1) input_ea eareshaped(:,47)]; 
input_wed = [wedreshaped(:,1) input_wed wedreshaped(:,47)]; 
input_bell = [bellreshaped(:,1) input_bell bellreshaped(:,47)]; 
input_tot = [totreshaped(:,1) input_tot totreshaped(:,47)]; 
%}

% Initialize a mask for NaN values
nanMask_ross = isnan(input_ross);
nanMask_kh = isnan(input_kh);
nanMask_ea = isnan(input_ea);
nanMask_wed = isnan(input_wed);
nanMask_bell = isnan(input_bell);
nanMask_tot = isnan(input_tot);


% Compute the sum along rows ignoring NaN values
rowSums_ross = sum(input_ross, 2, 'omitnan');
rowSums_kh = sum(input_kh, 2, 'omitnan');
rowSums_ea = sum(input_ea, 2, 'omitnan');
rowSums_wed = sum(input_wed, 2, 'omitnan');
rowSums_bell = sum(input_bell, 2, 'omitnan');
rowSums_tot = sum(input_tot, 2, 'omitnan');

% Count non-NaN elements along rows
counts_ross = sum(~nanMask_ross, 2);
counts_kh = sum(~nanMask_kh, 2);
counts_ea = sum(~nanMask_ea, 2);
counts_wed = sum(~nanMask_wed, 2);
counts_bell = sum(~nanMask_bell, 2);
counts_tot = sum(~nanMask_tot, 2);

% Calculate the mean row-wise
input_ross = rowSums_ross ./ counts_ross;
input_kh = rowSums_kh ./ counts_kh;
input_ea = rowSums_ea ./ counts_ea;
input_wed = rowSums_wed ./ counts_wed;
input_bell = rowSums_bell ./ counts_bell;
input_tot = rowSums_tot ./ counts_tot;







[v_ross,p_ross,V_ross,VAR_ross,CI_ross]  = csapsGCV(s,input_ross); 
[v_kh,p_kh,V_kh,VAR_kh,CI_kh]  = csapsGCV(s,input_kh);
[v_ea,p_ea,V_ea,VAR_ea,CI_ea]  = csapsGCV(s,input_ea);
[v_wed,p_wed,V_wed,VAR_wed,CI_wed]  = csapsGCV(s,input_wed);
[v_bell,p_bell,V_bell,VAR_bell,CI_bell]  = csapsGCV(s,input_bell);
[v_tot,p_tot,V_tot,VAR_tot,CI_tot]  = csapsGCV(s,input_tot);



aA_ross = v_ross*(max(rossreshaped(:,year_index))-min(rossreshaped(:,year_index)))+min(rossreshaped(:,year_index));
aA_kh = v_kh*(max(khreshaped(:,year_index))-min(khreshaped(:,year_index)))+min(khreshaped(:,year_index)); 
aA_ea = v_ea*(max(eareshaped(:,year_index))-min(eareshaped(:,year_index)))+min(eareshaped(:,year_index)); 
aA_wed = v_wed*(max(wedreshaped(:,year_index))-min(wedreshaped(:,year_index)))+min(wedreshaped(:,year_index)); 
aA_bell = v_bell*(max(bellreshaped(:,year_index))-min(bellreshaped(:,year_index)))+min(bellreshaped(:,year_index)); 
aA_tot = v_tot*(max(totreshaped(:,year_index))-min(totreshaped(:,year_index)))+min(totreshaped(:,year_index)); 


% Find the index of the minima for all regions 
[~, minIndex_ross] = min(aA_ross);
[~, minIndex_kh] = min(aA_kh);
[~, minIndex_ea] = min(aA_ea);
[~, minIndex_wed] = min(aA_wed);
[~, minIndex_bell] = min(aA_bell);
[~, minIndex_tot] = min(aA_tot);


% Reorganize a[s] so day-of-cycle = 0 corresponds to annual minimum 
aA_ross_cycle = [aA_ross(minIndex_ross:end) ; aA_ross(1:minIndex_ross-1)];
aA_kh_cycle = [aA_kh(minIndex_kh:end) ; aA_kh(1:minIndex_kh-1)];
aA_ea_cycle = [aA_ea(minIndex_ea:end) ; aA_ea(1:minIndex_ea-1)];
aA_wed_cycle = [aA_wed(minIndex_wed:end) ; aA_wed(1:minIndex_wed-1)];
aA_bell_cycle = [aA_bell(minIndex_bell:end) ; aA_bell(1:minIndex_bell-1)];
aA_tot_cycle = [aA_tot(minIndex_tot:end) ; aA_tot(1:minIndex_tot-1)];


%{
figure(9)
plot(day_of_cycle,a_ross_cycle,'linewidth',4,'color','red', 'LineStyle', '-', 'Marker', 'none')
hold on 
plot(day_of_cycle,ross_SIE(:,year_index),'linewidth',4,'color','black', 'LineStyle', '-', 'Marker', 'none')
hold on 
plot(day_of_cycle,aA_ross_cycle,'linewidth',2,'color','green', 'LineStyle', '-', 'Marker', 'none')
xlabel('Day of Cycle','fontsize',25)
ylabel('Regional Sea Ice Extent [10^6 km^2]','fontsize',25)
title('Ross','fontsize',25)
xlim([0 365])
set(gcf,'color','white')
set(gca,'fontsize',20)
legend('Traditional', 'Recorded in 2016','Amplitude-adjusted')
grid on 


figure(10)
plot(day_of_cycle,a_kh_cycle,'linewidth',4,'color','red', 'LineStyle', '-', 'Marker', 'none')
hold on 
plot(day_of_cycle,kh_SIE(:,year_index),'linewidth',4,'color','black', 'LineStyle', '-', 'Marker', 'none')
hold on 
plot(day_of_cycle,aA_kh_cycle,'linewidth',2,'color','green', 'LineStyle', '-', 'Marker', 'none')
xlabel('Day of Cycle','fontsize',25)
ylabel('Regional Sea Ice Extent [10^6 km^2]','fontsize',25)
title('King Hakon','fontsize',25)
xlim([0 365])
set(gcf,'color','white')
set(gca,'fontsize',20)
legend('Traditional', 'Recorded in 2016','Amplitude-adjusted')
grid on 

figure(11)
plot(day_of_cycle,a_ea_cycle,'linewidth',4,'color','red', 'LineStyle', '-', 'Marker', 'none')
hold on 
plot(day_of_cycle,ea_SIE(:,year_index),'linewidth',4,'color','black', 'LineStyle', '-', 'Marker', 'none')
hold on 
plot(day_of_cycle,aA_ea_cycle,'linewidth',2,'color','green', 'LineStyle', '-', 'Marker', 'none')
xlabel('Day of Cycle','fontsize',25)
ylabel('Regional Sea Ice Extent [10^6 km^2]','fontsize',25)
title('East Antarctica','fontsize',25)
xlim([0 365])
set(gcf,'color','white')
set(gca,'fontsize',20)
legend('Traditional', 'Recorded in 2016','Amplitude-adjusted')
grid on 


figure(12)
plot(day_of_cycle,a_wed_cycle,'linewidth',4,'color','red', 'LineStyle', '-', 'Marker', 'none')
hold on 
plot(day_of_cycle,wed_SIE(:,year_index),'linewidth',4,'color','black', 'LineStyle', '-', 'Marker', 'none')
hold on 
plot(day_of_cycle,aA_wed_cycle,'linewidth',2,'color','green', 'LineStyle', '-', 'Marker', 'none')
xlabel('Day of Cycle','fontsize',25)
ylabel('Regional Sea Ice Extent [10^6 km^2]','fontsize',25)
title('Weddell','fontsize',25)
xlim([0 365])
set(gcf,'color','white')
set(gca,'fontsize',20)
legend('Traditional', 'Recorded in 2016','Amplitude-adjusted')
grid on 


figure(13)
plot(day_of_cycle,a_bell_cycle,'linewidth',4,'color','red', 'LineStyle', '-', 'Marker', 'none')
hold on 
plot(day_of_cycle,bell_SIE(:,year_index),'linewidth',4,'color','black', 'LineStyle', '-', 'Marker', 'none')
hold on 
plot(day_of_cycle,aA_bell_cycle,'linewidth',2,'color','green', 'LineStyle', '-', 'Marker', 'none')
xlabel('Day of Cycle','fontsize',25)
ylabel('Regional Sea Ice Extent [10^6 km^2]','fontsize',25)
title('Amundsen-Bellinghausen','fontsize',25)
xlim([0 365])
set(gcf,'color','white')
set(gca,'fontsize',20)
legend('Traditional', 'Recorded in 2016','Amplitude-adjusted')
grid on


figure(14)
plot(day_of_cycle,a_tot_cycle,'linewidth',4,'color','red', 'LineStyle', '-', 'Marker', 'none')
hold on 
plot(day_of_cycle,tot_SIE(:,year_index),'linewidth',4,'color','black', 'LineStyle', '-', 'Marker', 'none')
hold on 
plot(day_of_cycle,aA_tot_cycle,'linewidth',2,'color','green', 'LineStyle', '-', 'Marker', 'none')
xlabel('Day of Cycle','fontsize',25)
ylabel('Total Sea Ice Extent [10^6 km^2]','fontsize',25)
xlim([0 365])
set(gcf,'color','white')
set(gca,'fontsize',20)
legend('Traditional', 'Recorded in 2016','Amplitude-adjusted')
grid on 
%}

%% Phase-adjusted cycle 


% Load time of advance and time of retreat data
toa = readtable('advance_decimalyears.csv');
toa = table2array(toa);

tor = readtable('retreat_decimalyears.csv');
tor = table2array(tor);

min_extent_day_ross = toa(:,2); 
min_extent_day_bell = toa(:,3);
min_extent_day_wed = toa(:,4); 
min_extent_day_kh = toa(:,5); 
min_extent_day_ea = toa(:,6); 
min_extent_day_tot = toa(:,7);

max_extent_day_ross = tor(:,2); 
max_extent_day_bell = tor(:,3);
max_extent_day_wed = tor(:,4); 
max_extent_day_kh = tor(:,5); 
max_extent_day_ea = tor(:,6); 
max_extent_day_tot = tor(:,7);

% Define phase(t) (AKA xx) 
xx_ross = zeros(1,365); 
xx_kh = zeros(1,365); 
xx_ea = zeros(1,365); 
xx_bell = zeros(1,365); 
xx_wed = zeros(1,365); 
xx_tot = zeros(1,365); 

% Initialize variables to store the first and last indices
ross_first_index = [];
ross_last_index = [];
kh_first_index = [];
kh_last_index = [];
ea_first_index = [];
ea_last_index = [];
wed_first_index = [];
wed_last_index = [];
bell_first_index = [];
bell_last_index = [];
tot_first_index = [];
tot_last_index = [];

year = year_index;  

for i = 1:365
    if timereshaped(i,year) >= min_extent_day_ross(year) && timereshaped(i,year) <= max_extent_day_ross(year)

        xx_ross(i) = (timereshaped(i,year) - min_extent_day_ross(year))/(max_extent_day_ross(year+1)- min_extent_day_ross(year));
  
        % Save the first index where the condition is met
        if isempty(ross_first_index)
            ross_first_index = i;
       
        end
        
        % Update the last index where the condition is met
        ross_last_index = i;
   
    else
        xx_ross(i) = 0;



    end
  
end


for i = 1:365
     if timereshaped(i,year) >= min_extent_day_kh(year) && timereshaped(i,year) <= max_extent_day_kh(year)

    
        xx_kh(i) = (timereshaped(i,year) - min_extent_day_kh(year))/(max_extent_day_kh(year+1)- min_extent_day_kh(year));
       
        % Save the first index where the condition is met
        if isempty(kh_first_index)
         
            kh_first_index = i;
       
        end
        
        % Update the last index where the condition is met
  
        kh_last_index = i;
 
    else
 
        xx_kh(i) = 0;
    


    end
       
end



for i = 1:365
     if timereshaped(i,year) >= min_extent_day_ea(year) && timereshaped(i,year) <= max_extent_day_ea(year)

    
        xx_ea(i) = (timereshaped(i,year) - min_extent_day_ea(year))/(max_extent_day_ea(year+1)- min_extent_day_ea(year));
      
        % Save the first index where the condition is met
        if isempty(ea_first_index)
         
            ea_first_index = i;
       
        end
        
        % Update the last index where the condition is met
  
        ea_last_index = i;
 
    else
 
        xx_ea(i) = 0;
    


    end
      
end


for i = 1:365
     if timereshaped(i,year) >= min_extent_day_wed(year) && timereshaped(i,year) <= max_extent_day_wed(year)

    
        xx_wed(i) = (timereshaped(i,year) - min_extent_day_wed(year))/(max_extent_day_wed(year+1)- min_extent_day_wed(year));
       
        % Save the first index where the condition is met
        if isempty(wed_first_index)
         
            wed_first_index = i;
       
        end
        
        % Update the last index where the condition is met
  
        wed_last_index = i;
 
    else
 
        xx_wed(i) = 0;
    


    end
       
end



for i = 1:365
     if timereshaped(i,year) >= min_extent_day_bell(year) && timereshaped(i,year) <= max_extent_day_bell(year)

    
        xx_bell(i) = (timereshaped(i,year) - min_extent_day_bell(year))/(max_extent_day_bell(year+1)- min_extent_day_bell(year));
       
        % Save the first index where the condition is met
        if isempty(bell_first_index)
         
            bell_first_index = i;
       
        end
        
        % Update the last index where the condition is met
  
        bell_last_index = i;
 
    else
 
        xx_bell(i) = 0;
    


    end
       
end



for i = 1:365
     if timereshaped(i,year) >= min_extent_day_tot(year) && timereshaped(i,year) <= max_extent_day_tot(year)

    
        xx_tot(i) = (timereshaped(i,year) - min_extent_day_tot(year))/(max_extent_day_tot(year+1)- min_extent_day_tot(year));
       
        
        % Save the first index where the condition is met
        if isempty(tot_first_index)
         
            tot_first_index = i;
       
        end
        
        % Update the last index where the condition is met
  
        tot_last_index = i;
 
    else
 
        xx_tot(i) = 0;
    


    end
        
end




[vP_ross,p_ross,V_ross,VAR_ross,CI_ross,beta1_ross,beta2_ross,xx_ross]  = csapsGCV_with_Phase(s,a_ross, [], xx_ross); 
[vP_kh,p_kh,V_kh,VAR_kh,CI_kh,beta1_kh,beta2_kh,xx_kh]  = csapsGCV_with_Phase(s,a_kh, [], xx_kh);
[vP_ea,p_ea,V_ea,VAR_ea,CI_ea,beta1_ea,beta2_ea,xx_ea]  = csapsGCV_with_Phase(s,a_ea, [], xx_ea);
[vP_wed,p_wed,V_wed,VAR_wed,CI_wed,beta1_wed,beta2_wed,xx_wed]  = csapsGCV_with_Phase(s,a_wed, [], xx_wed);
[vP_bell,p_bell,V_bell,VAR_bell,CI_bell,beta1_bell,beta2_bell,xx_bell]  = csapsGCV_with_Phase(s,a_bell, [], xx_bell);
[vP_tot,p_tot,V_tot,VAR_tot,CI_tot,beta1_tot,beta2_tot,xx_tot]  = csapsGCV_with_Phase(s,a_tot, [], xx_tot);


% Ross
% Re-perform GCV with extrapolated matrices (define smooth function of time)  
m_ross =  diff(xx_ross(ross_last_index-1:ross_last_index)); 
new_values_first = xx_ross(ross_first_index):-m_ross:0; 
new_values_first = [new_values_first 0]; 
new_values_first = flip(new_values_first); 
new_values_end = xx_ross(ross_last_index):m_ross:364;


if length(xx_ross(1:ross_first_index-1)) < length(new_values_first)
    xx_ross_new = [new_values_first(1:end), xx_ross(ross_first_index:ross_last_index), new_values_end(2:end)];
    % Re-define s for phase adjustment
    s_ross = 0:length(xx_ross_new)-1;
    cutoff_ross1 = length(new_values_first(1:end-1)) - length(a_ross(1:ross_first_index-1));
    cutoff_ross2 = length(new_values_end(2:end))  - length(a_ross(ross_last_index+1:end));
    a_ross_phase = [a_ross(end-cutoff_ross1:end); a_ross ; a_ross(1:cutoff_ross2)];
    [vP_ross,p_ross,V_ross,VAR_ross,CI_ross]  = csapsGCV(s_ross,a_ross_phase, [], xx_ross_new);






else
    xx_ross_new = [xx_ross(1:ross_last_index), new_values_end(2:end)];
    xx_ross_new(ross_first_index-length(new_values_first(1:end-1)):ross_first_index-1) = new_values_first(1:end-1);
    % Re-define s for phase adjustment
    s_ross = 0:length(xx_ross_new)-1;
    cutoff_ross = length(new_values_end(2:end))  - length(a_ross(ross_last_index+1:end));
    a_ross_phase = [a_ross ; a_ross(1:cutoff_ross)  ];
    [vP_ross,p_ross,V_ross,VAR_ross,CI_ross]  = csapsGCV(s_ross,a_ross_phase, [], xx_ross_new);


    index_ross = find(xx_ross_new, 1);
    xx_ross_new = xx_ross_new(index_ross-1:end);
    vP_ross = vP_ross(index_ross-1:end);


end








% King Hakon
% Re-perform GCV with extrapolated matrices (define smooth function of time)  
m_kh =  diff(xx_kh(kh_last_index-1:kh_last_index)); 
new_values_first = xx_kh(kh_first_index):-m_kh:0; 
new_values_first = [new_values_first 0]; 
new_values_first = flip(new_values_first); 
new_values_end = xx_kh(kh_last_index):m_kh:364;


if length(xx_kh(1:kh_first_index-1)) < length(new_values_first)
    xx_kh_new = [new_values_first(1:end), xx_kh(kh_first_index:kh_last_index), new_values_end(2:end)];
    % Re-define s for phase adjustment
    s_kh = 0:length(xx_kh_new)-1;
    cutoff_kh1 = length(new_values_first(1:end-1)) - length(a_kh(1:kh_first_index-1));
    cutoff_kh2 = length(new_values_end(2:end))  - length(a_kh(kh_last_index+1:end));
    a_kh_phase = [a_kh(end-cutoff_kh1:end); a_kh ; a_kh(1:cutoff_kh2)];
    [vP_kh,p_kh,V_kh,VAR_kh,CI_kh]  = csapsGCV(s_kh,a_kh_phase, [], xx_kh_new);






else
    xx_kh_new = [xx_kh(1:kh_last_index), new_values_end(2:end)];
    xx_kh_new(kh_first_index-length(new_values_first(1:end-1)):kh_first_index-1) = new_values_first(1:end-1);
    % Re-define s for phase adjustment
    s_kh = 0:length(xx_kh_new)-1;
    cutoff_kh = length(new_values_end(2:end))  - length(a_kh(kh_last_index+1:end));
    a_kh_phase = [a_kh ; a_kh(1:cutoff_kh)  ];
    [vP_kh,p_kh,V_kh,VAR_kh,CI_kh]  = csapsGCV(s_kh,a_kh_phase, [], xx_kh_new);


    index_kh = find(xx_kh_new, 1);
    xx_kh_new = xx_kh_new(index_kh-1:end);
    vP_kh = vP_kh(index_kh-1:end);


end



% East Antarctica
% Re-perform GCV with extrapolated matrices (define smooth function of time)  
m_ea =  diff(xx_ea(ea_last_index-1:ea_last_index)); 
new_values_first = xx_ea(ea_first_index):-m_ea:0; 
new_values_first = [new_values_first 0]; 
new_values_first = flip(new_values_first); 
new_values_end = xx_ea(ea_last_index):m_ea:364;


if length(xx_ea(1:ea_first_index-1)) < length(new_values_first)
    xx_ea_new = [new_values_first(1:end), xx_ea(ea_first_index:ea_last_index), new_values_end(2:end)];
    % Re-define s for phase adjustment
    s_ea = 0:length(xx_ea_new)-1;
    cutoff_ea1 = length(new_values_first(1:end-1)) - length(a_ea(1:ea_first_index-1));
    cutoff_ea2 = length(new_values_end(2:end))  - length(a_ea(ea_last_index+1:end));
    a_ea_phase = [a_ea(end-cutoff_ea1:end); a_ea ; a_ea(1:cutoff_ea2)];
    [vP_ea,p_ea,V_ea,VAR_ea,CI_ea]  = csapsGCV(s_ea,a_ea_phase, [], xx_ea_new);






else
    xx_ea_new = [xx_ea(1:ea_last_index), new_values_end(2:end)];
    xx_ea_new(ea_first_index-length(new_values_first(1:end-1)):ea_first_index-1) = new_values_first(1:end-1);
    % Re-define s for phase adjustment
    s_ea = 0:length(xx_ea_new)-1;
    cutoff_ea = length(new_values_end(2:end))  - length(a_ea(ea_last_index+1:end));
    a_ea_phase = [a_ea ; a_ea(1:cutoff_ea)  ];
    [vP_ea,p_ea,V_ea,VAR_ea,CI_ea]  = csapsGCV(s_ea,a_ea_phase, [], xx_ea_new);

    index_ea = find(xx_ea_new, 1);
    xx_ea_new = xx_ea_new(index_ea-1:end);
    vP_ea = vP_ea(index_ea-1:end);


end


% Weddell
% Re-perform GCV with extrapolated matrices (define smooth function of time)  
m_wed =  diff(xx_wed(wed_last_index-1:wed_last_index)); 
new_values_first = xx_wed(wed_first_index):-m_wed:0; 
new_values_first = [new_values_first 0]; 
new_values_first = flip(new_values_first); 
new_values_end = xx_wed(wed_last_index):m_wed:364;


if length(xx_wed(1:wed_first_index-1)) < length(new_values_first)
    xx_wed_new = [new_values_first(1:end), xx_wed(wed_first_index:wed_last_index), new_values_end(2:end)];
    % Re-define s for phase adjustment
    s_wed = 0:length(xx_wed_new)-1;
    cutoff_wed1 = length(new_values_first(1:end-1)) - length(a_wed(1:wed_first_index-1));
    cutoff_wed2 = length(new_values_end(2:end))  - length(a_wed(wed_last_index+1:end));
    a_wed_phase = [a_wed(end-cutoff_wed1:end); a_wed ; a_wed(1:cutoff_wed2)];
    [vP_wed,p_wed,V_wed,VAR_wed,CI_wed]  = csapsGCV(s_wed,a_wed_phase, [], xx_wed_new);







else
    xx_wed_new = [xx_wed(1:wed_last_index), new_values_end(2:end)];
    xx_wed_new(wed_first_index-length(new_values_first(1:end-1)):wed_first_index-1) = new_values_first(1:end-1);
    % Re-define s for phase adjustment
    s_wed = 0:length(xx_wed_new)-1;
    cutoff_wed = length(new_values_end(2:end))  - length(a_wed(wed_last_index+1:end));
    a_wed_phase = [a_wed ; a_wed(1:cutoff_wed)  ];
    [vP_wed,p_wed,V_wed,VAR_wed,CI_wed]  = csapsGCV(s_wed,a_wed_phase, [], xx_wed_new);


    index_wed = find(xx_wed_new, 1);
    xx_wed_new = xx_wed_new(index_wed-1:end);
    vP_wed = vP_wed(index_wed-1:end);


end



% Amundsen-Bellingshausen
% Re-perform GCV with extrapolated matrices (define smooth function of time)  
m_bell =  diff(xx_bell(bell_last_index-1:bell_last_index)); 
new_values_first = xx_bell(bell_first_index):-m_bell:0; 
new_values_first = [new_values_first 0]; 
new_values_first = flip(new_values_first); 
new_values_end = xx_bell(bell_last_index):m_bell:364;


if length(xx_bell(1:bell_first_index-1)) < length(new_values_first)
    xx_bell_new = [new_values_first(1:end), xx_bell(bell_first_index:bell_last_index), new_values_end(2:end)];
    % Re-define s for phase adjustment
    s_bell = 0:length(xx_bell_new)-1;
    cutoff_bell1 = length(new_values_first(1:end-1)) - length(a_bell(1:bell_first_index-1));
    cutoff_bell2 = length(new_values_end(2:end))  - length(a_bell(bell_last_index+1:end));
    a_bell_phase = [a_bell(end-cutoff_bell1:end); a_bell ; a_bell(1:cutoff_bell2)];
    [vP_bell,p_bell,V_bell,VAR_bell,CI_bell]  = csapsGCV(s_bell,a_bell_phase, [], xx_bell_new);






else
    xx_bell_new = [xx_bell(1:bell_last_index), new_values_end(2:end)];
    xx_bell_new(bell_first_index-length(new_values_first(1:end-1)):bell_first_index-1) = new_values_first(1:end-1);
    % Re-define s for phase adjustment
    s_bell = 0:length(xx_bell_new)-1;
    cutoff_bell = length(new_values_end(2:end))  - length(a_bell(bell_last_index+1:end));
    a_bell_phase = [a_bell ; a_bell(1:cutoff_bell)  ];
    [vP_bell,p_bell,V_bell,VAR_bell,CI_bell]  = csapsGCV(s_bell,a_bell_phase, [], xx_bell_new);


    index_bell = find(xx_bell_new, 1);
    xx_bell_new = xx_bell_new(index_bell-1:end);
    vP_bell = vP_bell(index_bell-1:end);


end



% Total
% Re-perform GCV with extrapolated matrices (define smooth function of time)
m_tot =  diff(xx_tot(tot_last_index-1:tot_last_index));
new_values_first = xx_tot(tot_first_index):-m_tot:0;
new_values_first = [new_values_first 0];
new_values_first = flip(new_values_first);
new_values_end = xx_tot(tot_last_index):m_tot:364;


if length(xx_tot(1:tot_first_index-1)) < length(new_values_first)
    xx_tot_new = [new_values_first(1:end), xx_tot(tot_first_index:tot_last_index), new_values_end(2:end)];
    % Re-define s for phase adjustment
    s_tot = 0:length(xx_tot_new)-1;
    cutoff_tot1 = length(new_values_first(1:end-1)) - length(a_tot(1:tot_first_index-1));
    cutoff_tot2 = length(new_values_end(2:end))  - length(a_tot(tot_last_index+1:end));
    a_tot_phase = [a_tot(end-cutoff_tot1:end); a_tot ; a_tot(1:cutoff_tot2)];
    [vP_tot,p_tot,V_tot,VAR_tot,CI_tot]  = csapsGCV(s_tot,a_tot_phase, [], xx_tot_new);







else
    xx_tot_new = [xx_tot(1:tot_last_index), new_values_end(2:end)];
    xx_tot_new(tot_first_index-length(new_values_first(1:end-1)):tot_first_index-1) = new_values_first(1:end-1);
    % Re-define s for phase adjustment
    s_tot = 0:length(xx_tot_new)-1;
    cutoff_tot = length(new_values_end(2:end))  - length(a_tot(tot_last_index+1:end));
    a_tot_phase = [a_tot ; a_tot(1:cutoff_tot)  ];
    [vP_tot,p_tot,V_tot,VAR_tot,CI_tot]  = csapsGCV(s_tot,a_tot_phase, [], xx_tot_new);


    index_tot = find(xx_tot_new, 1);
    xx_tot_new = xx_tot_new(index_tot-1:end);
    vP_tot = vP_tot(index_tot-1:end);



end



% Find the index of the minima for all regions 
[~, minIndex_ross] = min(vP_ross);
[~, minIndex_kh] = min(vP_kh);
[~, minIndex_ea] = min(vP_ea);
[~, minIndex_wed] = min(vP_wed);
[~, minIndex_bell] = min(vP_bell);
[~, minIndex_tot] = min(vP_tot);


% Reorganize a[s] so day-of-cycle = 0 corresponds to annual minimum 
vP_ross_cycle = [vP_ross(minIndex_ross:end)  vP_ross(1:minIndex_ross-1)];
vP_kh_cycle = [vP_kh(minIndex_kh:end)  vP_kh(1:minIndex_kh-1)];
vP_ea_cycle = [vP_ea(minIndex_ea:end)  vP_ea(1:minIndex_ea-1)];
vP_wed_cycle = [vP_wed(minIndex_wed:end)  vP_wed(1:minIndex_wed-1)];
vP_bell_cycle = [vP_bell(minIndex_bell:end)  vP_bell(1:minIndex_bell-1)];
vP_tot_cycle = [vP_tot(minIndex_tot:end)  vP_tot(1:minIndex_tot-1)];




% Save xx data 
% Create years array 
phase_years_ross = years(year_index)*ones(1,length(xx_ross_new)); 
phase_years_kh = years(year_index)*ones(1,length(xx_kh_new)); 
phase_years_ea = years(year_index)*ones(1,length(xx_ea_new)); 
phase_years_wed = years(year_index)*ones(1,length(xx_wed_new));
phase_years_bell = years(year_index)*ones(1,length(xx_bell_new)); 
phase_years_tot = years(year_index)*ones(1, length(xx_tot_new)); 

phase_years_ross_all = [phase_years_ross_all, phase_years_ross]; 
phase_years_kh_all = [phase_years_kh_all, phase_years_kh]; 
phase_years_ea_all = [phase_years_ea_all, phase_years_ea]; 
phase_years_wed_all = [phase_years_wed_all, phase_years_wed]; 
phase_years_bell_all = [phase_years_bell_all, phase_years_bell]; 
phase_years_tot_all = [phase_years_tot_all, phase_years_tot]; 

% Append the current array to the combined array
xx_phase_adjusted_cycle_ross_all = [xx_phase_adjusted_cycle_ross_all, xx_ross_new];
xx_phase_adjusted_cycle_kh_all = [xx_phase_adjusted_cycle_kh_all, xx_kh_new]; 
xx_phase_adjusted_cycle_ea_all = [xx_phase_adjusted_cycle_ea_all, xx_ea_new];
xx_phase_adjusted_cycle_wed_all = [xx_phase_adjusted_cycle_wed_all, xx_wed_new];
xx_phase_adjusted_cycle_bell_all = [xx_phase_adjusted_cycle_bell_all, xx_bell_new]; 
xx_phase_adjusted_cycle_tot_all = [xx_phase_adjusted_cycle_tot_all, xx_tot_new];





%{


figure(15)
plot(day_of_cycle,a_ross_cycle,'linewidth',4,'color','red', 'LineStyle', '-', 'Marker', 'none')
hold on 
plot(day_of_cycle,ross_SIE(:,year_index),'linewidth',4,'color','black', 'LineStyle', '-', 'Marker', 'none')
hold on 
plot(xx_ross_new,vP_ross_cycle,'linewidth',2,'color','green', 'LineStyle', '-', 'Marker', 'none')
xlabel('Julian Day','fontsize',25)
ylabel('Regional Sea Ice Extent [10^6 km^2]','fontsize',25)
title('Ross','fontsize',25)
xlim([0 365])
set(gcf,'color','white')
set(gca,'fontsize',20)
legend('Traditional', 'Recorded in 2016','Phase-adjusted')
grid on 


figure(16)
plot(day_of_cycle,a_kh_cycle,'linewidth',4,'color','red', 'LineStyle', '-', 'Marker', 'none')
hold on 
plot(day_of_cycle,kh_SIE(:,year_index),'linewidth',4,'color','black', 'LineStyle', '-', 'Marker', 'none')
hold on 
plot(xx_kh_new,vP_kh_cycle,'linewidth',2,'color','green', 'LineStyle', '-', 'Marker', 'none')
xlabel('Julian Day','fontsize',25)
ylabel('Regional Sea Ice Extent [10^6 km^2]','fontsize',25)
title('King Hakon','fontsize',25)
xlim([0 365])
set(gcf,'color','white')
set(gca,'fontsize',20)
legend('Traditional', 'Recorded in 2016','Phase-adjusted')
grid on 

figure(17)
plot(day_of_cycle,a_ea_cycle,'linewidth',4,'color','red', 'LineStyle', '-', 'Marker', 'none')
hold on 
plot(day_of_cycle,ea_SIE(:,year_index),'linewidth',4,'color','black', 'LineStyle', '-', 'Marker', 'none')
hold on 
plot(xx_ea_new,vP_ea_cycle,'linewidth',2,'color','green', 'LineStyle', '-', 'Marker', 'none')
xlabel('Julian Day','fontsize',25)
ylabel('Regional Sea Ice Extent [10^6 km^2]','fontsize',25)
title('East Antarctica','fontsize',25)
xlim([0 365])
set(gcf,'color','white')
set(gca,'fontsize',20)
legend('Traditional', 'Recorded in 2016','Phase-adjusted')
grid on 


figure(18)
plot(day_of_cycle,a_wed_cycle,'linewidth',4,'color','red', 'LineStyle', '-', 'Marker', 'none')
hold on 
plot(day_of_cycle,wed_SIE(:,year_index),'linewidth',4,'color','black', 'LineStyle', '-', 'Marker', 'none')
hold on 
plot(xx_wed_new,vP_wed_cycle,'linewidth',2,'color','green', 'LineStyle', '-', 'Marker', 'none')
xlabel('Julian Day','fontsize',25)
ylabel('Regional Sea Ice Extent [10^6 km^2]','fontsize',25)
title('Weddell','fontsize',25)
xlim([0 365])
set(gcf,'color','white')
set(gca,'fontsize',20)
legend('Traditional', 'Recorded in 2016','Phase-adjusted')
grid on 


figure(19)
plot(day_of_cycle,a_bell_cycle,'linewidth',4,'color','red', 'LineStyle', '-', 'Marker', 'none')
hold on 
plot(day_of_cycle,bell_SIE(:,year_index),'linewidth',4,'color','black', 'LineStyle', '-', 'Marker', 'none')
hold on 
plot(xx_bell_new,vP_bell_cycle,'linewidth',2,'color','green', 'LineStyle', '-', 'Marker', 'none')
xlabel('Julian Day','fontsize',25)
ylabel('Regional Sea Ice Extent [10^6 km^2]','fontsize',25)
title('Bellinghausen','fontsize',25)
xlim([0 365])
set(gcf,'color','white')
set(gca,'fontsize',20)
legend('Traditional', 'Recorded in 2016','Phase-adjusted')
grid on


figure(20)
plot(day_of_cycle,a_tot_cycle,'linewidth',4,'color','red', 'LineStyle', '-', 'Marker', 'none')
hold on 
plot(day_of_cycle,tot_SIE(:,year_index),'linewidth',4,'color','black', 'LineStyle', '-', 'Marker', 'none')
hold on 
plot(xx_tot_new,vP_tot_cycle,'linewidth',2,'color','green', 'LineStyle', '-', 'Marker', 'none')
xlabel('Julian Day','fontsize',25)
ylabel('Total Sea Ice Extent [10^6 km^2]','fontsize',25)
xlim([0 365])
set(gcf,'color','white')
set(gca,'fontsize',20)
legend('Traditional', 'Recorded in 2016','Phase-adjusted')
grid on 
%}

%% Phase and Amplitude Adjusted 

% Define variables
s = day_of_cycle(:); 

% Scale inputs by annual maxima and minimum
for j = 1:length(years)
    input_ross(:,j) = (rossreshaped(:,j)-min(rossreshaped(:,j)))/(max(rossreshaped(:,j))-min(rossreshaped(:,j)));
    input_kh(:,j) = (khreshaped(:,j)-min(khreshaped(:,j)))/(max(khreshaped(:,j))-min(khreshaped(:,j))); 
    input_ea(:,j) = (eareshaped(:,j)-min(eareshaped(:,j)))/(max(eareshaped(:,j))-min(rossreshaped(:,j))); 
    input_wed(:,j) = (wedreshaped(:,j)-min(wedreshaped(:,j)))/(max(wedreshaped(:,j))-min(wedreshaped(:,j))); 
    input_bell(:,j) = (bellreshaped(:,j)-min(bellreshaped(:,j)))/(max(bellreshaped(:,j))-min(bellreshaped(:,j))); 
    input_tot(:,j) = (totreshaped(:,j)-min(totreshaped(:,j)))/(max(totreshaped(:,j))-min(totreshaped(:,j))); 
end

%{
input_ross = [rossreshaped(:,1) input_ross rossreshaped(:,47)]; 
input_kh = [khreshaped(:,1) input_kh khreshaped(:,47)]; 
input_ea = [eareshaped(:,1) input_ea eareshaped(:,47)]; 
input_wed = [wedreshaped(:,1) input_wed wedreshaped(:,47)]; 
input_bell = [bellreshaped(:,1) input_bell bellreshaped(:,47)]; 
input_tot = [totreshaped(:,1) input_tot totreshaped(:,47)]; 
%}

% Initialize a mask for NaN values
nanMask_ross = isnan(input_ross);
nanMask_kh = isnan(input_kh);
nanMask_ea = isnan(input_ea);
nanMask_wed = isnan(input_wed);
nanMask_bell = isnan(input_bell);
nanMask_tot = isnan(input_tot);


% Compute the sum along rows ignoring NaN values
rowSums_ross = sum(input_ross, 2, 'omitnan');
rowSums_kh = sum(input_kh, 2, 'omitnan');
rowSums_ea = sum(input_ea, 2, 'omitnan');
rowSums_wed = sum(input_wed, 2, 'omitnan');
rowSums_bell = sum(input_bell, 2, 'omitnan');
rowSums_tot = sum(input_tot, 2, 'omitnan');

% Count non-NaN elements along rows
counts_ross = sum(~nanMask_ross, 2);
counts_kh = sum(~nanMask_kh, 2);
counts_ea = sum(~nanMask_ea, 2);
counts_wed = sum(~nanMask_wed, 2);
counts_bell = sum(~nanMask_bell, 2);
counts_tot = sum(~nanMask_tot, 2);

% Calculate the mean row-wise
input_ross = rowSums_ross ./ counts_ross;
input_kh = rowSums_kh ./ counts_kh;
input_ea = rowSums_ea ./ counts_ea;
input_wed = rowSums_wed ./ counts_wed;
input_bell = rowSums_bell ./ counts_bell;
input_tot = rowSums_tot ./ counts_tot;


% Define phase(t) (AKA xx) for 2020 
xx_ross = zeros(1,365); 
xx_kh = zeros(1,365); 
xx_ea = zeros(1,365); 
xx_bell = zeros(1,365); 
xx_wed = zeros(1,365); 
xx_tot = zeros(1,365); 

% Initialize variables to store the first and last indices
ross_first_index = [];
ross_last_index = [];
kh_first_index = [];
kh_last_index = [];
ea_first_index = [];
ea_last_index = [];
wed_first_index = [];
wed_last_index = [];
bell_first_index = [];
bell_last_index = [];
tot_first_index = [];
tot_last_index = [];

year = year_index;  

for i = 1:365
    if timereshaped(i,year) >= min_extent_day_ross(year) && timereshaped(i,year) <= max_extent_day_ross(year)

        xx_ross(i) = (timereshaped(i,year) - min_extent_day_ross(year))/(max_extent_day_ross(year+1)- min_extent_day_ross(year));
  
        % Save the first index where the condition is met
        if isempty(ross_first_index)
            ross_first_index = i;
       
        end
        
        % Update the last index where the condition is met
        ross_last_index = i;
   
    else
        xx_ross(i) = 0;



    end
  
end


for i = 1:365
     if timereshaped(i,year) >= min_extent_day_kh(year) && timereshaped(i,year) <= max_extent_day_kh(year)

    
        xx_kh(i) = (timereshaped(i,year) - min_extent_day_kh(year))/(max_extent_day_kh(year+1)- min_extent_day_kh(year));
       
        % Save the first index where the condition is met
        if isempty(kh_first_index)
         
            kh_first_index = i;
       
        end
        
        % Update the last index where the condition is met
  
        kh_last_index = i;
 
    else
 
        xx_kh(i) = 0;
    


    end
       
end



for i = 1:365
     if timereshaped(i,year) >= min_extent_day_ea(year) && timereshaped(i,year) <= max_extent_day_ea(year)

    
        xx_ea(i) = (timereshaped(i,year) - min_extent_day_ea(year))/(max_extent_day_ea(year+1)- min_extent_day_ea(year));
      
        % Save the first index where the condition is met
        if isempty(ea_first_index)
         
            ea_first_index = i;
       
        end
        
        % Update the last index where the condition is met
  
        ea_last_index = i;
 
    else
 
        xx_ea(i) = 0;
    


    end
      
end


for i = 1:365
     if timereshaped(i,year) >= min_extent_day_wed(year) && timereshaped(i,year) <= max_extent_day_wed(year)

    
        xx_wed(i) = (timereshaped(i,year) - min_extent_day_wed(year))/(max_extent_day_wed(year+1)- min_extent_day_wed(year));
       
        % Save the first index where the condition is met
        if isempty(wed_first_index)
         
            wed_first_index = i;
       
        end
        
        % Update the last index where the condition is met
  
        wed_last_index = i;
 
    else
 
        xx_wed(i) = 0;
    


    end
       
end



for i = 1:365
     if timereshaped(i,year) >= min_extent_day_bell(year) && timereshaped(i,year) <= max_extent_day_bell(year)

    
        xx_bell(i) = (timereshaped(i,year) - min_extent_day_bell(year))/(max_extent_day_bell(year+1)- min_extent_day_bell(year));
       
        % Save the first index where the condition is met
        if isempty(bell_first_index)
         
            bell_first_index = i;
       
        end
        
        % Update the last index where the condition is met
  
        bell_last_index = i;
 
    else
 
        xx_bell(i) = 0;
    


    end
       
end



for i = 1:365
     if timereshaped(i,year) >= min_extent_day_tot(year) && timereshaped(i,year) <= max_extent_day_tot(year)

    
        xx_tot(i) = (timereshaped(i,year) - min_extent_day_tot(year))/(max_extent_day_tot(year+1)- min_extent_day_tot(year));
       
        
        % Save the first index where the condition is met
        if isempty(tot_first_index)
         
            tot_first_index = i;
       
        end
        
        % Update the last index where the condition is met
  
        tot_last_index = i;
 
    else
 
        xx_tot(i) = 0;
    


    end
        
end






[vPa_ross,p_ross,V_ross,VAR_ross,CI_ross,beta1_ross,beta2_ross,xx_ross]  = csapsGCV_with_Phase(s,input_ross, [], xx_ross); 
[vPa_kh,p_kh,V_kh,VAR_kh,CI_kh,beta1_kh,beta2_kh,xx_kh]  = csapsGCV_with_Phase(s,input_kh, [], xx_kh);
[vPa_ea,p_ea,V_ea,VAR_ea,CI_ea,beta1_ea,beta2_ea,xx_ea]  = csapsGCV_with_Phase(s,input_ea, [], xx_ea);
[vPa_wed,p_wed,V_wed,VAR_wed,CI_wed,beta1_wed,beta2_wed,xx_wed]  = csapsGCV_with_Phase(s,input_wed, [], xx_wed);
[vPa_bell,p_bell,V_bell,VAR_bell,CI_bell,beta1_bell,beta2_bell,xx_bell]  = csapsGCV_with_Phase(s,input_bell, [], xx_bell);
[vPa_tot,p_tot,V_tot,VAR_tot,CI_tot,beta1_tot,beta2_tot,xx_tot]  = csapsGCV_with_Phase(s,input_tot, [], xx_tot);




% Ross
% Re-perform GCV with extrapolated matrices (define smooth function of time)  
m_ross =  diff(xx_ross(ross_last_index-1:ross_last_index)); 
new_values_first = xx_ross(ross_first_index):-m_ross:0; 
new_values_first = [new_values_first 0]; 
new_values_first = flip(new_values_first); 
new_values_end = xx_ross(ross_last_index):m_ross:364;


if length(xx_ross(1:ross_first_index-1)) < length(new_values_first)
    xx_ross_new = [new_values_first(1:end), xx_ross(ross_first_index:ross_last_index), new_values_end(2:end)];
    % Re-define s for phase adjustment
    s_ross = 0:length(xx_ross_new)-1;
    cutoff_ross1 = length(new_values_first(1:end-1)) - length(input_ross(1:ross_first_index-1));
    cutoff_ross2 = length(new_values_end(2:end))  - length(input_ross(ross_last_index+1:end));
    input_ross_phase = [input_ross(end-cutoff_ross1:end); input_ross ; input_ross(1:cutoff_ross2)];
    [vPa_ross,p_ross,V_ross,VAR_ross,CI_ross]  = csapsGCV(s_ross,input_ross_phase, [], xx_ross_new);






else
    xx_ross_new = [xx_ross(1:ross_last_index), new_values_end(2:end)];
    xx_ross_new(ross_first_index-length(new_values_first(1:end-1)):ross_first_index-1) = new_values_first(1:end-1);
    % Re-define s for phase adjustment
    s_ross = 0:length(xx_ross_new)-1;
    cutoff_ross = length(new_values_end(2:end))  - length(input_ross(ross_last_index+1:end));
    input_ross_phase = [input_ross ; input_ross(1:cutoff_ross)  ];
    [vPa_ross,p_ross,V_ross,VAR_ross,CI_ross]  = csapsGCV(s_ross,input_ross_phase, [], xx_ross_new);


    index_ross = find(xx_ross_new, 1);
    xx_ross_new = xx_ross_new(index_ross-1:end);
    vPa_ross = vPa_ross(index_ross-1:end);


end








% King Hakon
% Re-perform GCV with extrapolated matrices (define smooth function of time)  
m_kh =  diff(xx_kh(kh_last_index-1:kh_last_index)); 
new_values_first = xx_kh(kh_first_index):-m_kh:0; 
new_values_first = [new_values_first 0]; 
new_values_first = flip(new_values_first); 
new_values_end = xx_kh(kh_last_index):m_kh:364;


if length(xx_kh(1:kh_first_index-1)) < length(new_values_first)
    xx_kh_new = [new_values_first(1:end), xx_kh(kh_first_index:kh_last_index), new_values_end(2:end)];
    % Re-define s for phase adjustment
    s_kh = 0:length(xx_kh_new)-1;
    cutoff_kh1 = length(new_values_first(1:end-1)) - length(input_kh(1:kh_first_index-1));
    cutoff_kh2 = length(new_values_end(2:end))  - length(input_kh(kh_last_index+1:end));
    input_kh_phase = [input_kh(end-cutoff_kh1:end); input_kh ; input_kh(1:cutoff_kh2)];
    [vPa_kh,p_kh,V_kh,VAR_kh,CI_kh]  = csapsGCV(s_kh,input_kh_phase, [], xx_kh_new);






else
    xx_kh_new = [xx_kh(1:kh_last_index), new_values_end(2:end)];
    xx_kh_new(kh_first_index-length(new_values_first(1:end-1)):kh_first_index-1) = new_values_first(1:end-1);
    % Re-define s for phase adjustment
    s_kh = 0:length(xx_kh_new)-1;
    cutoff_kh = length(new_values_end(2:end))  - length(input_kh(kh_last_index+1:end));
    input_kh_phase = [input_kh ; input_kh(1:cutoff_kh)  ];
    [vPa_kh,p_kh,V_kh,VAR_kh,CI_kh]  = csapsGCV(s_kh,input_kh_phase, [], xx_kh_new);


    index_kh = find(xx_kh_new, 1);
    xx_kh_new = xx_kh_new(index_kh-1:end);
    vPa_kh = vPa_kh(index_kh-1:end);


end



% East Antarctica
% Re-perform GCV with extrapolated matrices (define smooth function of time)  
m_ea =  diff(xx_ea(ea_last_index-1:ea_last_index)); 
new_values_first = xx_ea(ea_first_index):-m_ea:0; 
new_values_first = [new_values_first 0]; 
new_values_first = flip(new_values_first); 
new_values_end = xx_ea(ea_last_index):m_ea:364;


if length(xx_ea(1:ea_first_index-1)) < length(new_values_first)
    xx_ea_new = [new_values_first(1:end), xx_ea(ea_first_index:ea_last_index), new_values_end(2:end)];
    % Re-define s for phase adjustment
    s_ea = 0:length(xx_ea_new)-1;
    cutoff_ea1 = length(new_values_first(1:end-1)) - length(input_ea(1:ea_first_index-1));
    cutoff_ea2 = length(new_values_end(2:end))  - length(input_ea(ea_last_index+1:end));
    input_ea_phase = [input_ea(end-cutoff_ea1:end); input_ea ; input_ea(1:cutoff_ea2)];
    [vPa_ea,p_ea,V_ea,VAR_ea,CI_ea]  = csapsGCV(s_ea,input_ea_phase, [], xx_ea_new);






else
    xx_ea_new = [xx_ea(1:ea_last_index), new_values_end(2:end)];
    xx_ea_new(ea_first_index-length(new_values_first(1:end-1)):ea_first_index-1) = new_values_first(1:end-1);
    % Re-define s for phase adjustment
    s_ea = 0:length(xx_ea_new)-1;
    cutoff_ea = length(new_values_end(2:end))  - length(input_ea(ea_last_index+1:end));
    input_ea_phase = [input_ea ; input_ea(1:cutoff_ea)  ];
    [vPa_ea,p_ea,V_ea,VAR_ea,CI_ea]  = csapsGCV(s_ea,input_ea_phase, [], xx_ea_new);

    index_ea = find(xx_ea_new, 1);
    xx_ea_new = xx_ea_new(index_ea-1:end);
    vPa_ea = vPa_ea(index_ea-1:end);


end


% Weddell
% Re-perform GCV with extrapolated matrices (define smooth function of time)  
m_wed =  diff(xx_wed(wed_last_index-1:wed_last_index)); 
new_values_first = xx_wed(wed_first_index):-m_wed:0; 
new_values_first = [new_values_first 0]; 
new_values_first = flip(new_values_first); 
new_values_end = xx_wed(wed_last_index):m_wed:364;


if length(xx_wed(1:wed_first_index-1)) < length(new_values_first)
    xx_wed_new = [new_values_first(1:end), xx_wed(wed_first_index:wed_last_index), new_values_end(2:end)];
    % Re-define s for phase adjustment
    s_wed = 0:length(xx_wed_new)-1;
    cutoff_wed1 = length(new_values_first(1:end-1)) - length(input_wed(1:wed_first_index-1));
    cutoff_wed2 = length(new_values_end(2:end))  - length(input_wed(wed_last_index+1:end));
    input_wed_phase = [input_wed(end-cutoff_wed1:end); input_wed ; input_wed(1:cutoff_wed2)];
    [vPa_wed,p_wed,V_wed,VAR_wed,CI_wed]  = csapsGCV(s_wed,input_wed_phase, [], xx_wed_new);







else
    xx_wed_new = [xx_wed(1:wed_last_index), new_values_end(2:end)];
    xx_wed_new(wed_first_index-length(new_values_first(1:end-1)):wed_first_index-1) = new_values_first(1:end-1);
    % Re-define s for phase adjustment
    s_wed = 0:length(xx_wed_new)-1;
    cutoff_wed = length(new_values_end(2:end))  - length(input_wed(wed_last_index+1:end));
    input_wed_phase = [input_wed ; input_wed(1:cutoff_wed)  ];
    [vPa_wed,p_wed,V_wed,VAR_wed,CI_wed]  = csapsGCV(s_wed,input_wed_phase, [], xx_wed_new);


    index_wed = find(xx_wed_new, 1);
    xx_wed_new = xx_wed_new(index_wed-1:end);
    vPa_wed = vPa_wed(index_wed-1:end);


end



% Amundsen-Bellingshausen
% Re-perform GCV with extrapolated matrices (define smooth function of time)  
m_bell =  diff(xx_bell(bell_last_index-1:bell_last_index)); 
new_values_first = xx_bell(bell_first_index):-m_bell:0; 
new_values_first = [new_values_first 0]; 
new_values_first = flip(new_values_first); 
new_values_end = xx_bell(bell_last_index):m_bell:364;


if length(xx_bell(1:bell_first_index-1)) < length(new_values_first)
    xx_bell_new = [new_values_first(1:end), xx_bell(bell_first_index:bell_last_index), new_values_end(2:end)];
    % Re-define s for phase adjustment
    s_bell = 0:length(xx_bell_new)-1;
    cutoff_bell1 = length(new_values_first(1:end-1)) - length(input_bell(1:bell_first_index-1));
    cutoff_bell2 = length(new_values_end(2:end))  - length(input_bell(bell_last_index+1:end));
    input_bell_phase = [input_bell(end-cutoff_bell1:end); input_bell ; input_bell(1:cutoff_bell2)];
    [vPa_bell,p_bell,V_bell,VAR_bell,CI_bell]  = csapsGCV(s_bell,input_bell_phase, [], xx_bell_new);






else
    xx_bell_new = [xx_bell(1:bell_last_index), new_values_end(2:end)];
    xx_bell_new(bell_first_index-length(new_values_first(1:end-1)):bell_first_index-1) = new_values_first(1:end-1);
    % Re-define s for phase adjustment
    s_bell = 0:length(xx_bell_new)-1;
    cutoff_bell = length(new_values_end(2:end))  - length(input_bell(bell_last_index+1:end));
    input_bell_phase = [input_bell ; input_bell(1:cutoff_bell)  ];
    [vPa_bell,p_bell,V_bell,VAR_bell,CI_bell]  = csapsGCV(s_bell,input_bell_phase, [], xx_bell_new);


    index_bell = find(xx_bell_new, 1);
    xx_bell_new = xx_bell_new(index_bell-1:end);
    vPa_bell = vPa_bell(index_bell-1:end);


end



% Total
% Re-perform GCV with extrapolated matrices (define smooth function of time)
m_tot =  diff(xx_tot(tot_last_index-1:tot_last_index));
new_values_first = xx_tot(tot_first_index):-m_tot:0;
new_values_first = [new_values_first 0];
new_values_first = flip(new_values_first);
new_values_end = xx_tot(tot_last_index):m_tot:364;


if length(xx_tot(1:tot_first_index-1)) < length(new_values_first)
    xx_tot_new = [new_values_first(1:end), xx_tot(tot_first_index:tot_last_index), new_values_end(2:end)];
    % Re-define s for phase adjustment
    s_tot = 0:length(xx_tot_new)-1;
    cutoff_tot1 = length(new_values_first(1:end-1)) - length(input_tot(1:tot_first_index-1));
    cutoff_tot2 = length(new_values_end(2:end))  - length(input_tot(tot_last_index+1:end));
    input_tot_phase = [input_tot(end-cutoff_tot1:end); input_tot ; input_tot(1:cutoff_tot2)];
    [vPa_tot,p_tot,V_tot,VAR_tot,CI_tot]  = csapsGCV(s_tot,input_tot_phase, [], xx_tot_new);







else
    xx_tot_new = [xx_tot(1:tot_last_index), new_values_end(2:end)];
    xx_tot_new(tot_first_index-length(new_values_first(1:end-1)):tot_first_index-1) = new_values_first(1:end-1);
    % Re-define s for phase adjustment
    s_tot = 0:length(xx_tot_new)-1;
    cutoff_tot = length(new_values_end(2:end))  - length(input_tot(tot_last_index+1:end));
    input_tot_phase = [input_tot ; input_tot(1:cutoff_tot)  ];
    [vPa_tot,p_tot,V_tot,VAR_tot,CI_tot]  = csapsGCV(s_tot,input_tot_phase, [], xx_tot_new);


    index_tot = find(xx_tot_new, 1);
    xx_tot_new = xx_tot_new(index_tot-1:end);
    vPa_tot = vPa_tot(index_tot-1:end);



end





aAp_ross = vPa_ross*(max(rossreshaped(:,year_index))-min(rossreshaped(:,year_index)))+min(rossreshaped(:,year_index)); %2019
aAp_kh = vPa_kh*(max(khreshaped(:,year_index))-min(khreshaped(:,year_index)))+min(khreshaped(:,year_index)); 
aAp_ea = vPa_ea*(max(eareshaped(:,year_index))-min(eareshaped(:,year_index)))+min(eareshaped(:,year_index)); 
aAp_wed = vPa_wed*(max(wedreshaped(:,year_index))-min(wedreshaped(:,year_index)))+min(wedreshaped(:,year_index)); 
aAp_bell = vPa_bell*(max(bellreshaped(:,year_index))-min(bellreshaped(:,year_index)))+min(bellreshaped(:,year_index)); 
aAp_tot = vPa_tot*(max(totreshaped(:,year_index))-min(totreshaped(:,year_index)))+min(totreshaped(:,year_index)); 







% Find the index of the minima for all regions 
[~, minIndex_ross] = min(aAp_ross);
[~, minIndex_kh] = min(aAp_kh);
[~, minIndex_ea] = min(aAp_ea);
[~, minIndex_wed] = min(aAp_wed);
[~, minIndex_bell] = min(aAp_bell);
[~, minIndex_tot] = min(aAp_tot);


% Reorganize a[s] so day-of-cycle = 0 corresponds to annual minimum 
vPa_ross_cycle = [aAp_ross(minIndex_ross:end)  aAp_ross(1:minIndex_ross-1)];
vPa_kh_cycle = [aAp_kh(minIndex_kh:end)  aAp_kh(1:minIndex_kh-1)];
vPa_ea_cycle = [aAp_ea(minIndex_ea:end)  aAp_ea(1:minIndex_ea-1)];
vPa_wed_cycle = [aAp_wed(minIndex_wed:end)  aAp_wed(1:minIndex_wed-1)];
vPa_bell_cycle = [aAp_bell(minIndex_bell:end)  aAp_bell(1:minIndex_bell-1)];
vPa_tot_cycle = [aAp_tot(minIndex_tot:end)  aAp_tot(1:minIndex_tot-1)];



% Save xx data 
% Create years array 
APAC_years_ross = years(year_index)*ones(1,length(xx_ross_new)); 
APAC_years_kh = years(year_index)*ones(1,length(xx_kh_new)); 
APAC_years_ea = years(year_index)*ones(1,length(xx_ea_new)); 
APAC_years_wed = years(year_index)*ones(1,length(xx_wed_new));
APAC_years_bell = years(year_index)*ones(1,length(xx_bell_new)); 
APAC_years_tot = years(year_index)*ones(1,length(xx_tot_new)); 

APAC_years_ross_all = [APAC_years_ross_all, APAC_years_ross]; 
APAC_years_kh_all = [APAC_years_kh_all, APAC_years_kh]; 
APAC_years_ea_all = [APAC_years_ea_all, APAC_years_ea]; 
APAC_years_wed_all = [APAC_years_wed_all, APAC_years_wed]; 
APAC_years_bell_all = [APAC_years_bell_all, APAC_years_bell]; 
APAC_years_tot_all = [APAC_years_tot_all, APAC_years_tot]; 


% Append the current array to the combined array
xx_APAC_adjusted_cycle_ross_all = [xx_APAC_adjusted_cycle_ross_all, xx_ross_new];
xx_APAC_adjusted_cycle_kh_all = [xx_APAC_adjusted_cycle_kh_all, xx_kh_new]; 
xx_APAC_adjusted_cycle_ea_all = [xx_APAC_adjusted_cycle_ea_all, xx_ea_new];
xx_APAC_adjusted_cycle_wed_all = [xx_APAC_adjusted_cycle_wed_all, xx_wed_new];
xx_APAC_adjusted_cycle_bell_all = [xx_APAC_adjusted_cycle_bell_all, xx_bell_new]; 
xx_APAC_adjusted_cycle_tot_all = [xx_APAC_adjusted_cycle_tot_all, xx_tot_new];


%{
figure;
subplot(2,3,1)
plot(day_of_cycle,a_ross_cycle,'linewidth',4,'color','black', 'LineStyle', ':', 'Marker', 'none')
hold on 
plot(day_of_cycle,ross_SIE(:,year_index),'linewidth',4,'color','blue', 'LineStyle', ':', 'Marker', 'none')
hold on 
plot(xx_ross_new,vPa_ross_cycle,'linewidth',2,'color','red', 'LineStyle', '-', 'Marker', 'none')
%xlabel('Day of Cycle','fontsize',25)
set(gca,'XTickLabel',[])
ylabel('Regional Sea Ice Extent [$10^{6}$ $\mathrm{km}^{2}$]','fontsize',25,'interpreter', 'latex')
title('Ross','fontsize',25,'interpreter','latex')
xlim([0 364])
set(gcf,'color','white')
ylim([0 6])
set(gca,'fontsize',20)
grid on 

subplot(2,3,2)
plot(day_of_cycle,a_bell_cycle,'linewidth',4,'color','black', 'LineStyle', ':', 'Marker', 'none')
hold on 
plot(day_of_cycle,bell_SIE(:,year_index),'linewidth',4,'color','blue', 'LineStyle', ':', 'Marker', 'none')
hold on 
plot(xx_bell_new,vPa_bell_cycle,'linewidth',2,'color','red', 'LineStyle', '-', 'Marker', 'none')
%xlabel('Day of Cycle','fontsize',25)
set(gca,'XTickLabel',[])
%ylabel('Regional Sea Ice Extent [10^6 km^2]','fontsize',25)
set(gca,'YTickLabel',[])
title('Amundsen-Bellingshausen','fontsize',25, 'interpreter','latex')
xlim([0 364])
set(gcf,'color','white')
ylim([0 6])
set(gca,'fontsize',20)
grid on

subplot(2,3,3)
plot(day_of_cycle,a_wed_cycle,'linewidth',4,'color','black', 'LineStyle', ':', 'Marker', 'none')
hold on 
plot(day_of_cycle,wed_SIE(:,year_index),'linewidth',4,'color','blue', 'LineStyle', ':', 'Marker', 'none')
hold on 
plot(xx_wed_new,vPa_wed_cycle,'linewidth',2,'color','red', 'LineStyle', '-', 'Marker', 'none')
%xlabel('Day of Cycle','fontsize',25)
set(gca,'XTickLabel',[])
%ylabel('Regional Sea Ice Extent [10^6 km^2]','fontsize',25)
set(gca,'YTickLabel',[])
title('Weddell','fontsize',25, 'interpreter','latex')
xlim([0 364])
set(gcf,'color','white')
ylim([0 6])
set(gca,'fontsize',20)
grid on 

subplot(2,3,4)
plot(day_of_cycle,a_kh_cycle,'linewidth',4,'color','black', 'LineStyle', ':', 'Marker', 'none')
hold on 
plot(day_of_cycle,kh_SIE(:,year_index),'linewidth',4,'color','blue', 'LineStyle', ':', 'Marker', 'none')
hold on 
plot(xx_kh_new,vPa_kh_cycle,'linewidth',2,'color','red', 'LineStyle', '-', 'Marker', 'none')
xlabel('Day of Cycle','fontsize',25, 'interpreter','latex')
ylabel('Regional Sea Ice Extent [$10^{6}$ $\mathrm{km}^{2}$]','fontsize',25, 'interpreter','latex')
title('King Hakon','fontsize',25, 'interpreter','latex')
xlim([0 364])
set(gcf,'color','white')
ylim([0 6])
set(gca,'fontsize',20)
grid on 

subplot(2,3,5)
plot(day_of_cycle,a_ea_cycle,'linewidth',4,'color','black', 'LineStyle', ':', 'Marker', 'none')
hold on 
plot(day_of_cycle,ea_SIE(:,year_index),'linewidth',4,'color','blue', 'LineStyle', ':', 'Marker', 'none')
hold on 
plot(xx_ea_new,vPa_ea_cycle,'linewidth',2,'color','red', 'LineStyle', '-', 'Marker', 'none')
xlabel('Day of Cycle','fontsize',25,'interpreter','latex')
%ylabel('Regional Sea Ice Extent [10^6 km^2]','fontsize',25)
set(gca,'YTickLabel',[])
title('East Antarctica','fontsize',25,'interpreter','latex')
xlim([0 364])
set(gcf,'color','white')
ylim([0 6])
set(gca,'fontsize',20)
grid on 

subplot(2,3,6)
plot(day_of_cycle,a_tot_cycle,'linewidth',4,'color','black', 'LineStyle', ':', 'Marker', 'none')
hold on 
plot(day_of_cycle,tot_SIE(:,year_index),'linewidth',4,'color','blue', 'LineStyle', ':', 'Marker', 'none')
hold on 
plot(xx_tot_new,vPa_tot_cycle,'linewidth',2,'color','red', 'LineStyle', '-', 'Marker', 'none')
xlabel('Day of Cycle','fontsize',25, 'interpreter','latex')
ylabel('Total Sea Ice Extent [$10^{6}$ $\mathrm{km}^{2}$]','fontsize',25,'interpreter','latex')
title('Total','fontsize',25,'interpreter','latex')
xlim([0 364])
set(gcf,'color','white')
set(gca,'fontsize',20)
legend('Traditional', ['Recorded in ', year_string],'Amp-Phase Adjusted (APAC)', 'interpreter','latex')
grid on 

%}

%{

 % Check if the current iteration is in the list of saveIterations
    if ismember(i, saveIterations)
        % Construct the filename
        filename = sprintf('APAC_iteration_%d.fig', i);
        % Save the figure
        savefig(filename);
    end
    
    % Close the figure to avoid having too many open figures
    close;


%}






%{
figure(21)
plot(day_of_cycle,a_ross_cycle,'linewidth',4,'color','red', 'LineStyle', '-', 'Marker', 'none')
hold on 
plot(day_of_cycle,ross_SIE(:,year_index),'linewidth',4,'color','black', 'LineStyle', '-', 'Marker', 'none')
hold on 
plot(xx_ross_new,vPa_ross_cycle,'linewidth',2,'color','green', 'LineStyle', '-', 'Marker', 'none')
xlabel('Julian Day','fontsize',25)
ylabel('Regional Sea Ice Extent [10^6 km^2]','fontsize',25)
title('Ross','fontsize',25)
xlim([0 365])
set(gcf,'color','white')
set(gca,'fontsize',20)
legend('Traditional', 'Recorded in 2016','Phase and Amp-adjusted')
grid on 


figure(22)
plot(day_of_cycle,a_kh_cycle,'linewidth',4,'color','red', 'LineStyle', '-', 'Marker', 'none')
hold on 
plot(day_of_cycle,kh_SIE(:,year_index),'linewidth',4,'color','black', 'LineStyle', '-', 'Marker', 'none')
hold on 
plot(xx_kh_new,vPa_kh_cycle,'linewidth',2,'color','green', 'LineStyle', '-', 'Marker', 'none')
xlabel('Julian Day','fontsize',25)
ylabel('Regional Sea Ice Extent [10^6 km^2]','fontsize',25)
title('King Hakon','fontsize',25)
xlim([0 365])
set(gcf,'color','white')
set(gca,'fontsize',20)
legend('Traditional', 'Recorded in 2016','Phase and Amp-adjusted')
grid on 

figure(23)
plot(day_of_cycle,a_ea_cycle,'linewidth',4,'color','red', 'LineStyle', '-', 'Marker', 'none')
hold on 
plot(day_of_cycle,ea_SIE(:,year_index),'linewidth',4,'color','black', 'LineStyle', '-', 'Marker', 'none')
hold on 
plot(xx_ea_new,vPa_ea_cycle,'linewidth',2,'color','green', 'LineStyle', '-', 'Marker', 'none')
xlabel('Julian Day','fontsize',25)
ylabel('Regional Sea Ice Extent [10^6 km^2]','fontsize',25)
title('East Antarctica','fontsize',25)
xlim([0 365])
set(gcf,'color','white')
set(gca,'fontsize',20)
legend('Traditional', 'Recorded in 2016','Phase and Amp-adjusted')
grid on 


figure(24)
plot(day_of_cycle,a_wed_cycle,'linewidth',4,'color','red', 'LineStyle', '-', 'Marker', 'none')
hold on 
plot(day_of_cycle,wed_SIE(:,year_index),'linewidth',4,'color','black', 'LineStyle', '-', 'Marker', 'none')
hold on 
plot(xx_wed_new,vPa_wed_cycle,'linewidth',2,'color','green', 'LineStyle', '-', 'Marker', 'none')
xlabel('Julian Day','fontsize',25)
ylabel('Regional Sea Ice Extent [10^6 km^2]','fontsize',25)
title('Weddell','fontsize',25)
xlim([0 365])
set(gcf,'color','white')
set(gca,'fontsize',20)
legend('Traditional', 'Recorded in 2016','Phase and Amp-adjusted')
grid on 


figure(25)
plot(day_of_cycle,a_bell_cycle,'linewidth',4,'color','red', 'LineStyle', '-', 'Marker', 'none')
hold on 
plot(day_of_cycle,bell_SIE(:,year_index),'linewidth',4,'color','black', 'LineStyle', '-', 'Marker', 'none')
hold on 
plot(xx_bell_new,vPa_bell_cycle,'linewidth',2,'color','green', 'LineStyle', '-', 'Marker', 'none')
xlabel('Julian Day','fontsize',25)
ylabel('Regional Sea Ice Extent [10^6 km^2]','fontsize',25)
title('Bellinghausen','fontsize',25)
xlim([0 365])
set(gcf,'color','white')
set(gca,'fontsize',20)
legend('Traditional', 'Recorded in 2016','Phase and Amp-adjusted')
grid on


figure(26)
plot(day_of_cycle,a_tot_cycle,'linewidth',4,'color','red', 'LineStyle', '-', 'Marker', 'none')
hold on 
plot(day_of_cycle,tot_SIE(:,year_index),'linewidth',4,'color','black', 'LineStyle', '-', 'Marker', 'none')
hold on 
plot(xx_tot_new,vPa_tot_cycle,'linewidth',2,'color','green', 'LineStyle', '-', 'Marker', 'none')
xlabel('Julian Day','fontsize',25)
ylabel('Total Sea Ice Extent [10^6 km^2]','fontsize',25)
xlim([0 365])
set(gcf,'color','white')
set(gca,'fontsize',20)
legend('Traditional', 'Recorded in 2016','Phase and Amp-adjusted')
grid on 

%}

%{

%% Compute root mean square error (RMSE) 
RMSE_trad_ross = rmse(a_ross_cycle,ross_SIE(year_index));
RMSE_inv_ross = rmse(vI_ross_cycle,ross_SIE(year_index));
RMSE_amp_ross = rmse(aA_ross_cycle,ross_SIE(year_index));
RMSE_phase_ross = rmse(vP_ross_cycle,ross_SIE(year_index));
RMSE_apac_ross = rmse(vPa_ross_cycle,ross_SIE(year_index));


RMSE_trad_kh = rmse(a_kh_cycle,kh_SIE(year_index));
RMSE_inv_kh = rmse(vI_kh_cycle,kh_SIE(year_index));
RMSE_amp_kh = rmse(aA_kh_cycle,kh_SIE(year_index));
RMSE_phase_kh = rmse(vP_kh_cycle,kh_SIE(year_index));
RMSE_apac_kh = rmse(vPa_kh_cycle,kh_SIE(year_index));


RMSE_trad_ea = rmse(a_ea_cycle,ea_SIE(year_index));
RMSE_inv_ea = rmse(vI_ea_cycle,ea_SIE(year_index));
RMSE_amp_ea = rmse(aA_ea_cycle,ea_SIE(year_index));
RMSE_phase_ea = rmse(vP_ea_cycle,ea_SIE(year_index));
RMSE_apac_ea = rmse(vPa_ea_cycle,ea_SIE(year_index));


RMSE_trad_wed = rmse(a_wed_cycle,wed_SIE(year_index));
RMSE_inv_wed = rmse(vI_wed_cycle,wed_SIE(year_index));
RMSE_amp_wed = rmse(aA_wed_cycle,wed_SIE(year_index));
RMSE_phase_wed = rmse(vP_wed_cycle,wed_SIE(year_index));
RMSE_apac_wed = rmse(vPa_wed_cycle,wed_SIE(year_index));


RMSE_trad_bell = rmse(a_bell_cycle,bell_SIE(year_index));
RMSE_inv_bell = rmse(vI_bell_cycle,bell_SIE(year_index));
RMSE_amp_bell = rmse(aA_bell_cycle,bell_SIE(year_index));
RMSE_phase_bell = rmse(vP_bell_cycle,bell_SIE(year_index));
RMSE_apac_bell = rmse(vPa_bell_cycle,bell_SIE(year_index));


RMSE_trad_tot = rmse(a_tot_cycle,tot_SIE(year_index));
RMSE_inv_tot = rmse(vI_tot_cycle,tot_SIE(year_index));
RMSE_amp_tot = rmse(aA_tot_cycle,tot_SIE(year_index));
RMSE_phase_tot = rmse(vP_tot_cycle,tot_SIE(year_index));
RMSE_apac_tot = rmse(vPa_tot_cycle,tot_SIE(year_index));



%}

%{


%% Calculate the daily rate of change (first derivative of our models) 
ROC_trad_ross = diff(a_ross_cycle); 
ROC_trad_kh = diff(a_kh_cycle); 
ROC_trad_ea = diff(a_ea_cycle); 
ROC_trad_wed = diff(a_wed_cycle); 
ROC_trad_bell = diff(a_bell_cycle); 
ROC_trad_tot = diff(a_tot_cycle); 

ROC_inv_ross = diff(vI_ross_cycle); 
ROC_inv_kh = diff(vI_kh_cycle); 
ROC_inv_ea = diff(vI_ea_cycle); 
ROC_inv_wed = diff(vI_wed_cycle); 
ROC_inv_bell = diff(vI_bell_cycle); 
ROC_inv_tot = diff(vI_tot_cycle); 

ROC_2016_ross = diff(ross_SIE(:,year_index));
ROC_2016_kh = diff(kh_SIE(:,year_index));
ROC_2016_ea = diff(ea_SIE(:,year_index));
ROC_2016_wed = diff(wed_SIE(:,year_index));
ROC_2016_bell = diff(bell_SIE(:,year_index));
ROC_2016_tot = diff(tot_SIE(:,year_index));


figure(21)
subplot(2,3,1)
%plot(0:363,ROC_2016_ross, 'linewidth',2,'color','green')
%hold on
plot(0:363,ROC_trad_ross,'linewidth',2,'color',[1, 0.5, 0])
hold on 
plot(0:363,ROC_inv_ross, 'linewidth',2,'color','black')
hold on

yline(0, 'r', 'LineWidth', 1) 
ylim([-0.3 0.2])


title('Ross', 'interpreter', 'latex')
ylabel('Day-to-Day Change', 'Fontsize', 14, 'interpreter', 'latex')
set(gca, 'FontSize', 20)
set(gcf,'color','white')

subplot(2,3,2)
%plot(0:363,ROC_2016_kh, 'linewidth',2,'color','green') 
%hold on
plot(0:363,ROC_trad_kh,'linewidth',2,'color',[1, 0.5, 0])
hold on 
plot(0:363,ROC_inv_kh, 'linewidth',2,'color','black')

yline(0, 'r', 'LineWidth', 1) 
ylim([-0.3 0.2])


title('King Hakon', 'interpreter', 'latex')
set(gca, 'FontSize', 20)
set(gcf,'color','white')

subplot(2,3,3)
%plot(0:363,ROC_2016_ea, 'linewidth',2,'color','green')
%hold on
plot(0:363,ROC_trad_ea,'linewidth',2,'color',[1, 0.5, 0])
hold on 
plot(0:363,ROC_inv_ea, 'linewidth',2,'color','black')
hold on

yline(0, 'r', 'LineWidth', 1) 
ylim([-0.3 0.2])


title('East Antarctica', 'interpreter', 'latex')
set(gca, 'FontSize', 20)
set(gcf,'color','white')

subplot(2,3,4)
%plot(0:363,ROC_2016_wed, 'linewidth',2,'color','green')
%hold on
plot(0:363,ROC_trad_wed,'linewidth',2,'color',[1, 0.5, 0])
hold on 
plot(0:363,ROC_inv_wed, 'linewidth',2,'color','black')

yline(0, 'r', 'LineWidth', 1) 
ylim([-0.3 0.2])


title('Wedddell', 'interpreter', 'latex')
ylabel('Day-to-Day Change', 'Fontsize', 14, 'interpreter', 'latex')
xlabel('Day of Cycle', 'Fontsize', 14, 'interpreter', 'latex')
set(gca, 'FontSize', 20)
set(gcf,'color','white')


subplot(2,3,5)
%plot(0:363,ROC_2016_bell, 'linewidth',2,'color','green')
%hold on
plot(0:363,ROC_trad_bell,'linewidth',2,'color',[1, 0.5, 0])
hold on 
plot(0:363,ROC_inv_bell, 'linewidth',2,'color','black')

yline(0, 'r', 'LineWidth', 1) 
ylim([-0.3 0.2])


title('Amundsen-Bellingshausen', 'interpreter', 'latex')
xlabel('Day of Cycle', 'Fontsize', 14, 'interpreter', 'latex')
set(gca, 'FontSize', 20)
set(gcf,'color','white')


subplot(2,3,6)
%plot(0:363,ROC_2016_tot, 'linewidth',2,'color','green')
%hold on
plot(0:363,ROC_trad_tot,'linewidth',2,'color',[1, 0.5, 0])
hold on 
plot(0:363,ROC_inv_tot, 'linewidth',2,'color','black')

yline(0, 'r', 'LineWidth', 1) 
ylim([-0.3 0.2])

title('Total', 'interpreter', 'latex')
xlabel('Day of Cycle', 'Fontsize', 14, 'interpreter', 'latex')
set(gca, 'FontSize', 20)
set(gcf,'color','white')
legend('Traditional','Invariant')



%}


% Append the current array to the combined array
amplitude_adjusted_cycle_ross_all = [amplitude_adjusted_cycle_ross_all, aA_ross_cycle];
amplitude_adjusted_cycle_kh_all = [amplitude_adjusted_cycle_kh_all, aA_kh_cycle]; 
amplitude_adjusted_cycle_ea_all = [amplitude_adjusted_cycle_ea_all, aA_ea_cycle];
amplitude_adjusted_cycle_wed_all = [amplitude_adjusted_cycle_wed_all, aA_wed_cycle];
amplitude_adjusted_cycle_bell_all = [amplitude_adjusted_cycle_bell_all, aA_bell_cycle]; 
amplitude_adjusted_cycle_tot_all = [amplitude_adjusted_cycle_tot_all, aA_tot_cycle];


% Append the current array to the combined array
phase_adjusted_cycle_ross_all = [phase_adjusted_cycle_ross_all, vP_ross_cycle];
phase_adjusted_cycle_kh_all = [phase_adjusted_cycle_kh_all, vP_kh_cycle]; 
phase_adjusted_cycle_ea_all = [phase_adjusted_cycle_ea_all, vP_ea_cycle];
phase_adjusted_cycle_wed_all = [phase_adjusted_cycle_wed_all, vP_wed_cycle];
phase_adjusted_cycle_bell_all = [phase_adjusted_cycle_bell_all, vP_bell_cycle]; 
phase_adjusted_cycle_tot_all = [phase_adjusted_cycle_tot_all, vP_tot_cycle];


% Append the current array to the combined array
APAC_adjusted_cycle_ross_all = [APAC_adjusted_cycle_ross_all, vPa_ross_cycle];
APAC_adjusted_cycle_kh_all = [APAC_adjusted_cycle_kh_all, vPa_kh_cycle]; 
APAC_adjusted_cycle_ea_all = [APAC_adjusted_cycle_ea_all, vPa_ea_cycle];
APAC_adjusted_cycle_wed_all = [APAC_adjusted_cycle_wed_all, vPa_wed_cycle];
APAC_adjusted_cycle_bell_all = [APAC_adjusted_cycle_bell_all, vPa_bell_cycle]; 
APAC_adjusted_cycle_tot_all = [APAC_adjusted_cycle_tot_all, vPa_tot_cycle];





end

% Reshape arrays
amplitude_adjusted_cycle_ross = reshape(amplitude_adjusted_cycle_ross_all , [], 1);
amplitude_adjusted_cycle_kh = reshape(amplitude_adjusted_cycle_kh_all , [], 1);
amplitude_adjusted_cycle_ea = reshape(amplitude_adjusted_cycle_ea_all , [], 1);
amplitude_adjusted_cycle_wed = reshape(amplitude_adjusted_cycle_wed_all , [], 1);
amplitude_adjusted_cycle_bell = reshape(amplitude_adjusted_cycle_bell_all , [], 1);
amplitude_adjusted_cycle_tot = reshape(amplitude_adjusted_cycle_tot_all , [], 1);


% Organize data 
%{
trad_ross_data = [day_of_cycle' ; a_ross']; 
trad_kh_data = [day_of_cycle' ; a_kh']; 
trad_ea_data = [day_of_cycle' ; a_ea']; 
trad_wed_data = [day_of_cycle' ; a_wed']; 
trad_bell_data = [day_of_cycle' ; a_bell']; 
trad_tot_data = [day_of_cycle' ; a_tot']; 

inv_ross_data = [day_of_cycle' ; vI_ross']; 
inv_kh_data = [day_of_cycle' ; vI_kh']; 
inv_ea_data = [day_of_cycle' ; vI_ea']; 
inv_wed_data = [day_of_cycle' ; vI_wed']; 
inv_bell_data = [day_of_cycle' ; vI_bell']; 
inv_tot_data = [day_of_cycle' ; vI_tot']; 

%}
% Define corresponding years and day of cycle for amplitude adjusted cycle 
amplitude_years_all = []; 
xx_amplitude_adjusted_all = []; 
for i = 1979:2022
    amplitude_year_matrix = i*ones(1,365); 
    amplitude_years_all = [amplitude_years_all, amplitude_year_matrix]; 
    amplitude_day_of_cycle = 0:364; 
    xx_amplitude_adjusted_all = [xx_amplitude_adjusted_all , amplitude_day_of_cycle]; 
end 

%{
amp_ross_data = [amplitude_years_all' ; xx_amplitude_adjusted_all ; amplitude_adjusted_cycle_ross_all']; 
amp_kh_data = [amplitude_years_all' ; xx_amplitude_adjusted_all ; amplitude_adjusted_cycle_kh_all']; 
amp_ea_data = [amplitude_years_all' ; xx_amplitude_adjusted_all ; amplitude_adjusted_cycle_ea_all']; 
amp_wed_data = [amplitude_years_all' ; xx_amplitude_adjusted_all ; amplitude_adjusted_cycle_wed_all']; 
amp_bell_data = [amplitude_years_all' ; xx_amplitude_adjusted_all ; amplitude_adjusted_cycle_bell_all']; 
amp_tot_data = [amplitude_years_all' ; xx_amplitude_adjusted_all ; amplitude_adjusted_cycle_tot_all']; 
%}




%% Write data to csv

% Combine day, year, and cycle data data into matrix

traditional_cycle = [day_of_cycle'  a_ross_cycle  a_bell_cycle  a_wed_cycle  a_kh_cycle  a_ea_cycle  a_tot_cycle];  
invariant_cycle = [day_of_cycle' vI_ross_cycle vI_bell_cycle vI_wed_cycle vI_kh_cycle vI_ea_cycle vI_tot_cycle]; 
amplitude_cycle = [amplitude_years_all' xx_amplitude_adjusted_all' amplitude_adjusted_cycle_ross amplitude_adjusted_cycle_bell amplitude_adjusted_cycle_wed amplitude_adjusted_cycle_kh amplitude_adjusted_cycle_ea amplitude_adjusted_cycle_tot]; 

phase_ross_data = [phase_years_ross_all'  xx_phase_adjusted_cycle_ross_all'  phase_adjusted_cycle_ross_all']; 
phase_kh_data = [phase_years_kh_all'  xx_phase_adjusted_cycle_kh_all'  phase_adjusted_cycle_kh_all']; 
phase_ea_data = [phase_years_ea_all'  xx_phase_adjusted_cycle_ea_all'  phase_adjusted_cycle_ea_all']; 
phase_wed_data = [phase_years_wed_all'  xx_phase_adjusted_cycle_wed_all'  phase_adjusted_cycle_wed_all']; 
phase_bell_data = [phase_years_bell_all'  xx_phase_adjusted_cycle_bell_all'  phase_adjusted_cycle_bell_all']; 
phase_tot_data = [phase_years_tot_all'  xx_phase_adjusted_cycle_tot_all'  phase_adjusted_cycle_tot_all']; 


APAC_ross_data = [APAC_years_ross_all'  xx_APAC_adjusted_cycle_ross_all'  APAC_adjusted_cycle_ross_all']; 
APAC_kh_data = [APAC_years_kh_all'  xx_APAC_adjusted_cycle_kh_all'  APAC_adjusted_cycle_kh_all']; 
APAC_ea_data = [APAC_years_ea_all'  xx_APAC_adjusted_cycle_ea_all'  APAC_adjusted_cycle_ea_all']; 
APAC_wed_data = [APAC_years_wed_all'  xx_APAC_adjusted_cycle_wed_all'  APAC_adjusted_cycle_wed_all']; 
APAC_bell_data = [APAC_years_bell_all'  xx_APAC_adjusted_cycle_bell_all'  APAC_adjusted_cycle_bell_all']; 
APAC_tot_data = [APAC_years_tot_all'  xx_APAC_adjusted_cycle_tot_all'  APAC_adjusted_cycle_tot_all']; 


% Write the data to a CSV file
csvwrite('traditional_cycle.csv', traditional_cycle);
csvwrite('invariant_cycle.csv', invariant_cycle);
csvwrite('amplitude_adjusted_cycle.csv', amplitude_cycle);


csvwrite('phase_adjusted_cycle_ross.csv', phase_ross_data);
csvwrite('phase_adjusted_cycle_bell.csv', phase_bell_data);
csvwrite('phase_adjusted_cycle_wed.csv', phase_wed_data);
csvwrite('phase_adjusted_cycle_kh.csv', phase_kh_data);
csvwrite('phase_adjusted_cycle_ea.csv', phase_ea_data);
csvwrite('phase_adjusted_cycle_tot.csv', phase_tot_data);


csvwrite('APAC_ross.csv', APAC_ross_data);
csvwrite('APAC_bell.csv', APAC_bell_data);
csvwrite('APAC_wed.csv', APAC_wed_data);
csvwrite('APAC_kh.csv', APAC_kh_data);
csvwrite('APAC_ea.csv', APAC_ea_data);
csvwrite('APAC_tot.csv', APAC_tot_data);



% Create a README file
readme_text = [
    "This dataset contains the Traditional, Invariant, Amplitude-Adjusted, Phase-Adjusted, and Amplitude-and-Phase-Adjusted Cycles of  " + ...
    "Antarctica sea ice extent, calculated using data spanning 1979 to 2023. See Handcock and Raphael 2019. "
    ""
    "Files included:"
    " 1. traditional_cycle.csv: Contains the following columns:"
    "    Column 1 - Day of cycle"
    "    Column 2 - Traditional Cycle of Antarctic sea ice extent in the Ross Sea"
    "    Column 3 - Traditional Cycle of Antarctic sea ice extent in the Amundsen-Bellingshausen Seas"
    "    Column 4 - Traditional Cycle of Antarctic sea ice extent in the Weddell Sea"
    "    Column 5 - Traditional Cycle of Antarctic sea ice extent in the King Hakon Sea"
    "    Column 6 - Traditional Cycle of Antarctic sea ice extent in East Antarctica"
    "    Column 7 - Traditional Cycle of total (continental) Antarctic sea ice extent"
    ""
    " 2. invariant_cycle.csv: Contains the following columns:"
    "    Column 1 - Day of cycle"
    "    Column 2 - Invariant Cycle of Antarctic sea ice extent in the Ross Sea"
    "    Column 3 - Invariant Cycle of Antarctic sea ice extent in the Amundsen-Bellingshausen Seas"
    "    Column 4 - Invariant Cycle of Antarctic sea ice extent in the Weddell Sea"
    "    Column 5 - Invariant Cycle of Antarctic sea ice extent in the King Hakon Sea"
    "    Column 6 - Invariant Cycle of Antarctic sea ice extent in East Antarctica"
    "    Column 7 - Invariant Cycle of total (continental) Antarctic sea ice extent"
    ""
    " 3. amplitude_adjusted_cycle.csv: Contains the following columns:"
    "    Column 1 - Years"
    "    Column 2 - Day of cycle"
    "    Column 3 - Amplitude-Adjusted Cycle of Antarctic sea ice extent in the Ross Sea"
    "    Column 4 - Amplitude-Adjusted Cycle of Antarctic sea ice extent in the Amundsen-Bellingshausen Seas"
    "    Column 5 - Amplitude-Adjusted Cycle of Antarctic sea ice extent in the Weddell Sea"
    "    Column 6 - Amplitude-Adjusted Cycle of Antarctic sea ice extent in the King Hakon Sea"
    "    Column 7 - Amplitude-Adjusted Cycle of Antarctic sea ice extent in East Antarctica"
    "    Column 8 - Amplitude-Adjusted Cycle of total (continental) Antarctic sea ice extent"
    ""
    " 4. phase_adjusted_cycle_ross.csv: Contains the following columns:"
    "    Column 1 - Years"
    "    Column 2 - Day of cycle"
    "    Column 3 - Phase-Adjusted Cycle of Antarctic sea ice extent in the Ross Sea"
    ""
    " 5. phase_adjusted_cycle_bell.csv: Contains the following columns:"
    "    Column 1 - Years"
    "    Column 2 - Day of cycle"
    "    Column 3 - Phase-Adjusted Aycle of Antarctic sea ice extent in the Amundsen-Bellingshausen Seas"
    ""
    " 6. phase_adjusted_cycle_wed.csv: Contains the following columns:"
    "    Column 1 - Years"
    "    Column 2 - Day of cycle"
    "    Column 3 - Phase-Adjusted Cycle of Antarctic sea ice extent in the Weddell Sea"
    ""
    " 7. phase_adjusted_cycle_kh.csv: Contains the following columns:"
    "    Column 1 - Years"
    "    Column 2 - Day of cycle"
    "    Column 3 - Phase-Adjusted Cycle of Antarctic sea ice extent in the King Hakon Sea"
    ""
    " 8. phase_adjusted_cycle_ea.csv: Contains the following columns:"
    "    Column 1 - Years"
    "    Column 2 - Day of cycle"
    "    Column 3 - Phase-Adjusted Cycle of Antarctic sea ice extent in East Antarctica"
    ""
    " 9. phase_adjusted_cycle_tot.csv: Contains the following columns:"
    "    Column 1 - Years"
    "    Column 2 - Day of cycle"
    "    Column 3 - Phase-Adjusted Cycle of total (continental) Antarctic sea ice extent"
    ""
    " 10. APAC_ross.csv: Contains the following columns:"
    "    Column 1 - Years"
    "    Column 2 - Day of cycle"
    "    Column 3 - Amplitude-and-Phase-adjusted (APAC) Cycle of Antarctic sea ice extent in the Ross Sea"
    ""
    " 11. APAC_bell.csv: Contains the following columns:"
    "    Column 1 - Years"
    "    Column 2 - Day of cycle"
    "    Column 3 - Amplitude-and-Phase-adjusted (APAC) Cycle of Antarctic sea ice extent in the Amundsen-Bellingshausen Sea"
    ""
    " 12. APAC_wed.csv: Contains the following columns:"
    "    Column 1 - Years"
    "    Column 2 - Day of cycle"
    "    Column 3 - Amplitude-and-Phase-adjusted (APAC) Cycle of Antarctic sea ice extent in the Weddell Sea"
    ""
    " 13. APAC_kh.csv: Contains the following columns:"
    "    Column 1 - Years"
    "    Column 2 - Day of cycle"
    "    Column 3 - Amplitude-and-Phase-Adjusted (APAC) Cycle of Antarctic sea ice extent in the King Hakon Sea"
    ""
    " 14. APAC_ea.csv: Contains the following columns:"
    "    Column 1 - Years"
    "    Column 2 - Day of cycle"
    "    Column 3 - Amplitude-and-Phase-Adjusted (APAC) Cycle of Antarctic sea ice extent in East Antarctica"
    ""
    " 15. APAC_tot.csv: Contains the following columns:"
    "    Column 1 - Years"
    "    Column 2 - Day of cycle"
    "    Column 3 - Amplitude-and-Phase-Adjusted (APAC) Cycle of total (continental) Antarctic sea ice extent"
    ""

];

% Write the README file
fid = fopen('README.txt', 'wt');
fprintf(fid, '%s\n', readme_text{:});
fclose(fid);

% Optionally, package the files into a ZIP file
zip('ASIE_annual_cycle_csv.zip', {'traditional_cycle.csv','invariant_cycle.csv' , 'amplitude_adjusted_cycle.csv', ...
    'phase_adjusted_cycle_ross.csv', 'phase_adjusted_cycle_bell.csv', 'phase_adjusted_cycle_wed.csv', 'phase_adjusted_cycle_kh.csv', 'phase_adjusted_cycle_ea.csv', 'phase_adjusted_cycle_tot.csv', ...
    'APAC_ross.csv', 'APAC_bell.csv', 'APAC_wed.csv', 'APAC_kh.csv', 'APAC_ea.csv', 'APAC_tot.csv', 'README.txt'});




%}




