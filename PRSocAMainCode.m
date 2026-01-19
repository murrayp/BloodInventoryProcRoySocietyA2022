close all;clear all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set plotting parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fsize=30;
axis_l_width=5;
l_width=4;
set(0, 'DefaultAxesFontSize', fsize,'DefaultAxesLineWidth',axis_l_width,'DefaultLineLineWidth',l_width);
list_factory = fieldnames(get(groot,'factory'));
index_interpreter = find(contains(list_factory,'Interpreter'));
for i = 1:length(index_interpreter)
    default_name = strrep(list_factory{index_interpreter(i)},'factory','default');
    set(groot, default_name,'latex');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Switch to parameters from one of three different clinical settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k_1=3/7; % supply frequency
beta=0.001;
clinical_setting=1; % toggle between one of three clinical settings
q_target_WAPI=0.2;

switch(clinical_setting)
    case(1)
        r_used_vec=[18 12 14 6 0 2 0 0]; % restock thresholds
        mean_age_transfused_units=[20.1976   26.7500   17.8427   22.2558       NaN   16.4643       NaN       NaN]; % mean age of transfused time measured from inventories
        k_2_vec=[  1.5724    0.5000    1.0658    0.2368         0    0.1513         0         0]; % demand frequencies - measured from inventory        
        % KPI target
        ISI_target=10;
        target_mean_transfusion_age=16;
        %  For plotting 
        Tag='Case1'; %
    case(2)
        mean_age_transfused_units=[15.8457 22.3980 15.8978 15.7945 8.000 11.4400  NaN    NaN];
        Tag='Case2';
        k_2_vec=[    1.6283    0.5131    1.1780    0.3822    0.0105    0.2618         0         0];
        r_used_vec=[18 8 14 4 0 2 0 0];
        ISI_target=10;
        target_mean_transfusion_age=16;
     case(3)
        k_2_vec=[    8.8724    2.8930    6.0576    1.2593    2.4198    0.3292    0.4979    0.1523];
        r_used_vec=[55 30 45 11 15 6 6 4];
        mean_age_transfused_units=[11.9944   23.6984   14.2181   17.5882   14.0595   23.5500   20.5455   27.9459];
        k_1=6/7;
        Tag='Case3';
        beta=0.001;
        ISI_target=6;
        target_mean_transfusion_age=11;
end


% Supply summary stats
mu=5.65; %mean supply  age 
a_0=2.0; % min supply age
a_1=2*mu-a_0; % upper bound for uniform approximation of supply age distribution
A=35; % max age
eta=1; % Num units transfused in each demand event
init_a_star_guess=25;

bloodtypesstr={'$O^+$','$O^-$','$A^+$','$A^-$','$B^+$','$B^-$','$AB^+$','$AB^-$'};
a_guess=32; % used in optimisation
beta_vec=0.000005:0.00001:0.18;% used in optimisation
a_vec=5.01:0.1:A-0.1; %used in plotting

% Storage containers for quantities computed for each blood group
r_E_vec=zeros(size(k_2_vec));
r_w_vec=zeros(size(k_2_vec));
r_min_vec=zeros(size(k_2_vec));
r_alpha_vec=zeros(size(k_2_vec));
r_isi=zeros(size(k_2_vec));
WAPI_exact=zeros(size(k_2_vec));
WAPI=zeros(size(k_2_vec));
r_star_target_wapi=zeros(size(k_2_vec));
a_star_used_vec=zeros(size(k_2_vec));
a_star_lin_approx=zeros(size(k_2_vec));
expected_mean_age_transfused_units=zeros(size(k_2_vec));
a_E_beta=zeros(size(k_2_vec));
n_E_beta=zeros(size(k_2_vec));
% Loop over blood groups
for i=1:length(k_2_vec)
    
    k_2=k_2_vec(i);
    target_a_star=target_mean_transfusion_age;
    
    % only compute if demand frequency is large enough
    if k_2>0.05
        % Solve transcendental equation for a^*
        a_star=fsolve(@(x) exp(-k_2*(A-x))-(1-k_2*(1/k_1+x-(a_1+a_0)/2)/(r_used_vec(i)/eta-1)),init_a_star_guess,[]); % eq 2.6 in paper

        expected_mean_age_transfused_units(i)=a_star+1/k_2;
    else
        expected_mean_age_transfused_units(i)=NaN;

    end
    % Steady state v  
    r_star_vec=eta+k_2*eta/k_1*(1+k_1*(a_vec-(a_1+a_0)/2))./(1-exp(-k_2*(A-a_vec))); % rearrange 
    % eq 2.6 in paper
    % wastage rate
    waste_vec=k_2*eta./(exp(k_2*(A-a_vec))-1); % see eq B23 in appendix
    
    % Plot a^* v r
    figure(21)
    hold on
    plot(r_star_vec,a_vec)


    % Plot wastage rate v  r
    figure(22)
    hold on
    plot(r_star_vec,waste_vec)

    % mean supply age
    mu=(a_1+a_0)/2;

    r_w=k_2*eta*(A+1/k_1-mu); % Eq B17 in appendix
    r_min=k_2*eta*(a_1-mu+1/k_1)+eta; % Eq C10 in appendix
    

    %% Find optimal a star
    [a_E,~,flag]=fminsearch(@ESingleBlood,a_guess,[],k_1,k_2,eta,A,beta,a_1,a_0);
    if flag==1
        n_E=eta+k_2*eta/k_1*(1+k_1*(a_E-mu))/(1-exp(-k_2*(A-a_E)));
    else % return NaN if does not converge
        a_E=NaN;
        n_E=NaN;
    end
    
    if k_2>0
        % heat maps v beta/a_star/n_star
        beta_vec=0.0001:0.0001:0.001;
        E_mat_heatmap=NaN*ones(length(beta_vec),length(a_vec));
        
        for j=1:length(beta_vec)
            E_vec=ESingleBlood(a_vec,k_1,k_2,eta,A,beta_vec(j),a_1,a_0); %=[0.0001:0.0001:0.001];
            E_mat_heatmap(j,:)=E_vec;
            a_star_opt=fminsearch(@ESingleBlood,a_guess,[],k_1,k_2,eta,A,beta_vec(j),a_1,a_0);
            a_E_beta(j)=a_star_opt;
            n_E_beta(j)=eta+k_2*eta/k_1*(1+k_1*(a_star_opt-mu))/(1-exp(-k_2*(A-a_star_opt)));
        end
    
    end

    %astar_data
    r_alpha=r_w-k_2*eta*(A-target_mean_transfusion_age); %k_2*eta*(1/k_1-mu+target_mean_transfusion_age);
    
    %Store various quantities in arrays for plotting later
    r_E_vec(i)=n_E;
    r_w_vec(i)=r_w;
    r_min_vec(i)=r_min;
    r_alpha_vec(i)=r_alpha;
    r_isi(i)=k_2*eta*(ISI_target+1/k_1); % Eq C12
    WAPI_exact(i)=100*exp(-k_2*(A-a_star)); %Eq C7
    WAPI(i)=100*exp(-k_2*(A+1/k_1-mu))*exp(r_used_vec(i)/eta-1);
    r_star_target_wapi(i)= r_w+eta*(1-log(100/q_target_WAPI));
    a_star_used_vec(i)=a_star;
    a_star_lin_approx(i)=1/k_2*(r_used_vec(i)/eta)-1/k_1+mu;
end

%% Steady state age distribution

figure(21);
box on
xlim([0 60])

xlabel('$r$')
ylabel('$a^*$')
legend(bloodtypesstr{[1:4 6]})
str=['Figures/Figure4a' Tag '.png'];
exportgraphics(gcf,str)

figure(22);
box on
xlim([0 60])

%bar([r_min_vec;r_used_vec;r_alpha_vec;]')
xlabel('$r$')
ylabel('$w$')
legend(bloodtypesstr{[1:4 6]})

%set(gca,'XTickLabel',bloodtypesstr,'TickLabelInterpreter','Latex')
str=['Figures/Figure4b' Tag '.png'];
exportgraphics(gcf,str)


figure;
v=[r_used_vec;r_isi;r_alpha_vec];
bar([r_used_vec;r_isi;r_alpha_vec]')
ylabel('$r$')
set(gca,'XTickLabel',bloodtypesstr,'TickLabelInterpreter','Latex')
legend('$r_{Used}$','$r_{ISI}$','$r_{age}$')
str=['Figures/Figure6' Tag '.png'];
ylim([0 1.05*max(v(:))])
exportgraphics(gcf,str)

figure;
bar([expected_mean_age_transfused_units;mean_age_transfused_units]')
ylabel('$a_{exp}$')
set(gca,'XTickLabel',bloodtypesstr,'TickLabelInterpreter','Latex')
legend('Model','Data')
str=['Figures/Figure3' Tag '.png'];
ylim([0 30])
exportgraphics(gcf,str)




%% Plot all data points togethers: gathers data from three different scenarios

figure
k_2_all=[ 1.5724    0.5000    1.0658    0.2368         0    0.1513         0         0  1.6283    0.5131    1.1780    0.3822    0.0105    0.2618         0         0  8.8724    2.8930    6.0576    1.2593    2.4198    0.3292    0.4979    0.1523];
mean_age_transfused_units_all=[20.1976   26.7500   17.8427   22.2558       NaN   16.4643       NaN       NaN 15.8457 22.3980 15.8978 15.7945 8.000 11.4400  NaN    NaN 11.9944   23.6984   14.2181   17.5882   14.0595   23.5500   20.5455   27.9459];
r_used_vec_all=[18 12 14 6 0 2 0 0 18 8 14 4 0 2 0 0 55 30 45 11 15 6 6 4 ];


% supply frequencies in three different scenarios
k_1_vec=1/7*[3*ones(1,8) 3*ones(1,8) 6*ones(1,8) ];

% Solve model to compute mean expected age of tranfused unit for each case
expected_mean_age_transfused_units_all=NaN*ones(size(mean_age_transfused_units_all));
for i=1:length(k_2_all)
    a_star=fsolve(@(x) exp(-k_2_all(i)*(A-x))-(1-k_2_all(i)*(1/k_1_vec(i)+x-(a_1+a_0)/2)/(r_used_vec_all(i)/eta-1)),25,[]);
    expected_mean_age_transfused_units_all(i)=a_star+1/k_2_all(i);   
end


% Gather together measurements and model predicts in single 2D array
data=[expected_mean_age_transfused_units_all'  mean_age_transfused_units_all' ];

% annotation
period={'1','1','1','1','1','1','1','1','2','2','2','2','2','2','2','2','3','3','3','3','3','3','3','3'};
group={'$O^+$','$O^-$','$A^+$','$A^-$','$B^+$','$B^-$','$AB^+$','$AB^-$','$O^+$','$O^-$','$A^+$','$A^-$','$B^+$','$B^-$','$AB^+$','$AB^-$','$O^+$','$O^-$','$A^+$','$A^-$','$B^+$','$B^-$','$AB^+$','$AB^-$'};

col1=data(:,1);
col2=data(:,2);

% remove NANs
notnanind=find((isnan(col1)==0) & (isinf(col1)==0)&col1<35&col2<35);
data=data(notnanind,:);
period=period(notnanind);
group=group(notnanind);


% Deviation between model and observation
z_all=data(:,1)-data(:,2);


Z_score=(z_all-mean(z_all))/(std(z_all)) ;

outliers=find(abs(Z_score)>2);
data_min_outliers=data;
data_min_outliers(outliers,:)=[];
data_outliers=data(outliers,:);

z=data_min_outliers(:,1)-data_min_outliers(:,2);
tbl = table(z);
mdl = fitlm(zeros(size(z)), z);   % intercept-only model
mdl = fitlm(tbl, 'z ~ 1');

y_hat = data_min_outliers(:,1) + mdl.Coefficients.Estimate(1);
y=data_min_outliers(:,1) ;
% R^2
SS_res = sum((y - y_hat).^2);
SS_tot = sum((y - mean(y)).^2);

R2 = 1 - SS_res / SS_tot;

figure;plot(Z_score,'x')
hold on;
plot(Z_score)
vec=[1 size(Z_score,1)];
plot(vec,-2*ones(size(vec)),vec,2*ones(size(vec)))
xlabel('Data index')
ylabel('Z score')

%% PLot non-outlier datapoints
figure; 
ax1 = axes;
hold on
box on
axis square

ylim([0 35])
xlim([7.5 35])
plot([0 35],[0 35],'LineWidth',1)
text(12,32,['$R^2$ = ' num2str(R2) ',   y=' num2str(mdl.Coefficients.Estimate(1)) '+ x' ],'Fontsize',16)
colorLevels  = unique(period);
markerLevels = bloodtypesstr;

% Plot regression line
x_plt=[0 35];
y_plt=-mdl.Coefficients.Estimate(1)+x_plt;
plot(x_plt,y_plt,':','LineWidth',1)


% Scatter plot each blood groups

colorMap = lines(numel(colorLevels));
markers  = {'o','s','^','d','v','>','p','h'};
x=data(:,1);
y=data(:,2);

for i = 1:numel(colorLevels)
    for j = 1:numel(markerLevels)
        idx = (strcmp(period,colorLevels{i}) & strcmp(group,markerLevels{j}));
        scatter(ax1,x(idx), y(idx),60,colorMap(i,:),markers{j},'filled');
    end
end


% plot outliers 
plot(data_outliers(:,1),data_outliers(:,2),'o','MarkerSize', 14, ...
    'LineWidth', 1.5, ...
    'MarkerEdgeColor', 'r', ...
    'LineStyle', 'none')

ylabel('${a}_{exp}$ (Data)')
xlabel('$a_{exp}$ (Model)')

%%  Add legend

axColor = axes('Position',ax1.Position,'Visible','off');
hold(axColor,'on')

hColor = gobjects(numel(colorLevels),1);
for i = 1:numel(colorLevels)
    hColor(i) = scatter(axColor,nan, nan, 60, colorMap(i,:), 'o', 'filled');
end

lgd1 = legend(axColor,hColor,{'DGH P1','DGH P2','RH'},'Location','northeastoutside','FontSize',20);
title(lgd1,'Dataset')
set(lgd1,'AutoUpdate','off')
hold on

%
axMarker = axes('Position',ax1.Position,'Visible','off');
hold(axMarker,'on')

hMarker = gobjects(numel(markerLevels),1);
for j = 1:numel(markerLevels)
    hMarker(j) = scatter(axMarker,nan, nan, 60, 'k', markers{j}, 'filled');
end

lgd2 = legend(axMarker,hMarker, markerLevels,'Location','southeastoutside','FontSize',20);
title(lgd2,'Group')
set(lgd2,'AutoUpdate','off')

drawnow                       % ensure positions are finalized

pos = lgd2.Position;          % [x y width height]
pos(2) = pos(2) - 0.06;       % move DOWN (adjust value as needed)
lgd2.Position = pos;

pos = lgd1.Position;          % [x y width height]
pos(2) = pos(2) + 0.06;       % move DOWN (adjust value as needed)
lgd1.Position = pos;



str='Figures/Figure3c.png';

exportgraphics(gcf,str)
pause(0.1)



ISI=r_used_vec./k_2_vec-1/k_1;

%%
x=r_w_vec -r_used_vec;
y=log(100./WAPI);
figure;
plot(x,x+1,'k',x,y,'rx','Linewidth',4)
xlabel('$n_{max}-r$')
ylabel('$\log\left(100/w_{API} \right)$')
str=['Figures/Figure7' Tag '.png'];
exportgraphics(gcf,str)

%%

Data_mat=[r_used_vec;k_2_vec;WAPI;ISI;r_min_vec;r_w_vec;r_E_vec;r_star_target_wapi;r_alpha_vec;r_isi];

(Data_mat(:,:));
%% 
%% 