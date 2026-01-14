close all;clear all;

% Set plotting parameters

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

clinical_setting=1;
k_1=3/7;
beta=0.001;

switch(clinical_setting)
    case(1)
        n_star_used_vec=[18 12 14 6 0 2 0 0];
        mean_age_transfused_units=[20.1976   26.7500   17.8427   22.2558       NaN   16.4643       NaN       NaN];
        Tag='Case1'
        k_2_vec=[  1.5724    0.5000    1.0658    0.2368         0    0.1513         0         0]
        ISI_target=10
        target_mean_transfusion_age=16;

        %k_2_vec=[1.5150 0.5030 1.0659 0.2575 0  0.1677  0     0];
    case(2)
        mean_age_transfused_units=[15.8457 22.3980 15.8978 15.7945 8.000 11.4400  NaN    NaN]
        Tag='Case2'
        k_2_vec=[    1.6283    0.5131    1.1780    0.3822    0.0105    0.2618         0         0];
        n_star_used_vec=[18 8 14 4 0 2 0 0];
        ISI_target=10
        target_mean_transfusion_age=16;


     case(3)
        k_2_vec=[    8.8724    2.8930    6.0576    1.2593    2.4198    0.3292    0.4979    0.1523];
        n_star_used_vec=[55 30 45 11 15 6 6 4];
        mean_age_transfused_units=[11.9944   23.6984   14.2181   17.5882   14.0595   23.5500   20.5455   27.9459]
        k_1=6/7;
        Tag='Case3'
        beta=0.001;
        ISI_target=6
        target_mean_transfusion_age=11;

end

% Load demand data
% Load supply data
mu=5.65;
a_0=2.0;
a_1=2*mu-a_0; 

A=35;
eta=1;

% Determine model params
bloodtypesstr={'$O^+$','$O^-$','$A^+$','$A^-$','$B^+$','$B^-$','$AB^+$','$AB^-$'}
a_guess=32;
beta_vec=[0.000005:0.00001:0.18];
a_vec=5.01:0.1:A-0.1;

% Loop over blood types
for i=1:length(k_2_vec)
    
    k_2=k_2_vec(i);
    target_a_star=target_mean_transfusion_age;

    a_1_factor=(a_1/(a_1-a_0));
    if k_2>0.05
        a_star=fsolve(@(x) exp(-k_2*(A-x))-(1-k_2*(1/k_1+x-(a_1+a_0)/2)/(n_star_used_vec(i)/eta-1)),25,[]);
        expected_mean_age_transfused_units(i)=a_star+1/k_2;
    else
        expected_mean_age_transfused_units(i)=NaN;

    end

    n_star_vec=eta+k_2*eta/k_1*(1+k_1*(a_vec-(a_1+a_0)/2))./(1-exp(-k_2*(A-a_vec)));
    waste_vec=k_2*eta./(exp(k_2*(A-a_vec))-1);
    figure(21)
    hold on
    plot(n_star_vec,a_vec)
    xlim([0 60])

    figure(22)
    hold on
    plot(n_star_vec,waste_vec)
    xlim([0 60])



    % PLot q against a - discrete v ctm

    %n_w

    mu=(a_1+a_0)/2;
    n_w=k_2*eta*(A+1/k_1-mu);
    %a_star_model
    %n_min
    n_min=k_2*eta*(a_1-mu+1/k_1)+eta;
    %
    %n_E
    %%
    [a_E,~,flag]=fminsearch(@ESingleBlood,a_guess,[],k_1,k_2,eta,A,beta,a_1,a_0);
    if flag==1
        n_E=eta+k_2*eta/k_1*(1+k_1*(a_E-mu))/(1-exp(-k_2*(A-a_E)));
    else
        a_E=NaN;
        n_E=NaN;
    end
    

    if k_2>0
        % heat maps v beta/a_star/n_star
        beta_vec=[0.0001:0.0001:0.001];
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
    n_alpha=n_w-k_2*eta*(A-target_mean_transfusion_age); %k_2*eta*(1/k_1-mu+target_mean_transfusion_age);

    n_E_vec(i)=n_E;
    n_w_vec(i)=n_w;
    n_min_vec(i)=n_min;
    n_alpha_vec(i)=n_alpha;
    n_isi(i)=k_2*eta*(ISI_target+1/k_1);
    WAPI_exact(i)=100*exp(-k_2*(A-a_star));
    WAPI(i)=100*exp(-k_2*(A+1/k_1-mu))*exp(n_star_used_vec(i)/eta-1);
    q_target_WAPI=0.2;
    %n_star_target_wapi(i)= eta*(1+k_2/k_1+k_2*(A-(a_0+a_1)/2)-log(100/q_target_WAPI))
    n_star_target_wapi(i)= n_w+eta*(1-log(100/q_target_WAPI));
    a_star_used(i)=a_star;
    a_star_lin_approx(i)=1/k_2*(n_star_used_vec(i)/eta-1)-1/k_1+mu;
end

% PLot Compare n stars
%% 
%% 

%%
figure(21);
box on
%bar([n_min_vec;n_star_used_vec;n_alpha_vec;]')
xlabel('$r$')
ylabel('$a^*$')
legend(bloodtypesstr{[1:4 6]})

%set(gca,'XTickLabel',bloodtypesstr,'TickLabelInterpreter','Latex')
str=['Figures/Figure4a' Tag '.png']
exportgraphics(gcf,str)
pause(0.1)

figure(22);
box on
%bar([n_min_vec;n_star_used_vec;n_alpha_vec;]')
xlabel('$r$')
ylabel('$w$')
legend(bloodtypesstr{[1:4 6]})

%set(gca,'XTickLabel',bloodtypesstr,'TickLabelInterpreter','Latex')
str=['Figures/Figure4b' Tag '.png']
exportgraphics(gcf,str)
pause(0.1)


figure;
v=[n_star_used_vec;n_isi;n_alpha_vec]
bar([n_star_used_vec;n_isi;n_alpha_vec]')
ylabel('$r$')
set(gca,'XTickLabel',bloodtypesstr,'TickLabelInterpreter','Latex')
legend('$r_{Used}$','$r_{ISI}$','$r_{age}$')
str=['Figures/Figure6' Tag '.png']
ylim([0 1.05*max(v(:))])
exportgraphics(gcf,str)
pause(0.1)

figure;
bar([expected_mean_age_transfused_units;mean_age_transfused_units]')
ylabel('$a_{exp}$')
set(gca,'XTickLabel',bloodtypesstr,'TickLabelInterpreter','Latex')
legend('Model','Data')
str=['Figures/Figure3' Tag '.png']
ylim([0 30])
exportgraphics(gcf,str)
pause(0.1)



%%
figure
k_2_all=[ 1.5724    0.5000    1.0658    0.2368         0    0.1513         0         0  1.6283    0.5131    1.1780    0.3822    0.0105    0.2618         0         0  8.8724    2.8930    6.0576    1.2593    2.4198    0.3292    0.4979    0.1523]
mean_age_transfused_units_all=[20.1976   26.7500   17.8427   22.2558       NaN   16.4643       NaN       NaN 15.8457 22.3980 15.8978 15.7945 8.000 11.4400  NaN    NaN 11.9944   23.6984   14.2181   17.5882   14.0595   23.5500   20.5455   27.9459]
n_star_used_vec_all=[18 12 14 6 0 2 0 0 18 8 14 4 0 2 0 0 55 30 45 11 15 6 6 4 ]


k_1_vec=1/7*[3*ones(1,8) 3*ones(1,8) 6*ones(1,8) ]


expected_mean_age_transfused_units_all=NaN*ones(size(mean_age_transfused_units_all))
for i=1:length(k_2_all)
    a_star=fsolve(@(x) exp(-k_2_all(i)*(A-x))-(1-k_2_all(i)*(1/k_1_vec(i)+x-(a_1+a_0)/2)/(n_star_used_vec_all(i)/eta-1)),25,[]);
    expected_mean_age_transfused_units_all(i)=a_star+1/k_2_all(i);
   
end


data=[expected_mean_age_transfused_units_all'  mean_age_transfused_units_all' ]

period={'1','1','1','1','1','1','1','1','2','2','2','2','2','2','2','2','3','3','3','3','3','3','3','3'}
group={'$O^+$','$O^-$','$A^+$','$A^-$','$B^+$','$B^-$','$AB^+$','$AB^-$','$O^+$','$O^-$','$A^+$','$A^-$','$B^+$','$B^-$','$AB^+$','$AB^-$','$O^+$','$O^-$','$A^+$','$A^-$','$B^+$','$B^-$','$AB^+$','$AB^-$'}

col1=data(:,1)
col2=data(:,2)

notnanind=find((isnan(col1)==0) & (isinf(col1)==0)&col1<35&col2<35)
data=data(notnanind,:)
period=period(notnanind)
group=group(notnanind)


data_min_outliers=data
outliers=[4,12]
data_min_outliers(outliers,:)=[]

data_outliers=data(outliers,:)


z_all=data(:,1)-data(:,2)
z=data_min_outliers(:,1)-data_min_outliers(:,2)
mdl = fitlm(zeros(size(z)), z);   % intercept-only model


Z_score=(z_all-mean(z_all))/(std(z_all)) 

R2 = mdl.Rsquared.Ordinary;
R2_adj = mdl.Rsquared.Adjusted

y_hat = data_min_outliers(:,1) + mdl.Coefficients.Estimate(1);
y=data_min_outliers(:,1) 
% R^2
SS_res = sum((y - y_hat).^2);
SS_tot = sum((y - mean(y)).^2);

R2 = 1 - SS_res / SS_tot;

figure;plot(Z_score,'x')
hold on;
plot(Z_score)
vec=[1 size(Z_score,1)]
plot(vec,-2*ones(size(vec)),vec,2*ones(size(vec)))

%%
figure; 
ax1 = axes;
hold on
box on
axis square

ylim([0 35])
xlim([7.5 35])
plot([0 35],[0 35],'LineWidth',1)
mdl.Coefficients
text(12,32,['$R^2$ = ' num2str(R2) ',   y=' num2str(mdl.Coefficients.Estimate(1)) '+ x' ],'Fontsize',16)
colorLevels  = unique(period);
markerLevels = bloodtypesstr;


x_plt=[0 35]
y_plt=-mdl.Coefficients.Estimate(1)+x_plt
plot(x_plt,y_plt,':','LineWidth',1)





colorMap = lines(numel(colorLevels));
markers  = {'o','s','^','d','v','>','p','h'};
%o  s  d  ^  v  <  >  p  h
x=data(:,1);
y=data(:,2);

for i = 1:numel(colorLevels)
    for j = 1:numel(markerLevels)
        idx = (strcmp(period,colorLevels{i}) & strcmp(group,markerLevels{j}));

        scatter(ax1,x(idx), y(idx),60,colorMap(i,:),markers{j},'filled');
    end
end

plot(data_outliers(:,1),data_outliers(:,2),'o','MarkerSize', 14, ...
    'LineWidth', 1.5, ...
    'MarkerEdgeColor', 'r', ...
    'LineStyle', 'none')



ylabel('${a}_{exp}$ (Data)')
xlabel('$a_{exp}$ (Model)')

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



str=['Figures/Figure3c.png']

exportgraphics(gcf,str)
pause(0.1)



ISI=n_star_used_vec./k_2_vec-1/k_1

%%
x=n_w_vec -n_star_used_vec;
y=log(100./WAPI);
figure;
plot(x,x+1,'k',x,y,'rx','Linewidth',4)
xlabel('$n_{max}-r$')
ylabel('$\log\left(100/w_{API} \right)$')
str=['Figures/Figure7' Tag '.png']
exportgraphics(gcf,str)

%%

Data_mat=[n_star_used_vec;k_2_vec;WAPI;ISI;n_min_vec;n_w_vec;n_E_vec;n_star_target_wapi;n_alpha_vec;n_isi]

(Data_mat(:,:))
%% 
%% 