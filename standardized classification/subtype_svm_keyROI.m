clc;
clear;

dataset={'train' ,'validate', 'test'}; % discovery / replication-1 / replication-2

Schaefer2018_416=readtable('Schaefer2018_416Parcels_17Networks.csv');
region_name=Schaefer2018_416.ROIName;
%key regions
load('smileGAN_generated\selectROIsByGradient.mat');
%subregion key ROI
KeyROI1=GM1_ROI{1};     
KeyROI2=GM2_ROI{1};
KeyROI2.ROIname{1}(strfind(KeyROI2.ROIname{1},'-'))='_';    
KeyROI=[KeyROI1 ; KeyROI2];


%------------------------------------discovery------------------------------------
demo_discovery=readtable('Demo_STB.xlsx');
cluster_result2=readtable('clustering_result.csv');
fc_discovery=readtable('FC_STB.xlsx');
       
[~,i1,i2]=intersect(demo_discovery.DP,cluster_result2.participant_id);
demo_discovery=demo_discovery(i1,:); fc_discovery=fc_discovery(i1,:);
cluster_result2=cluster_result2(i2,:);
demo_discovery.cluster_labelbyWT=cluster_result2.cluster_label;
demo_discovery.P1byWT=cluster_result2.p1;
demo_discovery.P2byWT=cluster_result2.p2;

demo1_discovery=demo_discovery(demo_discovery.cluster_labelbyWT==1,:);
fc1_discovery=fc_discovery(demo_discovery.cluster_labelbyWT==1,:);
demo2_discovery=demo_discovery(demo_discovery.cluster_labelbyWT==2,:);
fc2_discovery=fc_discovery(demo_discovery.cluster_labelbyWT==2,:);


clear X_disc y_disc cov_disc;
%X
X_disc=[fc1_discovery(: , KeyROI.index+2); fc2_discovery(: , KeyROI.index+2)];
X_disc=table2array(X_disc);
%Y
y_disc =[ones(height(demo1_discovery),1)-1; ones(height(demo2_discovery),1)] ;
%covariate
cov_disc=[[demo1_discovery.age(:), demo1_discovery.sex(:) , demo1_discovery.education(:)] ; ...
    [demo2_discovery.age(:), demo2_discovery.sex(:) , demo2_discovery.education(:)] ];

fprintf('Discovery sample size: %d\n', size(X_disc, 1));
fprintf('Number of SVM features: %d\n', size(X_disc, 2));


%------------------------------------replication-1------------------------------------
clear X_rep1 y_rep1 cov_rep1;
demo1_replication1=readtable('replication\siemens\demo\Visual_STB.xlsx');
demo2_replication1=readtable('replication\siemens\demo\DMN_CEN_STB.xlsx');

fc1_replication1=readtable('replication\siemens\FC\Visual_STB.xlsx');
fc2_replication1=readtable('replication\siemens\FC\DMN_CEN_STB.xlsx');
X_rep1=[fc1_replication1(: , KeyROI.index+2) ; fc2_replication1(: , KeyROI.index+2)];
X_rep1=table2array(X_rep1);

%Y
y_rep1 =[ones(height(demo1_replication1),1)-1; ones(height(demo2_replication1),1)] ;
%covariate
cov_rep1=[[demo1_replication1.age(:), demo1_replication1.sex(:) , demo1_replication1.education(:)] ; ...
    [demo2_replication1.age(:), demo2_replication1.sex(:) , demo2_replication1.education(:)] ];

fprintf('Replication-1 sample size: %d\n', size(X_rep1, 1));


%------------------------------------replication-2------------------------------------
clear X_rep2;
demo1_replication2=readtable('replication\GE\demo\Visual_STB.xlsx');
demo2_replication2=readtable('replication\GE\demo\DMN_CEN_STB.xlsx');

fc1_replication2=readtable('replication\GE\FC\Visual_STB.xlsx');
fc2_replication2=readtable('replication\GE\FC\DMN_CEN_STB.xlsx');
X_rep2=[fc1_replication2(: , KeyROI.index+2) ; fc2_replication2(: , KeyROI.index+2)];
X_rep2=table2array(X_rep2);

%Y
y_rep2 =[ones(height(demo1_replication2),1)-1; ones(height(demo2_replication2),1)] ;
%covariate
cov_rep2=[[demo1_replication2.age(:), demo1_replication2.sex(:) , demo1_replication2.education(:)] ; ...
    [demo2_replication2.age(:), demo2_replication2.sex(:) , demo2_replication2.education(:)] ];

fprintf('Replication-2 sample size: %d\n', size(X_rep2, 1));

%%  SVM settings
C_grid = [0.001, 0.01, 0.1, 1, 10, 100];
outerK = 5;
innerK = 5;
% For speed, you may first set nPerm = 1000.
% For manuscript-level final analysis, use nPerm = 10000.
nPerm = 10000;


%%   Nested CV in discovery cohort
[cv_results, cv_summary] = nested_cv_svm(X_disc, y_disc, cov_disc, C_grid, outerK, innerK);

disp(' ');
disp('===== Discovery nested cross-validation results =====');
disp(cv_results);

disp(' ');
disp('===== Discovery nested cross-validation summary =====');
disp(cv_summary);

%%  5. Train final locked classifier in full discovery cohort

best_C = choose_best_C_inner_cv(X_disc, y_disc, cov_disc, C_grid, innerK);
fprintf('\nBest C selected from full discovery cohort: %.4f\n', best_C);
final_model = train_preprocessed_svm(X_disc, y_disc, cov_disc, best_C);


%%  6. Evaluate final classifier

result_disc = evaluate_preprocessed_svm( final_model, X_disc, y_disc, cov_disc, nPerm, 'Discovery full sample');
result_rep1 = evaluate_preprocessed_svm( final_model, X_rep1, y_rep1, cov_rep1, nPerm, 'Replication-1');
result_rep2 = evaluate_preprocessed_svm( final_model, X_rep2, y_rep2, cov_rep2, nPerm, 'Replication-2');


%%  7. Save results
%  -----------------------------

all_results = struct();
all_results.discovery_nested_cv = cv_results;
all_results.discovery_nested_cv_summary = cv_summary;
all_results.discovery_full = result_disc;
all_results.replication1 = result_rep1;
all_results.replication2 = result_rep2;
all_results.best_C = best_C;
all_results.final_model = final_model;

save('subtype_svm_results.mat', 'all_results');

fprintf('\nSaved results to subtype_svm_results.mat\n');



