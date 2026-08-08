clc;
clear;

dataset={'discovery'  , 'replication1' ,  'replication2'};

clear regional_mytstat regional_mypval regional_result_beforeFDR regional_result_FDR sigregs_SANSA;
for d=1:length(dataset)
    
    demo_STB=readtable('Demo_STB.xlsx');
    fc_STB=readtable('FC_STB.xlsx'); 
    cluster_result2=readtable('clustering_result.csv');
    
    demo_nSTB=readtable('Demo_nSTB.xlsx');
    fc_nSTB=readtable('FC_nSTB.xlsx');   

        
    %将smile-GAN的分类结果融进demo中
    [~,i1,i2]=intersect(demo_STB.DP,cluster_result2.participant_id);
    demo_STB=demo_STB(i1,:); fc_STB=fc_STB(i1,:);
    cluster_result2=cluster_result2(i2,:);
    demo_STB.cluster_labelbyWT=cluster_result2.cluster_label;
    demo_STB.P1byWT=cluster_result2.p1;
    demo_STB.P2byWT=cluster_result2.p2;

    %分簇的个数
    nSubtype=height(tabulate(demo_STB.cluster_labelbyWT));
    
    clear SAindex SAdemo SAfc NSAindex NSAdemo NSAfc;
    for c=1:nSubtype
        %分簇后的人口学
        clear DEMO_STB_c FC_STB_c;
        switch c
            case 1
                DEMO_STB_c=demo_STB(demo_STB.cluster_labelbyWT==1,:);
                FC_STB_c=fc_STB(demo_STB.cluster_labelbyWT==1,:);
            case 2
                DEMO_STB_c=demo_STB(demo_STB.cluster_labelbyWT==2,:);
                FC_STB_c=fc_STB(demo_STB.cluster_labelbyWT==2,:);
           
        end
        

        

        clear convaribales Group y;
        convaribales.v1=[demo_nSTB.age ; DEMO_STB_c.age];
        convaribales.v2=[demo_nSTB.sex ; DEMO_STB_c.sex];
        convaribales.v3=[demo_nSTB.education ; DEMO_STB_c.education];
        Group=[ones(height(demo_nSTB),1); 2*ones(height(DEMO_STB_c),1)];%LS:1 ;HS:2
        y=cell2mat(table2cell(vertcat(fc_nSTB(:,3:end),FC_STB_c(:,3:end))));
        
        [regional_mytstat{d}{c},regional_mypval{d}{c},regional_result_beforeFDR{d}{c},regional_result_FDR{d}{c},sigregs_SANSA{d}{c}]=regionaltmap(convaribales,Group,...
            y,region_name);


    end

end


%% discovery与replication的t-map的相关性
clear rp;
for c=1:nSubtype
    [rp(c,1),rp(c,2)] = perm_sphere_p(regional_mytstat{1}{c}',regional_mytstat{2}{c}',perm_id,'Pearson');
end





















