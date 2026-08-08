clc;
clear;

dataset={'discovery'  , 'replication1' ,  'replication2'};
load('spintest.mat');%spin test 的perm_id
roiXgene=load('roiXgene.mat');
X=roiXgene.expressionROI(:,2:end);
X=zscore(X);% Predictors
gene_Symbol=roiXgene.probeInformation.GeneSymbol;
%% gene analysis
for d=1:length(dataset)
    for c=1:nSubtype
        
        Y=regional_mytstat{d}{c}'; % Response variable
        Y=zscore(Y);
        %% 查看每个pls component对Y的explained 占比情况
        %perform full PLS and plot variance in Y explained by top 15 components typically top 2 or 3 components will explain a large part of the variance
        %(hopefully!)
        [XL,YL,XS,YS,BETA,PCTVAR,MSE,stats]=plsregress(X,Y);
        dim=15;
        figure;
        plot(1:dim,cumsum(100*PCTVAR(2,1:dim)),'-o','LineWidth',1.5,'Color',[140/255,0,0]);
        ylim([0 100])
        set(gca,'LineWidth',2.7,'FontSize',14,'FontName','Times New Roman');
        xlabel('Number of PLS components','FontSize',20,'Fontname','Times New Roman');
        ylabel('Percent Variance Explained in Y','FontSize',20,'Fontname','Times New Roman');
        grid on
        set(gcf,'color','w');
        box off  ;
    
        % permutation testing to assess significance of PLS result as a function of the number of components (dim) included:
        R_p=zeros(dim,2);
        rep=10000;
        for dim=1:2
            [~,~,~,~,~,PCTVAR,~,~]=plsregress(X,Y,dim);
            temp=cumsum(100*PCTVAR(2,1:dim));
            Rsquared = temp(dim);
            clear Rsq;
            for j=1:rep
                %j
                order=randperm(size(Y,1));
                Yp=Y(order,:);
        
                [~,~,~,~,~,PCTVAR,~,~]=plsregress(X,Yp,dim);
                temp=cumsum(100*PCTVAR(2,1:dim));
                Rsq(j) = temp(dim);
            end
            R_perm{dim}=Rsq;
            
            dim
            R_p(dim,1)=Rsquared
            R_p(dim,2)=length(find(Rsq>=Rsquared))/rep
        end
        %% Do PLS in 2 dimensions (with 2 components):
        dim=2;
        [XL,YL,XS,YS,BETA,PCTVAR,MSE,stats]=plsregress(X,Y,dim);
        %store regions' IDs and weights in descending order of weight for both components:
        % PLS1   vs.  case-control MSN t-values
        [R1,p1]=corr(XS(:,1),Y)%主成分XS
        % [R1,p1]=corr([XS(:,1),XS(:,2)],Y);
        [rho_emp,p_perm] = perm_sphere_p(XS(:,1),Y,perm_id,'Pearson');
    
        figure
        f1=scatter(XS(:,1),Y,45);
        f1.MarkerFaceColor=[0.7451,0.7451,0.7451];
        f1.MarkerEdgeColor=[0.7451,0.7451,0.7451];
        set(gca,'LineWidth',3,'FontSize',26,'FontName','Arial');
        %xlabel('PLS1 score ','FontSize',20,'Fontname','Arial');
        %ylabel('SA-NSA t-value in lh','FontSize',20,'Fontname','Arial');
        [p,S] = polyfit(XS(:,1),Y,1);
        [y_fit,delta] = polyval(p,XS(:,1),S);
        hold on;
        plot(XS(:,1),y_fit,'LineWidth',3,'Color','r','LineStyle','-');%拟合曲线
        grid on
        set(gcf,'color','w');
        box off  ;
        xlim([-0.3 0.3])
        set(gca,'XTick',-0.3:0.1:0.3);
        ylim([-6 6])
        set(gca,'YTick',-6:3:6);
        hold on;
        line_horizontal=refline(0,0);%斜率为0，截距为0的参考线，即y=0
        line_horizontal.Color='k';
        line_horizontal.LineStyle=':';
        line_horizontal.LineWidth=2.7;
        hold on
        % plot([0 0], ylim)%x=0轴线
        line_vertical=plot([0 0], ylim,'-k');%x=0轴线设为蓝色
        line_vertical.LineStyle=':';
        line_vertical.LineWidth=2.7;
        %% Bootstrap to get the gene list:自举法(求小样本的统计参数)，目的：求pls1对应的genes weight的std
        %依据自举法求PLS1的weight的统计参数
        geneindex=1:length(X);
        %number of bootstrap iterations:
        bootnum=10000;
        %align PLS components with desired direction for interpretability 
        if R1(1,1)<0  %this is specific to the data shape we were using - will need ammending
            stats.W(:,1)=-1*stats.W(:,1);
            XS(:,1)=-1*XS(:,1);
        end
        
        %依据weight将15043genes降序
        [PLS1w,x1] = sort(stats.W(:,1),'descend');
        PLS1_Symbol=gene_Symbol(x1);%
        geneindex1=geneindex(x1);
        PLS1_X=X(:,x1);%regions*gene_expression:按照weight排序对应symbols的expression
        % 
        % %print out results
        % csvwrite('PLS1_ROIscores.csv',XS(:,1));
        
        %define variables for storing the (ordered) weights from all bootstrap runs
        PLS1weights=[];
        %start bootstrap
        clear res
        for i=1:10000%1-10000
            
            myresample = randsample(size(X,1),size(X,1),1);
            res(i,:)=myresample; %store resampling out of interest
            Xr=X(myresample,:); % define X for resampled subjects
            Yr=Y(myresample,:); % define X for resampled subjects
            [~,~,~,~,~,~,~,stats]=plsregress(Xr,Yr,dim); %perform PLS for resampled data
            clear temp 
            temp=stats.W(:,1);%extract PLS1 weights
            newW=temp(x1); %order the newly obtained weights the same way as initial PLS 
            if corr(PLS1w,newW)<0 % the sign of PLS components is arbitrary - make sure this aligns between runs
                newW=-1*newW;
            end
            PLS1weights=[PLS1weights,newW];%store (ordered) weights from this bootstrap run
            
        end
        
        %get standard deviation of weights from bootstrap runs
        PLS1sw=std(PLS1weights');
        %get bootstrap weights,即Z-scores=每个gene的原始weight/bootstrap_std
        temp1=PLS1w./PLS1sw';
        
        %order bootstrap weights (Z) and names of regions
        [Z1_weight ind1]=sort(temp1,'descend');
        PLS1_Symbol=PLS1_Symbol(ind1);%按照pls1对应于weight的排序：geneSmbol
        geneindex1=geneindex1(ind1);%目前gene在10027里的index
        PLS1_X=PLS1_X(:,ind1);%ROI*gene expression
        XL1=XL(geneindex1,1);%PLS1 loading
        
        clear PLS1_geneWeights;
        PLS1_geneWeights=table(PLS1_Symbol,Z1_weight,XL1,geneindex1','VariableNames',{'geneSymbol','Z_score','PLS1_loading','index'}); 

        yvalues=Y;
        %导入FIQT后的结果
        Z_weight_FIQT=readmatrix(strcat('FIQT_result_',dataset{d},'LesionTrajectory',num2str(c),'_20241223.txt'));%all subjects：划窗

        %选取FDR<0.05的值
        pos=find(Z_weight_FIQT(:,2)>1.96);
        neg=find(Z_weight_FIQT(:,2)<-1.96);
        
        index_pos=Z_weight_FIQT(pos,1);
        index_neg=Z_weight_FIQT(neg,1);
        %找到z-socres中FDR<0.05的z-scores边界值
        limit1=min(PLS1_geneWeights.Z_score(index_pos))
        limit2=max(PLS1_geneWeights.Z_score(index_neg))
        
        figure
        h = histogram(PLS1_geneWeights.Z_score);
        h.Normalization = 'probability';
        h.BinWidth = 0.5;
        h.FaceColor=[95 158 160]./255; 
        h.EdgeColor=[95 158 160]./255; 
        xlim([-10 10]);
        ylim([0 0.15]);
        set(gca,'XTick',-10:5:10);
        set(gca,'YTick',0:0.05:0.15);
        hold on;
        line1=plot([limit1 limit1],ylim,'LineStyle','-','LineWidth',2.7,'Color','k');
        line2=plot([limit2 limit2],ylim,'LineStyle','-','LineWidth',2.7,'Color','k');
        hold on
        line3=plot([3 3],ylim,'LineStyle','-','LineWidth',2.7,'Color',[255 193 193]./255);
        line4=plot([-3 -3],ylim,'LineStyle','-','LineWidth',2.7,'Color',[255 193 193]./255);
        hold on
        line5=plot([4 4],ylim,'LineStyle','-','LineWidth',2.7,'Color',[255 110 180]./255);
        line6=plot([-4 -4],ylim,'LineStyle','-','LineWidth',2.7,'Color',[255 110 180]./255);
        hold on
        line7=plot([5 5],ylim,'LineStyle','-','LineWidth',2.7,'Color',[205 0 0]./255);
        line8=plot([-5 -5],ylim,'LineStyle','-','LineWidth',2.7,'Color',[205 0 0]./255);
        
        set(gca,'LineWidth',2.7,'FontSize',24,'FontName','Times New Roman');
        xlabel('Z-score','FontSize',28,'Fontname','Times New Roman');
        ylabel('Relative frequency','FontSize',28,'Fontname','Times New Roman');
        box off  ;
        set(gcf,'color','w');
        leg=legend([line1,line3,line5,line7],'p_F_D_R=0.05','|Z|=3','|Z|=4','|Z|=5');
        set(leg,'box','off');
        
       

    end
    
end

