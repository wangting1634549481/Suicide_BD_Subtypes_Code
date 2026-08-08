
function [tstat,pval,regional_result_beforeFDR,regional_result_FDR,sigregs]=regionaltmap(convaribales,Group,y,region_name)

    v1=convaribales.v1;
    v2=convaribales.v2;
    v3=convaribales.v3;
    v4=Group;
    nregs=size(y,2);

    clear tstat pval;
    for region=1:nregs

       
        tbl = table(v1,v2,v3,v4,y(:,region),...
        'VariableNames',{'age','sex','education','group','regional'});
        tbl.sex = categorical(tbl.sex);
        tbl.group = categorical(tbl.group);
        lm = fitlm(tbl,'regional~age+sex+education+group');
        tstat(region)=lm.Coefficients{5,3};
        pval(region)=lm.Coefficients{5,4};
        
        

    end

    
    id=find(pval<0.05);
    if ~isempty(id)
        regional_result_beforeFDR=table(region_name(id,:),tstat(id)',pval(id)',...
                'VariableNames',{'regions','t_statistics','p_value'});
    else
        regional_result_beforeFDR=NaN;
    end


    pvalue_fdr = mafdr(pval,'BHFDR',1); % FDR corrected p-values
    sigregs=find(pvalue_fdr<0.05); % list of the statistically significant regions
    if ~isempty(sigregs)
        regional_result_FDR=table(region_name(sigregs,:),tstat(sigregs)',pvalue_fdr(sigregs)',...
            'VariableNames',{'regions','t_statistics','adjusted_p_value'});
    else
        regional_result_FDR=NaN;
    end


end






