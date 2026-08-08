clc;
clear;

Nperm=10000;
nSubtype=2;
dataset={'discovery','replication1','replication2'};
Clinical={'Anxiety','Weight','Cognition','Hysteresis','Dysomnia'}';
path='your URL';
clear RESULT R  pvalue;
for d=1:length(dataset)
    for c=1:nSubtype
        clear DEMO FC fl;
        %DEMO与FC一一对应
        DEMO=readtable(strcat(path,dataset{d},'\Your demographics',num2str(c),'.xlsx'));
       
        switch c
            case 1
                fl=DEMO.P1byWT(:);
            case 2
                fl=DEMO.P2byWT(:);
        end

        %% Regress out age, sex, motion and sites from factor loading
        conf=[DEMO.age,DEMO.sex,DEMO.education];
        fl_d = fl-conf*(pinv(conf)*fl); % Regress out regressors from factor loading
        %% Regress out age, sex, motion and sites from behavioral scores
        scores=[DEMO.Anxiety,DEMO.Weight,DEMO.Cognition,DEMO.Hysteresis,DEMO.Dysomnia];
        scores_d = zeros(size(scores));
        for i=1:size(scores,2)
            % Regress out regressors from behavioral scores
            scores_d(:,i) = scores(:,i)-conf*(pinv(conf)*scores(:,i)); 
        end
        %% CCA
        % Age, sex, motion & sites regressed on both sides
        [A, B, R, U, V, ~] = canoncorr(scores_d,fl_d); 
        %% Permutation test
        Rp=zeros(Nperm,size(fl,2)); 
        clear pVal PAPset;
        for j=1:Nperm
            PAPset(:,j)=randperm(length(fl),length(fl));
            %fprintf('Permutation No.%d\n',j);
            % Age, sex, motion & sites regressed on both sides
            [Ap,Bp,Rp(j,:),Up,Vp] = canoncorr(scores_d,fl_d(PAPset(:,j),:)); 
        end
        
        pVal = (1+sum(Rp(2:end,1)>=R))/Nperm;
        fprintf('r is %f ; p value is %f\n',R,pVal);
        %Ncca = sum(pVal<0.05)

        % Always keep B (canonical coefficient for factor loading) positive 
        if B < 0
            A = -A;
            U = -U;
            B = -B;
            V = -V;
        end
        
        % Compute structure coefficient for behavioral scores
        strucCorr_score(:,1) = corr(U(:,1),scores,'rows','complete'); 
        
    end%for c=1:ncluster
end%for d=1:length(dataset)
