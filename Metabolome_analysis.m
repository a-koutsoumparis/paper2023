clc;
close all;
clear;

% Read in data
% Adjust indexes according to the number of samples
a = readtable('FILE NAME');
Intensity = table2array(a(:,3:14));
Proteins = table2array(a(:,1:2));

% Quantile normalization
n_Int = quantilenorm(Intensity,'Median','true');

% Remove high variance metabolites
% Calculate coefficient of vaiation
cv1 = std(log2(n_Int(:,1:6)),[],2)./mean(log2(n_Int(:,1:6)),2);
cv2 = std(log2(n_Int(:,7:12)),[],2)./mean(log2(n_Int(:,7:12)),2);

% Accept CVs < 97.5% of max
accepted_Int = n_Int(cv1<quantile(cv1,0.975)&cv2<quantile(cv2,0.975),:);
accepted_proteins = Proteins(cv1<quantile(cv1,0.975)&cv2<quantile(cv2,0.975),:);

% Calculate statistical parameters
% Calculate the mean value
m_wt = mean(accepted_Int(:,1:6),2);
m_m = mean(accepted_Int(:,7:12),2);

% Fold-change
log2FC = log2(m_m./m_wt);

% Significance assessment
for i=1:size(accepted_Int,1)
[h(i),p_value(i)] = ttest2(accepted_Int(i,1:6),accepted_Int(i,7:12),'tail','both','vartype','unequal');
end
p_value = p_value';

% False discovery rate
FDR = mafdr(p_value,'BHFDR',true);