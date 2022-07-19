function [mdl]= mymultisvm2(featuretrain1,featuretrain2,featuretrain3)
%ONE VS ONE
% 1 versus 2  
options.MaxIter = 1000000;
data1_a= featuretrain1;
data2_a= featuretrain2;

datatrain_a=[data1_a,data2_a];
labeltrain_a = [ones(1,size(data1_a,2)),2*ones(1,size(data2_a,2))];
mdl.svm1=svmtrain(datatrain_a',labeltrain_a, 'kernel_function','rbf', 'Options', options);

% 2 versus 3 
data1_b= featuretrain2;
data2_b= featuretrain3;
datatrain_b=[data1_b,data2_b];
labeltrain_b= [ones(1,size(data1_b,2)),2*ones(1,size(data2_b,2))];
mdl.svm2= svmtrain(datatrain_b',labeltrain_b, 'kernel_function','rbf','Options', options);

% 3 versus 1  
data1_c= featuretrain3;
data2_c= featuretrain1;
datatrain_c=[data1_c,data2_c];
labeltrain_c = [ones(1,size(data1_c,2)),2*ones(1,size(data2_c,2))];
mdl.svm3= svmtrain(datatrain_c',labeltrain_c, 'kernel_function','rbf','Options', options);
end

