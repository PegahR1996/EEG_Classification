function [w]= mymulticlasscsp(nch,gg,tr1,tr2,tr3,datatrain1,datatrain2,datatrain3)
%on vs all
%% 1 versus 2 3 4
data1= datatrain1;
data2= [datatrain2,datatrain3];
[w1] = mycsp(nch,gg,tr1,tr2+tr3,data1,data2);
%% 2 versus 1 3 4
data1= datatrain2;
data2= [datatrain1,datatrain3];
[w2] = mycsp(nch,gg,tr2,tr1+tr3,data1,data2);
%% 3 versus 1 2 4
data1= datatrain3;
data2= [datatrain1,datatrain2];
[w3] = mycsp(nch,gg,tr3,tr1+tr2,data1,data2);


w= cat(3,w1,w2,w3);

end