function [cspdata] =myapplaycsp(f,w,data)

for k = 0:((size(data,2)/f)-1)
 jk = f*(k+1);
    
 tjj=data(:,k*f+1:jk); 
      s1=w'*tjj;
    d = f*(k+1);
  cspdata(:,k*f+1:d)= s1;
end
end