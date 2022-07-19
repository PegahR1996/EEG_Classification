function [w] =mycsp(nch,ff,tr1,tr2,dataa,datab)
  m=1;
 Rh=0;

for k1 = 0:tr1-1
 jk1 = ff*(k1+1);
    
 t1=dataa(:,k1*ff+1:jk1); 
    for i1= 1:nch
       q1(i1,:)=( t1(i1,:) - mean(t1(i1,:)));
    end
    rh=(q1*q1')/(trace(q1*q1'));
    Rh=Rh+rh;
end
Rh=Rh/tr1;

Rf=0;
for k = 0:tr2-1
 jk = ff*(k+1);
    
 t=datab(:,k*ff+1:jk); 
    for i= 1:nch
       q(i,:)= (t(i,:) - mean(t(i,:)));
    end
    rf=(q*q')/(trace(q*q'));
    Rf=Rf+rf;
end
Rf=Rf/tr2;

[u,v]=eig(Rh,Rf);
v=diag(v);
[v,id]=sort(v,'descend');
u=real(u(:,id));
w=real([u(:,1:m),u(:,end-m+1:end)]);

% end
end


