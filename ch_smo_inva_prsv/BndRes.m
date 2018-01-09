function [ res1,res2 ] = BndRes( Ord,Nelm,Time,period )

elm_size=Ord+1;

uf=zeros(elm_size,2);
uf(:,1)=basis_1d(Ord,-1);
uf(:,2)=basis_1d(Ord,1);

res1=zeros(Nelm*elm_size,1);
res2=zeros(Nelm*elm_size,1);

XL=-period/2;
XR=period/2;

uLp=0.25*exp(XL-0.25*Time);
uRp=0.25*exp(-1*(XR-0.25*Time));
rLm=0.25*exp(XL-0.25*Time);
rRm=-0.25*exp(-1*(XR-0.25*Time));       
        
for i=1:elm_size
    res1(i,1)=rLm*uf(i,1);
    res1((Nelm-1)*elm_size+i,1)=rRm*uf(i,2);
    
    res2(i,1)=uLp*uf(i,1);
    res2((Nelm-1)*elm_size+i,1)=uRp*uf(i,2);
end

end

