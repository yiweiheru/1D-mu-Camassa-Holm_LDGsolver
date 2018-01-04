function [ Residue4 ] = residue4( x,Nelm,Ord,qh,rh,mh )

elm_size=Ord+1;

npt_quad=Ord+2;
[qpt, qwt] = QuadLG(npt_quad);

un=zeros(elm_size,npt_quad);
un_der=zeros(elm_size,npt_quad);
uf=zeros(elm_size,2);

for k = 1 : npt_quad
    un(:,k)      =basis_1d(Ord,qpt(k));
    un_der(:,k)=basisDer_1d(Ord,qpt(k));
end

uf(:,1)=basis_1d(Ord,-1); %left boundary of elment
uf(:,2)=basis_1d(Ord,1);  %right boundary of element

Residue4=zeros(Nelm*elm_size,1);
Intergral=zeros(Nelm*elm_size,1);
Flux_L=zeros(Nelm*elm_size,1);
Flux_R=zeros(Nelm*elm_size,1);

%imply the periodic boundary condition to get the neighbourhood.
nei=zeros(Nelm,2);
for ne=1:Nelm
    if ne==1
        nei(ne,1)=Nelm;
        nei(ne,2)=2;
    elseif ne==Nelm
        nei(ne,1)=Nelm-1;
        nei(ne,2)=1;
    else
        nei(ne,1)=ne-1;
        nei(ne,2)=ne+1;
    end
end

for ne=1:Nelm
    Jaco=(x(ne+1)-x(ne))/2;
    for ik=1:npt_quad
        q=0;m=0;r=0;
        for j=1:elm_size
            q=q+un(j,ik)*qh(ne,j);
            m=m+un(j,ik)*mh(ne,j);
            r=r+un(j,ik)*rh(ne,j);
        end

        for i=1:elm_size
            num=(ne-1)*elm_size+i;
            Intergral(num)=Intergral(num) - qwt(ik)*un_der(i,ik)*q ...
                + 2*Jaco*qwt(ik)*un(i,ik)*(r*m);
        end
    end

    qLm=0;qLp=0;qRm=0;qRp=0;
    for j=1:elm_size
        qLm=qLm+uf(j,2)*qh(nei(ne,1),j);
        qLp=qLp+uf(j,1)*qh(ne,j);
        qRm=qRm+uf(j,2)*qh(ne,j);
        qRp=qRp+uf(j,1)*qh(nei(ne,2),j);
    end

    qhat_L=1/2*(qLm+qLp);
    qhat_R=1/2*(qRm+qRp);

    for i=1:elm_size
        num=(ne-1)*elm_size+i;
        Flux_L(num)=Flux_L(num)+uf(i,1)*qhat_L;
        Flux_R(num)=Flux_R(num)+uf(i,2)*qhat_R;
    end
end

Residue4=Residue4-(Intergral+Flux_R-Flux_L);

end

