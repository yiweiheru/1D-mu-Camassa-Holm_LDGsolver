function [ L2err ] = L2_error( Uh,Time,Ord,Nelm,x )

elm_size=Ord+1;
[Ue,~]= getFinalExactSol( Nelm,elm_size,x,Time );
U=Ue-Uh;
u=uhTransform( Nelm,elm_size,U );

npt_quad=Ord+2;
[qpt, qwt] = QuadLG(npt_quad);

un=zeros(elm_size,npt_quad);
for k = 1 : npt_quad
    un(:,k)=basis_1d(Ord,qpt(k));
end

Val=0;
for ne=1:Nelm
    if x(ne)>=15
    Jaco=(x(ne+1)-x(ne))/2;
    for ik=1:npt_quad
        Val=Val+(u(ne,:)*un(:,ik))^2*Jaco*qwt(ik);
    end
    end
end

L2err=sqrt(Val);

end

