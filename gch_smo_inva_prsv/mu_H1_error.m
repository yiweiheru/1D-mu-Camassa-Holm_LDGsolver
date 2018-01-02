function [ muH1_err ] = mu_H1_error( Uh,Rh,Time,Ord,Nelm,x,Xexc,uexc,rexc ,period,c)

elm_size=Ord+1;
[Ue,Re]= getFinalExactSol( Nelm,elm_size,x,Time,Xexc,uexc,rexc,period,c );
U=Ue-Uh;
R=Re-Rh;
u=uhTransform( Nelm,elm_size,U );
r=uhTransform( Nelm,elm_size,R );

npt_quad=Ord+2;
[qpt, qwt] = QuadLG(npt_quad);

un=zeros(elm_size,npt_quad);
for k = 1 : npt_quad
    un(:,k)=basis_1d(Ord,qpt(k));
end

Val=0;
for ne=1:Nelm
    Jaco=(x(ne+1)-x(ne))/2;
    for ik=1:npt_quad
        Val=Val+(r(ne,:)*un(:,ik))^2*Jaco*qwt(ik);
    end
end


mu_un=zeros(elm_size,1);
for ik=1:npt_quad
    mu_un(:,1)=mu_un(:,1)+qwt(ik)*un(:,ik);
end
mu_u=0;
for ne=1:Nelm
    Jaco=(x(ne+1)-x(ne))/2;
    for j=1:elm_size
        mu_u=mu_u+mu_un(j,1)*u(ne,j)*Jaco;
    end
end

muH1_err=sqrt(Val+mu_u^2);
% muH1_err=sqrt(Val);

end

