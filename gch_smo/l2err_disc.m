function [ L2err ] = L2_error( Uh,Time,Ord,Nelm,x,Xexc,uexc,rexc,period,c )

elm_size=Ord+1;
[Ue,~]= getFinalExactSol( Nelm,elm_size,x,Time,Xexc,uexc,rexc,period,c );
U_diff=Ue-Uh;
u_diff=uhTransform( Nelm,elm_size,U_diff );

  NN = 20;
  xi = linspace(-1,1,NN+1);
  un_at_xi = zeros(elm_size,NN+1);
  for k = 1 : NN+1
    un_at_xi(:,k) = basis_1d(Ord,xi(k));
  end

  l2err = 0;
  for ne = 1 : Nelm
    for k = 1 : NN+1
      l2err = l2err + (u_diff(ne,:)*un_at_xi(:,k))^2;
    end
  end
	l2err = sqrt(l2err/((NN+1)*Nelm));
end
