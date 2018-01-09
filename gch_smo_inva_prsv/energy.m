
n_RK=4;
Tfinal=10;
ord_num=2;
ir_num=4;
cfl=0.1;
% 1 is Alternating; 2 is Central. Here
global type_of_flux
type_of_flux = 1;

for Ord=ord_num:ord_num
    for ir=ir_num:ir_num        
        Nelm=10*2^(ir-1);
        dx=period/Nelm;

        dt=dx*cfl;
        Tsteps=floor((Tfinal-0.000000001)/dt)+1;
        
        x=-period/2:dx:period/2;
        elm_size=Ord+1;
        U0 = setInitial(Nelm,elm_size,x,Xexc,uexc);
        [ Amat,Dumat,massMat,massMat_inv ] = getAmat( Ord,Nelm,x );
        U=U0;
        nn = 0;
        Time=0;

        for nt=1:Tsteps
            if nt == Tsteps
               dt = Tfinal - Time;
            end
               
            U = RKn( x,Ord,Nelm,U,Amat,massMat_inv,Dumat,massMat,n_RK,dt );
            R = Dumat*U;
            Ey(nt) = getEnergy( U,R,Ord,Nelm,x); 
            
            
            Time=Time+dt;
        end

  
    end
end
hold on
Ey_pt = 0.1*linspace(0,Tfinal,size(Ey,2));
plot(Ey_pt,Ey,'m')

% Ey_pt = 0:0.1:Tfinal;
% NN = size(Ey_pt,2);
% for nn=1:NN
%     Time=Ey_pt(nn);
%     [ U_final,R_final ] =getFinalExactSol( Nelm,elm_size,x,Time,Xexc,uexc,rexc ,period,c);
%     Ey_exc(nn)=getEnergy( U_final,R_final,Ord,Nelm,x); 
% end
% hold on 
% plot(Ey_pt,Ey_exc,'y')

