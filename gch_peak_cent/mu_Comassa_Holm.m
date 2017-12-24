clear
close all

Tfinal =500;
figure_at_time=[10,50,100,250,500];
% figure_at_time=[1,2,3];

ord_num = 3;
ir_num = 4;

n_RK  = 4;
period = 1;
CS = 1;	% indicator of the initial data

P0 = zeros(CS,1);
Q0 = zeros(CS,1);
switch CS
    case 1
        P0(1) = 0.3;        Q0(1) = -0.5;  %P0 is speed of wave; Q0 is location of peak.
    case 2
        P0(1) = 0.1;      Q0(1) = 0.2;
        P0(2) = 0.08;     Q0(2) = 0.1;
    case 3
        P0(1) = 0.1;      Q0(1) = 0.2;
        P0(2) = 0.08;     Q0(2) = 0.1;
        P0(3) = 0.12;     Q0(3) = 0.05;
end

% U1Store = zeros(size(figure_at_time),(ord_num+1)*(10*2^(ir_num-1)));
% U2Store = zeros(size(figure_at_time),(ord_num+1)*(10*2^(ir_num-1)));
% UexcStore = zeros(size(figure_at_time),(ord_num+1)*(10*2^(ir_num-1)));


fig = 0;
for Ord = ord_num:ord_num
    for ir = ir_num:ir_num
        
        Nelm = 10*2^(ir-1);
        dx = period/Nelm;
        if CS ==1
            x = -period/2:dx:period/2;
        elseif CS == 2 || CS == 3
            x = 0:dx:period;
        end
        elm_size = Ord+1;
        
        cfl = 0.05;
        dt = cfl * dx;
        Tsteps = floor((Tfinal-0.1*dt)/dt)+1;
        dt_final = Tfinal - (Tsteps-1) * dt;
        
        U0 = setInitial(Nelm,elm_size,x,CS,period,P0,Q0);
        
        [ Amat,massMat_inv ] = getAmat( Ord,Nelm,x );
        UC = U0;
        UA = U0;
        

        
        
        
        Time = 0;
        for nt = 1:Tsteps
            
%             if Time == 0
%                     fig = fig+1;
%                     figure(fig)
%                     [ Uexc,~,~ ] = multi_pkns_solu( P0,Q0,Nelm,elm_size, x ,period,Time,CS );
%                     plot_uh( Uexc,Ord,Nelm,x ,"exact")
%                     grid on
%                     xlabel('x')
%                     ylabel('u')
%                     Title_str = strcat('t=',num2str(Time));
%                     title(Title_str);
%             end
%             
            if nt == Tsteps-1
                dt = dt_final;
            end
            
            UC = RKn_C( Ord,x,Nelm,UC,Amat,massMat_inv,n_RK,dt,Time );
%             UA = RKn_A( Ord,x,Nelm,UA,Amat,massMat_inv,n_RK,dt,Time );
            Time = Time+dt;
            
            for p=1:size(figure_at_time,2)
                if figure_at_time(p)-dt/2 < Time && Time <= figure_at_time(p)+dt/2
                    fig = fig+1;
                    figure(fig)
                    plot_uh( UC,Ord,Nelm,x ,"Cnumerical")
                    hold on
%                     plot_uh( UA,Ord,Nelm,x ,"Anumerical")
                    [ Uexc,~,~ ] = multi_pkns_solu( P0,Q0,Nelm,elm_size, x ,period,Time,CS );
                    plot_uh( Uexc,Ord,Nelm,x ,"exact")
                    grid on
                    xlabel('x')
                    ylabel('u')
                    Title_str = strcat('t=',num2str(Time));
                    title(Title_str);
%                     legend('Central','Alternating','Exact')
                end
            end
        end
    end
end





