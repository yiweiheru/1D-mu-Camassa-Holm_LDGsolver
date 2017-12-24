clear
close all

Tfinal =15;
% figure_at_time=[0,Tfinal];
% figure_at_time=[0,1,5,10];
figure_at_time=[0,1,3,5,10,15];
ord_num = 3;
ir_num = 5;

n_RK  = 4;
period = 1;
CS = 1;	% indicator of the initial data

P0 = zeros(CS,1);
Q0 = zeros(CS,1);
switch CS
    case 1
        P0(1) = 0.1;        Q0(1) = -0.5;
    case 2
        P0(1) = 0.1;      Q0(1) = 0.2;
        P0(2) = 0.08;     Q0(2) = 0.1;
    case 3
        P0(1) = 0.1;      Q0(1) = 0.2;
        P0(2) = 0.08;     Q0(2) = 0.1;
        P0(3) = 0.12;     Q0(3) = 0.05;
end

UStore = zeros((ord_num+1)*(10*2^(ir_num))+1,ir_num,ord_num);
UexcStore = zeros((ord_num+1)*(10*2^(ir_num))+1,ir_num,ord_num);

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
        
        cfl = 0.1;
        dt = cfl * dx;
        Tsteps = floor((Tfinal-0.1*dt)/dt)+1;
        dt_final = Tfinal - (Tsteps-1) * dt;
        
        U0 = setInitial(Nelm,elm_size,x,CS,period,P0,Q0);
        
        [ Amat,massMat_inv ] = getAmat( Ord,Nelm,x );
        U = U0;
        
        Time = 0;
        for nt = 1:Tsteps
            
            if Time == 0
                    fig = fig+1;
                    figure(fig)
%                     plot_uh( U,Ord,Nelm,x ,"numerical")
%                     hold on
                    [ Uexc,~,~ ] = multi_pkns_solu( P0,Q0,Nelm,elm_size, x ,period,Time,CS );
                    plot_uh( Uexc,Ord,Nelm,x ,"exact")
                    grid on
                    xlabel('x')
                    ylabel('u')
%                     Title_str = strcat('P^',num2str(Ord),'polynomial,t = ',num2str(Time));
                    Title_str = strcat('t=',num2str(Time));
                    title(Title_str);
            end
            
            if nt == Tsteps-1
                dt = dt_final;
            end
            
            U = RKn( Ord,x,Nelm,U,Amat,massMat_inv,n_RK,dt,Time );
            Time = Time+dt;
            
            for p=1:size(figure_at_time,2)
                if figure_at_time(p)-dt/2 < Time && Time <= figure_at_time(p)+dt/2
                    fig = fig+1;
                    figure(fig)
                    plot_uh( U,Ord,Nelm,x ,"numerical")
                    hold on
                    [ Uexc,~,~ ] = multi_pkns_solu( P0,Q0,Nelm,elm_size, x ,period,Time,CS );
                    plot_uh( Uexc,Ord,Nelm,x ,"exact")
                    grid on
                    xlabel('x')
                    ylabel('u')
%                     Title_str = strcat('P^',num2str(Ord),'polynomial,T = ',num2str(Time));
                    Title_str = strcat('t=',num2str(Time));
                    title(Title_str);
                    legend('LDG','Exact')
                end
            end
        end
    end
end





