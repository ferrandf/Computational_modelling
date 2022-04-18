function [] = plot_time_signal(parameter,cellinfo)

X = cellinfo{1}.meshparam.X;
indR = cellinfo{1}.meshparam.indR;
indS = cellinfo{1}.meshparam.indS;
n = parameter.step;
Dt = parameter.Dt;
tEnd = parameter.tEnd / 2.;
sol = cellinfo{1}.meshparam.vec_RS;
stimulus = cellinfo{1}.meshparam.stimulus;


time = Dt*n;
R = sol(indR);
S = sol(indS);

if parameter.chemotaxis
    if (time == Dt)
        figure(30);
        subplot(1,3,1);
        plot(X,sol(indR),X,sol(indS),X,stimulus,'--','LineWidth',2)
        xlabel('x [\mu m]')
        lg = legend('R','S','stimulus');
        xlim([X(1) X(end)]);
        ylim([0 2]);
        set(lg,'FontSize',16,'Location','NorthEast')
        set(gca,'FontSize',16)
        title(['t = ',num2str(time-Dt,'%d'), ' s'])
    elseif (time==parameter.t_chemotaxis)
        figure(30);
        subplot(1,3,2);
        plot(X,R,X,S,X,stimulus,'--','LineWidth',2)
        xlabel('x [\mu m]')
        lg = legend('R','S','stimulus');
        xlim([X(1) X(end)]);
        ylim([0 2]);
        set(lg,'FontSize',16,'Location','NorthEast')
        set(gca,'FontSize',16)
        title(['t = ',num2str(time,'%d'), ' s'])
    elseif (time == tEnd*2.)
        figure(30);
        subplot(1,3,3);
        plot(X,R,X,S,X,stimulus,'--','LineWidth',2)
        xlabel('x [\mu m]')
        lg = legend('R','S','stimulus');
        xlim([X(1) X(end)]);
        ylim([0 2]);
        set(lg,'FontSize',16,'Location','NorthEast')
        set(gca,'FontSize',16)
        title(['t = ',num2str(time,'%d'), ' s'])
    end
else
    if (time == Dt)
        figure(30);
        subplot(1,3,1);
        plot(X,sol(indR),X,sol(indS),X,stimulus,'--','LineWidth',2)
        xlabel('x [\mu m]')
        lg = legend('R','S','stimulus');
        xlim([X(1) X(end)]);
        ylim([0 2]);
        set(lg,'FontSize',16,'Location','NorthEast')
        set(gca,'FontSize',16)
        title(['t = ',num2str(time-Dt,'%d'), ' s'])
    elseif (time == parameter.init_signal+parameter.time_signal)
        figure(30);
        subplot(1,3,2);
        plot(X,R,X,S,X,stimulus,'--','LineWidth',2)
        xlabel('x [\mu m]')
        lg = legend('R','S','stimulus');
        xlim([X(1) X(end)]);
        ylim([0 2]);
        set(lg,'FontSize',16,'Location','NorthEast')
        set(gca,'FontSize',16)
        title(['t = ',num2str(time,'%d'), ' s'])
    elseif (time == tEnd)
        figure(30);
        subplot(1,3,3);
        plot(X,R,X,S,X,stimulus,'--','LineWidth',2)
        xlabel('x [\mu m]')
        lg = legend('R','S','stimulus');
        xlim([X(1) X(end)]);
        ylim([0 2]);
        set(lg,'FontSize',16,'Location','NorthEast')
        set(gca,'FontSize',16)
        title(['t = ',num2str(time,'%d'), ' s'])
    end
    
end

end