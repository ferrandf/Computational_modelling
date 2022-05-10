%Solution of the equations for an AC grain (local problem)
% Explicit method (forward Euler)
% 1) Run the code and analyze the plots
% 2) Solve also the unloading problem: load for 200s, and unload (c_ext=0)
% for 100s

%BAX 1500
epse = 0.37; epsp = 0.8; rhoS = 2.22e6; K = 3.e-3; qm = 0.4; %Porosity constants, Laumuir's isoterm constants.
Dp=[1.e-9, 1.e-8, 1.e-7]; %Difusion intraparticular
B=[0.5, 1.35, 3];  %The constant related to the boundary condition (B.C) mass transfer. (small B, absorb slowly, big B, faster absortion)
R=0.002; %Radius of the bullet.



L=@(q) q./(K*(qm-q)); %Laumir's isoterm.
dL=@(q) qm./(K*(qm-q).^2);


k = 1;
for i = 1:length(Dp)
    for j = 1:length(B)

        c_ext = 1300; %Exterior concentration
        Tfinal = 200;
        dt = 0.01;  
        nOfTimeSteps = round(Tfinal/dt);

        alpha = Dp(i)*B(j)/R;
        beta=35*Dp(i)/(R^2);
        qb=zeros(nOfTimeSteps+1,1);
        qR=qb;
        
        %Loop for time steps
        for n=1:nOfTimeSteps
            cmL=c_ext-L(qR(n));
            qb(n+1)=qb(n)+dt*3*alpha*cmL;
            qR(n+1)=qR(n)+dt*(beta*(qb(n)-qR(n))+10*alpha*cmL);
        end

        ts = [0:nOfTimeSteps]*dt;
        %fig = figure(k);
        %plot(ts,[qR,qb],'-'), legend('qr','qb'),xlabel('t'), ylabel('q'),title("Dp = " + Dp(i) + " B = " + B(j))
        %SaveFileName=strcat('load_Dp_', num2str(Dp(i)), '_B_', num2str(B(j)));
        %saveas(fig,SaveFileName,'png')

        Tfinal = 100;
        dt = 0.01;
        nOfTimeStepsUnloading = round(Tfinal/dt);
        qbLoading=qb; qRLoading=qR;
        qbU = zeros(nOfTimeStepsUnloading+1,1); qRU = qbU;
        qbU(1) = qbLoading(end); qRU(1) = qRLoading(end);
        c_ext = 0;
        
        %Loop for time steps
        for n=1:nOfTimeStepsUnloading
            cmL=c_ext-L(qRU(n));
            qbU(n+1)=qbU(n)+dt*3*alpha*cmL;
            qRU(n+1)=qRU(n)+dt*(beta*(qbU(n)-qRU(n))+10*alpha*cmL);
        end
        
        
        ts = [0:nOfTimeStepsUnloading]*dt;
        %fig = figure(k+1);
        %plot(ts,[qRU,qbU],'-'), legend('qr','qb'),xlabel('t'), ylabel('q'),title("Dp = " + Dp(i) + " B = " + B(j))
        %SaveFileName=strcat('unload_Dp_', num2str(Dp(i)), '_B_', num2str(B(j)));
        %saveas(fig,SaveFileName,'png')
        ts = [0:nOfTimeSteps + nOfTimeStepsUnloading+1]*dt;
        
        fig = figure(k+2);
        qRT = [qRLoading;qRU];
        qbT = [qbLoading;qbU];
        plot(ts,[qRT,qbT],'-','LineWidth',8), legend('qr', 'qb'),xlabel('t'), ylabel('q'), title("Dp = " + Dp(i) + " B = " + B(j))
        
        saveas(fig,strcat('both_Dp_', num2str(Dp(i)), '_B_', num2str(B(j)*100)),'png')
        k = k+3;

    end
end

% sigma = ( ((1-epse)/epse)*3*Dp*B/R )*((1-epsp)*rhoS-epsp*dL(qb));
% figure(2)
% plot(ts(2:end),sigma(2:end)), ylabel('sigma'), xlabel('t')
% 
