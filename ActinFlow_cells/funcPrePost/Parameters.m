% Physical Property parameters - Table values

% parameter.visc = 100;                                         
% parameter.eta = 5;  
% parameter.k = 2.5;                                    
% parameter.chi = 1.75;  
% parameter.diff = 0.15;

parameter.visc      = 10;    
parameter.chi       = 0.05;
parameter.eta       = 0.1;


parameter.kd   = 0.1;  
parameter.kp   = 0.1;
parameter.kdM  = 0.1;  % Set to 0
parameter.kpM  = 0.1;  % Set to 0

parameter.diff  = 0.3;
parameter.diffG = 0.8;
parameter.diffM = 0.4;
parameter.diffm = 0.8;

parameter.k = 0.05; 

buffer = 0.3;

%-----------------------------------------------------------------------
% ALE parameters

parameter.V0 = 0.55;
parameter.coeff1 = 0.12;

%-----------------------------------------------------------------------
% Eta parameters

parameter.eta_xlim = 100;
parameter.eta_xE = linspace(-parameter.eta_xlim,parameter.eta_xlim,1000);
parameter.eta_E = linspace(100,0.1,1000);

%-----------------------------------------------------------------------
% Signaling parameters

parameter.DR = 0.1; 
parameter.DS = 10;  
parameter.time_signal = 20;
parameter.stimulus_orientation = false; %false: Right; true: Left;
parameter.t_chemotaxis = 450;

%-----------------------------------------------------------------------
% Time-space discretizacion parameters

parameter.theta = 1/2;
parameter.Dt = 0.4;
parameter.tEnd = 600; 
parameter.nSteps= parameter.tEnd / parameter.Dt;  
parameter.nDib = 10;
parameter.t = 0:parameter.Dt:parameter.tEnd;

%-------------------------------------------------------------------------
parameter.L0 = 10;  %Initial length
parameter.dom = [-parameter.L0/2,parameter.L0/2];                  %Initial domain
parameter.nx = 21;                       %Number of elements
parameter.Lpre = parameter.L0 * buffer; %Guathier give a buffer area of 50%

%-----------------------------------------------------------------------
% On-line plotting Parameter

parameter.intermplotflag = 1;     % plotting simulation while running the code
parameter.nplotSteps     = 10;    % plot every nplotSteps
parameter.nsavesteps     = 1000; 
parameter.step = 1;

%-----------------------------------------------------------------------
% Polarization Parameter

parameter.zz           = 8;           % charges per molecules
parameter.kbT          = 1.38*10^(-11); % um2 kg / s
parameter.e            = 1.6*10^(-19); % charge [C]
parameter.sign         = 1;        % sign of the charge molecule we want to transport
parameter.EF           = 10^(7);   % Electric field [pN/C]
parameter.Dm           = 0.01;   % effective diffusion [KbT/D--pN.s/um]   
parameter.D            = parameter.kbT/parameter.Dm;

parameter.Velect = parameter.zz * parameter.e * parameter.sign * parameter.EF/parameter.D;


