% 1D FDTD simulation with an additive source.

SIZE = 200;
LOSS = 0.01; % = sigma*delta_t/2/epsilon, 1D FDTD simulation of a lossy dielectric region
LOSS = 0; % = sigma*delta_t/2/epsilon, 1D FDTD simulation of a lossy dielectric region
imp0 = 377.0; % = delta_t/delta_x/epsilon, impedance
maxTime = 700;
LOSS_LAYER = 180;
LOSS_LAYER = SIZE;

delta_x = 1;
delta_t = 1;
tau_p = 10;
tau_p = 20;
t0 = 4*tau_p;

additiveSource = 50;

%% initialize electric field
ez = zeros(1,SIZE);
% initialize magnetic field
%hy = zeros(1,SIZE);
hy = zeros(1,SIZE-1);

% set relative permittivity
epsR = 1^2*ones(1,SIZE);
epsR(1,100:end) = 3^2;
% or set electric-field update coefficients
ceze = ones(1,SIZE);
cezh = imp0*ones(1,SIZE);
ceze(1,100:LOSS_LAYER) = 1;
cezh(1,100:LOSS_LAYER) = imp0/3^2;
% lossy dielectric
ceze(1,LOSS_LAYER:end) = (1-LOSS)/(1+LOSS);
cezh(1,LOSS_LAYER:end) = imp0/3^2/(1+LOSS);


%filename = 'FDTD out.dat';
%FileOut = fopen(filename, 'w');


figure
waterfall(1:SIZE,0,ez)
hf = gca;
hold on
view(-30,76)
xlabel('Space [spatial index]')
ylabel('Time [frame number]')
%zlabel('Electric field')

% do time stepping
mmH = 1:SIZE-1;
%mmE = 2:SIZE;
mmE = 2:SIZE-1;
for Time =0:maxTime-1
    %% update magnetic field
    % simple ABC for hy(1) and hy(SIZE)
%	hy(1) = hy(2); 
% 	hy(SIZE) = hy(SIZE-1);
    hy(1,mmH) = hy(1,mmH) + (ez(1,mmH+1) - ez(1,mmH)) / imp0;
%	hy(1,mmH) += (ez(1,mmH+1) - ez(1,mmH)) / imp0;
    % correction for Hy adjacent to TFSF boundary
    hy(additiveSource-1) = hy(additiveSource-1) - exp(-((Time-t0)/tau_p)^2) / imp0;
%	hy(additiveSource-1) -= exp(-((Time-t0)/tau_p)^2) / imp0;
    ez(additiveSource) = ez(additiveSource) + exp(-((Time+1-t0)/tau_p)^2);
% 	ez(additiveSource) += exp(-((Time-t0)/tau_p)^2 );


    %% update electric field
    % simple ABC for ez(1) and ez(SIZE)
    ez(1) = ez(2); 
% 	ez(SIZE) = ez(SIZE-1);
%	ez(1,mmE) += (hy(1,mmE) - hy(1,mmE-1)) * imp0 ./ epsR(mmE);
    ez(1,mmE) = ceze(1,mmE).*ez(1,mmE) + cezh(1,mmE).*(hy(1,mmE) - hy(1,mmE-1));
  
%     % use additive source at node $additiveSource
%     ez(additiveSource) = ez(additiveSource) + exp(-((Time+1-t0)/tau_p)^2);
% % 	ez(additiveSource) += exp(-((Time-t0)/tau_p)^2 );
  
    % show/write snapshot if time a multiple of 10
    if (mod(Time,10) == 0)
%       fprintf(FileOut,'%s\n', mat2str(ez));
      waterfall(hf,1:SIZE,Time,ez);
    end
%     plot(1:SIZE,ez,'.-');
end %end of time-stepping

%fclose(FileOut);
