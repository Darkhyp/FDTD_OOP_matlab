function fdtd1D()

global eps0 mu0 c nm fs
eps0 = 8.8541878176209821e-12;
mu0  = 1.2566370614359173e-06;
c    = 299792458; % 1/sqrt(eps0*mu0)
nm   = 1e-9;
fs   = 1e-15;


isWaterFall = true;
% isWaterFall = false;

source.ppw          = 20;
source.timeDelay	= source.ppw;
% sourceType        = 'gauss';
source.Type         = 'Ricker';
source.isAdd        = true;
% source.isAdd        = false;

sourceG.ppw          = 20;
sourceG.timeDelay	= 4*source.ppw;
sourceG.Type        = 'gauss';


% FileName = 'init.dat';
% loadDataFromFile(FileName);

MaxTime = 1000; % duration of simulation

additiveSource = 50;
additiveSource = 3;

%% initialization
SizeX = 200; % x size of domain
FDTDregion = [0,400*nm];
calFDTD = FDTD(SizeX,FDTDregion);           % creat 1D FDTD grid (isTE for 2D)
calFDTD.initGrid([],[]);                    % initialize the 1D FDTD grid
calFDTD.initABC('2abc');                    % initialize ABC
calFDTD.initTFSF(additiveSource,source);	% initialize TFSF

% graphic output initialization
x  = calFDTD.getGridXYZ('Ez')/nm; % in nm
dt = calFDTD.get_dt()/fs; % in fs
figure;
if isWaterFall
    waterfall(x,0,calFDTD.snapshotE('z').');
else
    plot3(x,0*ones(1,SizeX),calFDTD.snapshotE().','color','blue');
end
hf = gca; hold on
view(-20,62)
xlabel('x, nm')
ylabel('Time, fs')
%zlabel('Electric field')


% do time stepping
for n_time = 1:MaxTime
    calFDTD.updateH();          % update magnetic field
    calFDTD.updateTFSF(n_time);	% apply TFSF boundary
%     calFDTD.add_Hy_inc(additiveSource-1, n_time, source); % add a source
    calFDTD.updateE();          % update electric field
%     calFDTD.add_Ez_inc(additiveSource, n_time, source); % add a source
    calFDTD.applyABC();         % apply ABC

    if mod(n_time,10)==0
        if isWaterFall
            waterfall(hf,x,n_time*dt,calFDTD.snapshotE().');
        else
            plot3(hf,x,n_time*ones(1,SizeX)*dt,calFDTD.snapshotE('z').','color','blue');
        end
    end
    
end % of time-stepping
 
end

