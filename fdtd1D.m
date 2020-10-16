function fdtd1D()


global eps0 mu0 imp0 c nm fs
load_global_constants();



isWaterFall = true;
% isWaterFall = false;

sourceR.ppw         = 50;
sourceR.omega       = 300e12; % 300/1micron = 300THz
sourceR.timeDelay	= sourceR.ppw;
sourceR.Type        = 'Ricker';
sourceR.isAdd       = true;
% sourceR.isAdd     = false;

sourceG.ppw         = 50;
sourceG.timeDelay	= 4*sourceG.ppw;
sourceG.Type        = 'gauss';
sourceG.isAdd       = true;
% sourceG.isAdd       = false;

source = sourceR;
source = sourceG;

PML = [];
TFSFregion = [];

% FileName = 'init.dat';
% loadDataFromFile(FileName);

MaxTime = 1000; % duration of simulation

additiveSource = 150;
% additiveSource = 3;

%% initialization
Objects = [];
domainSize = 200; % x size of domain
% FDTDregion = [0,400]*nm;
FDTDregion = [-200,200]*nm;
PML = struct('layers',10, 'm',3.5, 'coef',0.8);
calFDTD = FDTD(domainSize,FDTDregion,PML);       % creat 1D FDTD grid for TM[Ez,Hy] (add boolean variable isTE for TE[])
calFDTD.initGrid(Objects,[]);                    % initialize the 1D FDTD grid
if isempty(PML)
    calFDTD.initABC('2abc');                % initialize ABC
end

% TFSFregion = [10;domainSize-10];
if ~isempty(TFSFregion)
    calFDTD.initTFSF(TFSFregion,source);        % initialize TFSF
end

% graphic output initialization
x  = calFDTD.getGridXYZ('Ez')/nm; % in nm
dt = calFDTD.get_dt()/fs; % in fs
figure;
if isWaterFall
    waterfall(x,zeros(size(x)),calFDTD.snapshotE('z').');
else
    plot3(x,zeros(size(x)),calFDTD.snapshotE('z').','color','blue');
end
hf = gca; hold on
view(-20,62)
xlabel('x, nm')
ylabel('Time, fs')
%zlabel('Electric field')


% do time stepping
X = calFDTD.getGridXYZ('Ex');
for n_time = 1:MaxTime
    calFDTD.updateH();          % update magnetic field
    if isempty(TFSFregion)
        calFDTD.add_Hy_inc(additiveSource-1, n_time-1, source); % add a source
    else
        calFDTD.updateTFSF(n_time);	% apply TFSF boundary
    end
    calFDTD.updateE();          % update electric field
    if isempty(TFSFregion)
        calFDTD.add_Ez_inc(additiveSource, n_time, source); % add a source
    end
    if isempty(PML)
        calFDTD.applyABC();         % apply ABC
    end
    
    %% Energy in the simulation volume
    if trapz(X, abs(calFDTD.Poynting(X)))<1e-30
        break
    end

    if mod(n_time,1)==0
        if isWaterFall
            waterfall(hf,x,n_time*dt,calFDTD.snapshotE().');
        else
            plot3(hf,x,n_time*ones(size(x))*dt,calFDTD.snapshotE('z').','color','blue');
        end
    end
    
end % of time-stepping
 
end

