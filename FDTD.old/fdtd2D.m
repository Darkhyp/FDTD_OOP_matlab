function fdtd2D()
%#codegen

global eps0 mu0 c nm fs
eps0 = 8.8541878176209821e-12;
mu0  = 1.2566370614359173e-06;
c    = 299792458; % 1/sqrt(eps0*mu0)
nm   = 1e-9;
fs   = 1e-15;

isTE = true;
% isTE = false;

source.ppw          = 20;
source.timeDelay	= source.ppw;
% sourceType        = 'gauss';
source.Type         = 'Ricker';
source.isAdd        = false;

sourceG.ppw          = 20;
sourceG.timeDelay	= 3*source.ppw;
sourceG.Type        = 'gauss';
sourceG.isAdd        = true;

% FileName = 'init.dat';
% loadDataFromFile(FileName);

MaxTime = 1000; % duration of simulation

% objects
Objects = {struct('type','sphere', 'position',[200,50]*nm, 'radius',[10,15]*nm, 'RI',3, 'loss',0*0.01)};
% Objects = {};


%% initialization
domainSize = [100,105]; % x,y size of domain
% FDTDregion = [0,400; 0,450]*nm;
% domainSize = [100,100]; % x,y size of domain
FDTDregion = [0,400; -200,200]*nm;
calFDTD = FDTD(domainSize,FDTDregion,isTE);	% creat 2D FDTD grid (isTM for 2D)
calFDTD.initGrid(Objects,[]);               % initialize the 2D FDTD grid
calFDTD.initABC('2abc');                    % initialize ABC
TFSFregion = [10,10;domainSize-[10,10]];
calFDTD.initTFSF(TFSFregion,source);        % initialize TFSF


% graphic output initialization
if isTE
    fstr = 'z';
else
    fstr = 'y';
end
[x,y] = calFDTD.getGridXYZ(['E',fstr]);
dt = calFDTD.get_dt()/fs; % in fs
figure
hs = surf(y/nm,x/nm,calFDTD.snapshotE(fstr),'linestyle','none'); % init a snapshot
colorbar
ht = title(sprintf('{E_%s}, Time = %g fs',fstr,0));
view(2); axis tight equal
% hold on
xlabel('y, nm')
ylabel('x, nm')

%% do time stepping
for n_time = 1:MaxTime
    tic
    calFDTD.updateH();          % update magnetic field
    calFDTD.updateTFSF(n_time);	% apply TFSF boundary
    calFDTD.updateE();          % update electric field
    calFDTD.applyABC();         % apply ABC
%     calFDTD.add_Ez_inc(round([SizeX/2, SizeY/2]), n_time, sourceG); % add a source
    toc

    if mod(n_time,10)==0 % take a snapshot
        set(hs,'Zdata',calFDTD.snapshotE(fstr));
        set(ht,'string',sprintf('{E_%s}, Time = %g fs',fstr,n_time*dt));
        pause(.2)
    end
    
end % of time-stepping
 
end
