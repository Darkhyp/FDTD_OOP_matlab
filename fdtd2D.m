function fdtd2D()
% %#codegen

global eps0 mu0 c imp0 nm fs StopSimulationsState PauseSimulationsState
load_global_constants();
StopSimulationsState	= false;
PauseSimulationsState	= false;

PML = [];
TFSFregion = [];
isTE = true;
% isTE = false;

sourceR.ppw         = 20;
sourceR.omega       = 300e12; % 300/1micron = 300THz
sourceR.timeDelay	= sourceR.ppw;
sourceR.Type        = 'Ricker';
sourceR.isAdd       = true;
% sourceR.isAdd       = false;

sourceG.ppw         = 20;
sourceG.omega       = 300e12; % 300/1micron = 300THz
sourceG.timeDelay	= 4*sourceG.ppw;
sourceG.Type        = 'gauss';
sourceG.isAdd       = true;
% sourceG.isAdd       = false;

% source = sourceR;
source = sourceG;


% FileName = 'init.dat';
% loadDataFromFile(FileName);

MaxTime = 500; % duration of simulation

% objects
Objects = {struct('type','sphere', 'position',[200,50]*nm, 'radius',[10,15]*nm, 'RI',3, 'loss',0*0.01)};
Objects = {};


%% initialization
% domainSize = [100;105]; % x,y size of domain
% FDTDregion = [0,400; 0,450]*nm;
domainSize = [100;100]; % x,y size of domain
% FDTDregion = [0,400; -200,200]*nm;
FDTDregion = [-200,200; -200,200]*nm;
PML = struct('layers',10, 'm',3.5, 'coef',0.8);
calFDTD = FDTD(domainSize,FDTDregion,PML,isTE);	% creat 2D FDTD grid (isTM for 2D)
calFDTD.initGrid(Objects,[]);               % initialize the 2D FDTD grid
if isempty(PML)
    calFDTD.initABC('2abc');                % initialize ABC
end
% TFSFregion = [[10;10], domainSize-[10;10]];
if ~isempty(TFSFregion)
    calFDTD.initTFSF(TFSFregion,source);        % initialize TFSF
end


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
uicontrol('Position',[10 10 100 40],'String','Break',...
    'Callback','StopSimulations');
uicontrol('Position',[110 10 100 40],'String','Pause',...
    'Callback','PauseSimulations');

%% do time stepping
for n_time = 1:MaxTime
    if(StopSimulationsState), break; end
    if(PauseSimulationsState), continue; end

    tic
    calFDTD.updateH();          % update magnetic field
    if ~isempty(TFSFregion)
        calFDTD.updateTFSF(n_time);	% apply TFSF boundary
    end
    calFDTD.updateE();          % update electric field
    if isempty(TFSFregion)
        calFDTD.add_Ez_inc(round(domainSize/2)+PML.layers, n_time, source); % add a source
    end
    if isempty(PML)
        calFDTD.applyABC();         % apply ABC
    end
    toc

    if mod(n_time,10)==0 % take a snapshot
        set(hs,'Zdata',calFDTD.snapshotE(fstr));
        set(ht,'string',sprintf('Time(%i of %i) = %g fs',n_time,MaxTime,n_time*dt));
        pause(.2)
    end
    
end % of time-stepping
 
end
