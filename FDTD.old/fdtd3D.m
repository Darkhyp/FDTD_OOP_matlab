function fdtd3D()
%#codegen

global eps0 mu0 c nm fs
eps0 = 8.8541878176209821e-12;
mu0  = 1.2566370614359173e-06;
c    = 299792458; % 1/sqrt(eps0*mu0)
nm   = 1e-9;
fs   = 1e-15;

source.ppw          = 20;
source.timeDelay	= source.ppw;
% sourceType        = 'gauss';
source.Type         = 'Ricker';
source.isAdd        = false;
isTE = false;
isTE = true;
    

sourceG.ppw          = 20;
sourceG.timeDelay	= 3*source.ppw;
sourceG.Type        = 'gauss';
sourceG.isAdd        = true;


MaxTime = 1000; % duration of simulation

% objects
Objects = {struct('type','sphere', 'position',[200,0,0]*nm, 'radius',[10*nm,15*nm,20*nm], 'RI',3, 'loss',0*0.01)};
Objects = {};


%% initialization
domainSize = [103, 100, 97]; % x,y,z sizes of domain
% domainSize = [100, 100, 100]; % x,y,z sizes of domain
FDTDregion = [0,400; -100,300; -200,200]*nm;
calFDTD = FDTD(domainSize,FDTDregion);	% creat 3D FDTD grid (isTE for 2D)
calFDTD.initGrid(Objects,[]);           % initialize the 3D FDTD grid
calFDTD.initABC('2abc');                % initialize ABC
TFSFregion = [10,10,10;domainSize(1)-[10,10,10]];
calFDTD.initTFSF(TFSFregion,source);	% initialize TFSF

% graphic output initialization
slice_pos = round(domainSize/2);
Ez = calFDTD.snapshotE('z');
hf = zeros(3);
[x,y,z] = calFDTD.getGridXYZ('Ez');
x = x/nm;
y = y/nm;
z = z/nm;
dt = calFDTD.get_dt()/fs; % in fs
figure
subplot(221)
    hslice = slice(y,x,z,Ez,y(slice_pos(2)),x(slice_pos(1)),z(slice_pos(3)));
    set(hslice,'linestyle','none');
    ht = title(sprintf('Time = %g fs',0));
    colorbar
    axis tight equal
    xlabel('x, nm')
    ylabel('y, nm')
    zlabel('z, nm')
subplot(222)
    hf(2) = surf(z,y,squeeze(Ez(slice_pos(1),:,:)),'linestyle','none'); % init a snapshot
    title('E_z y0z');
    colorbar
    view(2); axis tight equal
    xlabel('z, nm')
    ylabel('y, nm')
subplot(223)
    hf(1) = surf(z,x,squeeze(Ez(:,slice_pos(2),:)),'linestyle','none'); % init a snapshot
    title('E_z x0z');
    colorbar
    view(2); axis tight equal
    xlabel('z, nm')
    ylabel('x, nm')
subplot(224)
    hf(3) = surf(y,x,Ez(:,:,slice_pos(3)),'linestyle','none'); % init a snapshot
    title('E_z x0y');
    colorbar
    view(2); axis tight equal
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
        Ez = calFDTD.snapshotE('z');
        set(hf(2),'Zdata',squeeze(Ez(slice_pos(1),:,:)));
        set(hf(1),'Zdata',squeeze(Ez(:,slice_pos(2),:))); 
        set(hf(3),'Zdata',Ez(:,:,slice_pos(3)));
        for n=1:3
            set(hslice(n),'Cdata',get(hf(n),'Zdata')); 
        end
        set(ht,'string',sprintf('Time = %g fs',n_time*dt));
        pause(0.2)
    end
    
end % of time-stepping
 
end
