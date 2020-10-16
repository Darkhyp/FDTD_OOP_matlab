function fdtd(Dimensionality)
% %#codegen
warning('off', 'MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');

if nargin==0
%     Dimensionality = 1;
    Dimensionality = 2;
%     Dimensionality = 3;
end

global h eV eps0 mu0 c imp0 nm fs StopSimulationsState PauseSimulationsState isVisualization
%% constants
eps0 = 8.8541878176209821e-12;
mu0  = 1.2566370614359173e-06;
imp0 = sqrt(mu0/eps0);
% c    = 299792458; % 1/sqrt(eps0*mu0)
c    = 1/sqrt(eps0*mu0);
nm   = 1e-9;
fs   = 1e-15;
eV   = 1.60217733e-19;
h    = 6.62606965729e-34;
% default variables
StopSimulationsState	= false;
PauseSimulationsState	= false;
PML             = [];
TFSFregion      = [];
Objects         = {};
isTE            = false;
isVisualization	= false;
isWaterFall     = true; % for 1D case

sourceR.ppw         = 20;
sourceR.omega       = 300e12; % 300/1micron = 300THz
sourceR.timeDelay	= sourceR.ppw;
sourceR.Type        = 'Ricker';
sourceR.isAdd       = true;
% sourceR.isAdd     = false;

sourceG.ppw         = 20;
sourceG.timeDelay	= 4*sourceG.ppw;
sourceG.Type        = 'gauss';
sourceG.isAdd       = true;
% sourceG.isAdd       = false;

sourceGm.ppw         = 20;
sourceGm.timeDelay	 = 4*sourceGm.ppw;
sourceGm.Type        = 'gauss_m';
sourceGm.omega       = 2*pi*300e-6/1e-6*1e12; % /1micron
sourceGm.isAdd       = true;
% sourceGm.isAdd       = false;

%%
% isWaterFall = false;
isVisualization	= true;

MaxTime = 1000; % duration of simulation
% MaxTime = 700; % duration of simulation
% MaxTime = 500; % duration of simulation

isTE = true;

% source = sourceR;
% source = sourceG;
source = sourceGm;

% PML = struct('type','split PML', 'layers',10, 'm',3.5, 'coef',0.8);
% PML = struct('type','UPML', 'layers',10, 'm',4, 'R_err',1e-16, 'ka_max',1, 'gamma',1);
PML = struct('type','UPML', 'layers',10, 'm',4, 'R_err',1e-10, 'ka_max',1.5, 'gamma',1);

domainSize = [110; 115; 120]; % x,y,z sizes of domain
% domainSize = [100; 100; 100]; % x,y,z sizes of domain
FDTDregion = [0,400; -100,300; -200,200]*nm;

% TFSFregion = [[15;15;15], domainSize-[15;15;15]];
TFSFregion = [[10;10;10], domainSize-[10;10;10]];

% objects
Objects = {struct('dispersive','Ag Palik', 'type','sphere', 'position',[200,100,0]*nm, 'radius',[24,30,25]*nm, 'RI',3, 'loss',0*0.01)};

% out objects
% OUTobjects{1} = [5,domainSize(1)-5; 5,domainSize(2)-5; 5,domainSize(3)-5];

% OUTobjects{1} = [domainSize(1)-5,domainSize(1)-5; 5,domainSize(2)-5; 5,domainSize(3)-5];

% OUTobjects{1} = [15,15; 10,domainSize(2)-10; 10,domainSize(3)-10];
% OUTobjects{2} = [domainSize(1)-15,domainSize(1)-15; 10,domainSize(2)-10; 10,domainSize(3)-10];

OUTobjects{1} = [25,25; 1,domainSize(2); 1,domainSize(3)];
OUTobjects{2} = [domainSize(1)-25,domainSize(1)-25; 1,domainSize(2); 1,domainSize(3)];


%% initialization
switch Dimensionality
    case 1
        calFDTD = FDTD1D(domainSize,FDTDregion,PML);        % creat 1D FDTD grid
    case 2
        calFDTD = FDTD2D(domainSize,FDTDregion,PML,isTE);	% creat 2D FDTD grid (isTE for 2D)
    case 3
        calFDTD = FDTD3D(domainSize,FDTDregion,PML);        % creat 3D FDTD grid
end

calFDTD.initGrid(Objects,[]);               % initialize the FDTD grid

if isempty(PML)
    calFDTD.initABC('2abc');                % initialize ABC
end

if ~isempty(TFSFregion)
    calFDTD.initTFSF(TFSFregion,source);	% initialize TFSF
end

calFDTD.initOUT(OUTobjects,MaxTime);	% initialize OUT data

% graphic output initialization
dt = calFDTD.get_dt()/fs; % in fs
if isVisualization
    hf = figure;

    colormap('jet')
    switch Dimensionality
        case 1
            x  = calFDTD.getGridXYZ('Ez')/nm; % in nm
            if isWaterFall
                waterfall(x,0,calFDTD.snapshotE('z').');
            else
                plot3(x,0*ones(1,domainSize),calFDTD.snapshotE().','color','blue');
            end
            ha = gca; hold on
            view(-20,62)
            xlabel('x, nm')
            ylabel('Time, fs')
            %zlabel('Electric field')
        case 2
            if isTE
                fstr = 'z';
            else
                fstr = 'y';
            end
            [x,y] = calFDTD.getGridXYZ(['E',fstr]);
                hs = surf(y/nm,x/nm,calFDTD.snapshotE(fstr),'linestyle','none'); % init a snapshot
                colorbar
                ht = title(sprintf('{E_%s}, Time = %g fs',fstr,0));
                view(2); axis tight equal
                % hold on
                xlabel('y, nm')
                ylabel('x, nm')
            hold on
            for n_out=1:length(calFDTD.OUT.data)
                plot(calFDTD.OUT.data{n_out}.Y/nm,calFDTD.OUT.data{n_out}.X/nm,'.-')
            end
            uicontrol('Position',[10 10 100 40],'String','Break',...
                'Callback','StopSimulations');
            uicontrol('Position',[110 10 100 40],'String','Pause',...
                'Callback','PauseSimulations');
        case 3
            if isempty(PML)
                slice_pos = round(domainSize/2);
            else
                slice_pos = round(domainSize/2+PML.layers);
            end
            Ez = calFDTD.snapshotE('z');
            ha = zeros(3);
            [x,y,z] = calFDTD.getGridXYZ('Ez');
            x = x/nm;
            y = y/nm;
            z = z/nm;
            subplot(2,2,1)
                hslice = slice(y,x,z,Ez,y(slice_pos(2)),x(slice_pos(1)),z(slice_pos(3)));
                set(hslice,'linestyle','none');
                ht = title(sprintf('Time = %g fs',0));
                colorbar
                axis tight equal
                xlabel('{\it x}, nm')
                ylabel('{\it y}, nm')
                zlabel('{\it z}, nm')
            subplot(2,2,2)
                ha(2) = surf(z,y,squeeze(Ez(slice_pos(1),:,:)),'linestyle','none'); % init a snapshot
                title('{\it E_z} y0z');
                colorbar
                view(2); axis tight equal
                xlabel('{\it z}, nm')
                ylabel('{\it y}, nm')
            subplot(2,2,3)
                ha(1) = surf(z,x,squeeze(Ez(:,slice_pos(2),:)),'linestyle','none'); % init a snapshot
                title('{\it E_z} x0z');
                colorbar
                view(2); axis tight equal
                xlabel('{\it z}, nm')
                ylabel('{\it x}, nm')
            subplot(2,2,4)
                ha(3) = surf(y,x,Ez(:,:,slice_pos(3)),'linestyle','none'); % init a snapshot
                title('{\it E_z} x0y');
                colorbar
                view(2); axis tight equal
                xlabel('{\it y}, nm')
                ylabel('{\it x}, nm')
            uicontrol('Position',[10 10 100 40],'String','Break',...
                'Callback','StopSimulations');
            uicontrol('Position',[110 10 100 40],'String','Pause',...
                'Callback','PauseSimulations');
    end
    jFrame = get(hf,'JavaFrame');
%     jFrame.getFigurePanelContainer.getTopLevelAncestor.setMaximized(true);    % maximize window [true/false]
else
    h_progress = waitbar(0, 'Calculations in progress...');
end

%% do time stepping
fullFDTDcalTime = tic;
for n_time = 1:MaxTime
    if(StopSimulationsState),  break;    end
    if(PauseSimulationsState), continue; end
    
    tic
    calFDTD.updateH();          % update magnetic field

    if isempty(TFSFregion)
        if Dimensionality==1
            calFDTD.add_Hy_inc(round(domainSize/2)+PML.layers-1, n_time-1, source); % add a source
        end
    else
        calFDTD.updateTFSF(n_time);	% apply TFSF boundary
    end
    
    calFDTD.updateE();          % update electric field
    
    if isempty(TFSFregion)
        calFDTD.add_Ez_inc(round(domainSize/2)+PML.layers, n_time, source); % add a source
    end
    
    if isempty(PML)
        calFDTD.applyABC();         % apply ABC
    end

    calFDTD.updateOUT(n_time);      % update OUT data
    toc

    if mod(n_time,10)==0 % take a snapshot
        if isVisualization
            switch Dimensionality
                case 1
                    if isWaterFall
                        waterfall(ha,x,n_time*dt,calFDTD.snapshotE().');
                    else
                        plot3(ha,x,n_time*dt*ones(1,domainSize),calFDTD.snapshotE('z').','color','blue');
                    end
                    for n_out=1:length(calFDTD.OUT.data)
                        scatter3(calFDTD.OUT.data{n_out}.X/nm,n_time*dt,calFDTD.OUT.data{n_out}.Sx(n_time),'o')
                    end
                case 2
                    set(hs,'Zdata',calFDTD.snapshotE(fstr));
                    set(ht,'string',sprintf('Time(%i of %i) = %g fs',n_time,MaxTime,n_time*dt));
                    pause(.2)
                case 3
                    Ez = calFDTD.snapshotE('z');
                    set(ha(2),'Zdata',squeeze(Ez(slice_pos(1),:,:)));
                    set(ha(1),'Zdata',squeeze(Ez(:,slice_pos(2),:))); 
                    set(ha(3),'Zdata',Ez(:,:,slice_pos(3)));
                    for n=1:3
                        set(hslice(n),'Cdata',get(ha(n),'Zdata')); 
                    end
                    set(ht,'string',sprintf('Time(%i of %i) = %g fs',n_time,MaxTime,n_time*dt));
                    pause(0.001)
            end
            %{
            figure, 
                hold on
                for n_out=1:length(calFDTD.OUT.data)
                    plot((0:n_time-1)*dt,calFDTD.OUT.data{n_out}.Sx_aver,'.-')
                end
            for n_out=1:length(calFDTD.OUT.data)
                figure, surf(squeeze(calFDTD.OUT.data{n_out}.Y),squeeze(calFDTD.OUT.data{n_out}.Z),calFDTD.OUT.data{n_out}.Sx{n_time},'linestyle','none')
            end
            %}
        else
            % Progress bar updates at some time
            waitbar(n_time/MaxTime, h_progress);
        end
    end    
end % of time-stepping
fprintf('\n Full calulation time is %gs.',toc(fullFDTDcalTime))
 
%% FFT
t = (0:MaxTime-1)*dt;
NFFT = 2^10;
NFFT = 2^nextpow2(MaxTime);
% Frequency = (-NFFT/2:NFFT/2-1)/NFFT/dt;	 	 
% Frequency = (0:NFFT-1)/NFFT/dt;	 	 
Frequency = linspace(0,1,NFFT/2)/2/dt;	 	 
% Courant = calFDTD.get_Courant(); Courant = Courant(diff(calFDTD.OUT.Objects{1},[],2)==0);

% Fourier transformations
for n_out=length(calFDTD.OUT.data):-1:1
%     Sout(:,n_out) = fftshift(fft(calFDTD.OUT.data{n_out}.Ez,NFFT));
%     Sout(:,n_out) = fftshift(fft(calFDTD.OUT.data{n_out}.Sx_aver,NFFT));
    Sout(:,n_out) = fft(calFDTD.OUT.data{n_out}.Sx_aver,NFFT);
end

figure
if Dimensionality<3
    subplot(4,1,1)
        hold on
        for n_out=1:length(calFDTD.OUT.data)
            plot(t,calFDTD.OUT.data{n_out}.Ez, '.-');
    %         plot(t,calFDTD.OUT.data{n_out}.Sx_aver, '.-');
        end
        legend('Incident field', 'Transmitted field');
        xlabel('Time, fs');
        ylabel('{\it E_z}({\it x_i},{\it t}), V/m');
        title('Temporal Fields');
    subplot(4,1,2)
else
    subplot(3,1,1)
end
%{
figure
    n_out = 1;
%     [~,hs] = contourf(squeeze(calFDTD.OUT.data{n_out}.Y)/nm, squeeze(calFDTD.OUT.data{n_out}.Z)/nm, calFDTD.OUT.data{n_out}.Sx(:,:,1), 100, 'linestyle','none')
    hs = surf(squeeze(calFDTD.OUT.data{n_out}.Y)/nm, squeeze(calFDTD.OUT.data{n_out}.Z)/nm, calFDTD.OUT.data{n_out}.Sx(:,:,1), 'linestyle','none')
%     hs = surf(squeeze(calFDTD.OUT.data{n_out}.Y)/nm, squeeze(calFDTD.OUT.data{n_out}.Z)/nm, calFDTD.OUT.data{n_out}.Ez(:,:,1), 'linestyle','none')
    view(2);
%     axis tight equal
    xlabel('{\it y}, nm')
    ylabel('{\it z}, nm')
    ht = title(sprintf('Time(%i of %i) = %g fs',1,MaxTime,1*dt));
    colorbar
    pause(.05)

for n_time=2:MaxTime
    hs(1).ZData = calFDTD.OUT.data{n_out}.Sx(:,:,n_time);
%     hs(1).ZData = calFDTD.OUT.data{n_out}.Ez(:,:,n_time);
    ht.String = sprintf('Time(%i of %i) = %g fs',n_time,MaxTime,n_time*dt);
    pause(.05)
end
%}

    hold on
    for n_out=1:length(calFDTD.OUT.data)
        plot(t,calFDTD.OUT.data{n_out}.Sx_aver, '.-');
    end
    legend('Incident field', 'Transmitted field');
    xlabel('Time, fs');
    ylabel('{\it E_z}({\it x_i},{\it t}), V/m');
    title('Temporal Poynting');
subplot(4,1,3)
    hold on
    for n_out=1:length(calFDTD.OUT.data)
        plot(Frequency/fs/1e12, abs(Sout(1:NFFT/2,n_out)), '.-');
    end
    legend('Incident field', 'Transmitted field');
    xlabel('Frequency, THz');
    ylabel('|{\it E_z}({\it x_i},\omega)|');
    title('Spectra');
    tmp = xlim;
    xlim([0,tmp(2)])
subplot(4,1,4)
    hold on
    plot(Frequency/fs/1e12, abs( Sout(1:NFFT/2,2)./Sout(1:NFFT/2,1) ),'.-');
    xlabel('Frequency, THz');
    ylabel('|{\it E_z}({\it x}_1,\omega)/{\it E^i_z}({\it x}_2,\omega)|');
    title('Normalized spectra');
    tmp = xlim;
    xlim([0,tmp(2)])
    ylim([0,1])

end


