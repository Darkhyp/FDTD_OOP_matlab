function fdtd3D()
% %#codegen

global eps0 mu0 c imp0 nm fs StopSimulationsState PauseSimulationsState
load_global_constants();
StopSimulationsState	= false;
PauseSimulationsState	= false;

sourceR.ppw         = 20;
sourceR.omega       = 300e12; % 300/1micron = 300THz
sourceR.timeDelay	= sourceR.ppw;
sourceR.Type        = 'Ricker';
sourceR.isAdd       = true;
% sourceR.isAdd     = false;

sourceG.ppw         = 20;
sourceG.omega       = 300e12; % 300/1micron = 300THz
sourceG.timeDelay	= 4*sourceG.ppw;
sourceG.Type        = 'gauss';
sourceG.isAdd       = true;
% sourceG.isAdd       = false;

% source = sourceR;
source = sourceG;

PML = [];
TFSFregion = [];
isTE = false;
isTE = true;

% MaxTime = 1000; % duration of simulation
MaxTime = 500; % duration of simulation

% objects
% Objects = {struct('type','sphere', 'position',[200,0,0]*nm, 'radius',[10,15,20]*nm, 'RI',3, 'loss',0*0.01)};
Objects = {};


%% initialization
% domainSize = [103; 100; 97]; % x,y,z sizes of domain
domainSize = [100; 100; 100]; % x,y,z sizes of domain
% FDTDregion = [0,400; -100,300; -200,200]*nm;
FDTDregion = [-200,200; -200,200; -200,200]*nm;
PML = struct('layers',10, 'm',3.5, 'coef',0.8);
calFDTD = FDTD(domainSize,FDTDregion,PML);	% creat 3D FDTD grid (isTE for 2D)
calFDTD.initGrid(Objects,[]);           % initialize the 3D FDTD grid
if isempty(PML)
    calFDTD.initABC('2abc');                % initialize ABC
end
TFSFregion = [[10;10;10], domainSize-[10;10;10]];
if ~isempty(TFSFregion)
    calFDTD.initTFSF(TFSFregion,source);	% initialize TFSF
end

% SCATregion = [5,5,5;domainSize-[5,5,5]];
% calFDTD.initSCAT(SCATregion);           % initialize SCAT data

% OUTobjects{1} = [5,domainSize(1)-5; 5,domainSize(2)-5; 5,domainSize(3)-5];

% OUTobjects{1} = [domainSize(1)-5,domainSize(1)-5; 5,domainSize(2)-5; 5,domainSize(3)-5];

% OUTobjects{1} = [15,15; 10,domainSize(2)-10; 10,domainSize(3)-10];
% OUTobjects{2} = [domainSize(1)-15,domainSize(1)-15; 10,domainSize(2)-10; 10,domainSize(3)-10];

OUTobjects{1} = [15,15; 1,domainSize(2); 1,domainSize(3)];
OUTobjects{2} = [domainSize(1)-15,domainSize(1)-15; 1,domainSize(2); 1,domainSize(3)];

calFDTD.initOUT(OUTobjects,MaxTime);	% initialize OUT data

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
    hf(2) = surf(z,y,squeeze(Ez(slice_pos(1),:,:)),'linestyle','none'); % init a snapshot
    title('{\it E_z} y0z');
    colorbar
    view(2); axis tight equal
    xlabel('{\it z}, nm')
    ylabel('{\it y}, nm')
subplot(2,2,3)
    hf(1) = surf(z,x,squeeze(Ez(:,slice_pos(2),:)),'linestyle','none'); % init a snapshot
    title('{\it E_z} x0z');
    colorbar
    view(2); axis tight equal
    xlabel('{\it z}, nm')
    ylabel('{\it x}, nm')
subplot(2,2,4)
    hf(3) = surf(y,x,Ez(:,:,slice_pos(3)),'linestyle','none'); % init a snapshot
    title('{\it E_z} x0y');
    colorbar
    view(2); axis tight equal
    xlabel('{\it y}, nm')
    ylabel('{\it x}, nm')
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
%     calFDTD.updateSCAT();       % update SCAT data

    calFDTD.updateOUT(n_time);	% update OUT data
    toc

    if mod(n_time,10)==0 % take a snapshot
        Ez = calFDTD.snapshotE('z');
        set(hf(2),'Zdata',squeeze(Ez(slice_pos(1),:,:)));
        set(hf(1),'Zdata',squeeze(Ez(:,slice_pos(2),:))); 
        set(hf(3),'Zdata',Ez(:,:,slice_pos(3)));
        for n=1:3
            set(hslice(n),'Cdata',get(hf(n),'Zdata')); 
        end
        set(ht,'string',sprintf('Time(%i of %i) = %g fs',n_time,MaxTime,n_time*dt));
        pause(0.001)
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
    end
    
end % of time-stepping
 
%% FFT
t = (0:MaxTime-1)*dt;
NFFT = 2^10;	 	 
Frequency = (-NFFT/2:NFFT/2-1)/NFFT/dt;	 	 
% Frequency = (0:NFFT-1)/NFFT/dt;	 	 
courantNumber = calFDTD.get_CourantNumber();
courantNumber = courantNumber(diff(calFDTD.OUT.Objects{1},[],2)==0);

% Fourier transformations
for n_out=length(calFDTD.OUT.data):-1:1
    Sout(:,n_out) = fftshift(fft(calFDTD.OUT.data{n_out}.Sx_aver,NFFT));
%     Sout(:,n_out) = fft(calFDTD.OUT.data{n_out}.Sx_aver,NFFT);
end

figure
subplot(3, 1, 1)
    hold on
    for n_out=1:length(calFDTD.OUT.data)
        plot(t,calFDTD.OUT.data{n_out}.Sx_aver, '.-');
    end
    legend('Incident field', 'Transmitted field');
    xlabel('Time, fs');
    ylabel('{\it E_z}({\it x_i},{\it t}), V/m');
    title('Temporal Fields');

subplot(3,1,2)
    hold on
    for n_out=1:length(calFDTD.OUT.data)
        plot(Frequency*2*pi/1000, abs(Sout(:,n_out)), '.-');
    end
    legend('Incident field', 'Transmitted field');
    xlabel('Frequency, THz');
    ylabel('|{\it E_z}({\it x_i},\omega)|');
    title('Spectra');
    tmp = xlim;
    xlim([0,tmp(2)])

subplot(3, 1, 3)
hold on
    plot(Frequency*2*pi/1000, abs( Sout(:,2)./Sout(:,1) ));
    xlabel('Frequency, THz');
    ylabel('|{\it E_z}({\it x}_1,\omega)/{\it E^i_z}({\it x}_2,\omega)|');
    title('Normalized spectra');
    tmp = xlim;
    xlim([0,tmp(2)])
    ylim([0,1])

end


