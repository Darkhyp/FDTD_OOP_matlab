%% Apply TFSF
function obj = updateTFSF(obj,time)
%#codegen

if isempty(obj.TFSF) % check if initTFSF() has been called
	fprintf('\nUse initTFSF() method before updateTFSF()!...')
%                 exit(-1);
end

if obj.isTE
    % correct H adjacent to TFSF boundary
    %**** constant y faces (scattered-field nodes) ****
    m = obj.TFSF.region(1,1):obj.TFSF.region(1,2);
    n_1 = obj.TFSF.region(2,1)-1;
    n_2 = obj.TFSF.region(2,2);
    % correct Hx at firstY-1/2 by adding Ez_inc
    obj.Hx(m,n_1) = obj.Hx(m,n_1) + obj.coefHx_z(m,n_1).*obj.TFSF.fields1D.Ez(m);
    % correct Hx at lastY+1/2 by subtracting Ez_inc
    obj.Hx(m,n_2) = obj.Hx(m,n_2) - obj.coefHx_z(m,n_2).*obj.TFSF.fields1D.Ez(m);

    %****** constant x faces (scattered-field nodes) ******
    m_1 = obj.TFSF.region(1,1)-1;
    m_2 = obj.TFSF.region(1,2);
    n = obj.TFSF.region(2,1):obj.TFSF.region(2,2);
    % correct Hy at firstX-1/2 by subtracting Ez_inc
    obj.Hy(m_1,n) = obj.Hy(m_1,n) - obj.coefHy_z(m_1,n).*obj.TFSF.fields1D.Ez(m_1+1);
    % correct Hy at lastX + 1/2 by adding Ez_inc
    obj.Hy(m_2,n) = obj.Hy(m_2,n) + obj.coefHy_z(m_2,n).*obj.TFSF.fields1D.Ez(m_2);

    % update 1D
    obj.TFSF.fields1D.updateH();         % update 1D magnetic field
    obj.TFSF.fields1D.updateTFSF(time);	% in 1D, apply TFSF
    obj.TFSF.fields1D.updateE();         % update 1D electric field
%                         obj.TFSF.fields1D.Ez(3) = SourceFunction(time, obj.TFSF.source,obj.CourantNumber); % set source node
    if isempty(obj.PML)
        obj.TFSF.fields1D.applyABC();        % apply ABC
    end

    % correct E adjacent to TFSF boundary

    %**** constant x faces (total-field nodes) ****
    m_1 = obj.TFSF.region(1,1);
    m_2 = obj.TFSF.region(1,2);
    % correct Ez at firstX face by subtracting Hy_inc
    obj.Ez(m_1,n) = obj.Ez(m_1,n) - obj.coefEz_y(m_1-1,n-1).*obj.TFSF.fields1D.Hy(m_1-1);
    % correct Ez at lastX face by adding Hy_inc
    obj.Ez(m_2,n) = obj.Ez(m_2,n) + obj.coefEz_y(m_2-1,n-1).*obj.TFSF.fields1D.Hy(m_2);

    % no need to correct Ez along top and bottom since
    % incident Hx is zero
else
    % correct Hz along x-axis (start)
    m_1 = obj.TFSF.region(1,1)-1;
    m_2 = obj.TFSF.region(1,2);
    n = obj.TFSF.region(2,1):obj.TFSF.region(2,2)-1;
    obj.Hz(m_1,n) = obj.Hz(m_1,n) + obj.coefHz_y(m_1,n).*(-obj.TFSF.fields1D.Ez(m_1+1)); % -Ez <- Ey
    % correct Hz along x-axis (end)
    obj.Hz(m_2,n) = obj.Hz(m_2,n) - obj.coefHz_y(m_2,n).*(-obj.TFSF.fields1D.Ez(m_2)); % -Ez <- Ey

    % update 1D
    obj.TFSF.fields1D.updateH();         % update 1D magnetic field
    obj.TFSF.fields1D.updateTFSF(time);	% in 1D, apply TFSF
    obj.TFSF.fields1D.updateE();         % update 1D electric field
%                         obj.TFSF.fields1D.Ey(1) = SourceFunction(time,obj.TFSF.source,obj.CourantNumber); %<-Ey in 1D, set source node
    if isempty(obj.PML)
        obj.TFSF.fields1D.applyABC();        % apply ABC
    end

    % correct Ey field along x-axis (start)
    m_1 = obj.TFSF.region(1,1);
    m_2 = obj.TFSF.region(1,2);
    obj.Ey(m_1,n) = obj.Ey(m_1,n) + obj.coefEy_z(m_1-1,n).*obj.TFSF.fields1D.Hy(m_1-1); % Hy <- Hz
    % correct Ey field along y-axis (end)
    obj.Ey(m_2,n) = obj.Ey(m_2,n) - obj.coefEy_z(m_2-1,n).*obj.TFSF.fields1D.Hy(m_2); % Hy <- Hz

    % correct Ex along y-axis (start)
    m = obj.TFSF.region(1,1):obj.TFSF.region(1,2)-1;
    n_1 = obj.TFSF.region(2,1);
    n_2 = obj.TFSF.region(2,2);
    obj.Ex(m,n_1) = obj.Ex(m,n_1) - obj.coefEx_z(m,n_1-1).*obj.TFSF.fields1D.Hy(m); % Hy <- Hz
    % correct Ex along y-axis (end)
    obj.Ex(m,n_2) = obj.Ex(m,n_2) + obj.coefEx_z(m,n_2-1).*obj.TFSF.fields1D.Hy(m); % Hy <- Hz

    % no need to correct Ex along top and bottom since
    % incident Ex is zero
end

end
