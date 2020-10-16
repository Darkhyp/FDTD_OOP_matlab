%% Apply TFSF 
function obj = updateTFSF(obj,time)
%#codegen
if isempty(obj.TFSF) % check if initTFSF() has been called
	fprintf('\nUse initTFSF() method before updateTFSF()!...')
%                 exit(-1);
end

% correct H adjacent to TFSF boundary

%**** constant y faces (scattered-field nodes) ****
m = obj.TFSF.region(1,1):obj.TFSF.region(1,2);
n_1 = obj.TFSF.region(2,1)-1;
n_2 = obj.TFSF.region(2,2);
for p = obj.TFSF.region(3,1):obj.TFSF.region(3,2)-1;
	% correct Hx at firstY-1/2 by adding Ez_inc
	obj.Hx(m,n_1,p) = obj.Hx(m,n_1,p) + obj.coefHx_z(m,n_1,p).*obj.TFSF.fields1D.Ez(m);
	% correct Hx at lastY+1/2 by subtracting Ez_inc
	obj.Hx(m,n_2,p) = obj.Hx(m,n_2,p) - obj.coefHx_z(m,n_2,p).*obj.TFSF.fields1D.Ez(m);
end

%****** constant x faces (scattered-field nodes) ******
m_1 = obj.TFSF.region(1,1)-1;
m_2 = obj.TFSF.region(1,2);
n = obj.TFSF.region(2,1):obj.TFSF.region(2,2);
p = obj.TFSF.region(3,1):obj.TFSF.region(3,2)-1;
% correct Hy at firstX-1/2 by subtracting Ez_inc
obj.Hy(m_1,n,p) = obj.Hy(m_1,n,p) - obj.coefHy_z(m_1,n,p).*obj.TFSF.fields1D.Ez(m_1+1);
% correct Hy at lastX + 1/2 by adding Ez_inc
obj.Hy(m_2,n,p) = obj.Hy(m_2,n,p) + obj.coefHy_z(m_2,n,p).*obj.TFSF.fields1D.Ez(m_2);

%**** constant z faces (scattered-field nodes) ****
% nothing to correct on this face

%**** update the fields in the auxiliary 1D grid ****
obj.TFSF.fields1D.updateH();         % update 1D magnetic field
obj.TFSF.fields1D.updateTFSF(time);	% in 1D, apply TFSF
obj.TFSF.fields1D.updateE();         % update 1D electric field
%                     obj.TFSF.fields1D.Ez(1) = SourceFunction(time,obj.TFSF.source,obj.CourantNumber); %<-Ey in 1D, set source node
if isempty(obj.PML)
	obj.TFSF.fields1D.applyABC();        % apply ABC
end

% correct E adjacent to TFSF boundary

%**** constant x faces (total-field nodes) ****
m_1 = obj.TFSF.region(1,1);
m_2 = obj.TFSF.region(1,2);
n = obj.TFSF.region(2,1):obj.TFSF.region(2,2);
p = obj.TFSF.region(3,1):obj.TFSF.region(3,2)-1;
% correct Ez at firstX face by subtracting Hy_inc
obj.Ez(m_1,n,p) = obj.Ez(m_1,n,p) - obj.coefEz_y(m_1-1,n-1,p).*obj.TFSF.fields1D.Hy(m_1-1);
% correct Ez at lastX face by adding Hy_inc
obj.Ez(m_2,n,p) = obj.Ez(m_2,n,p) + obj.coefEz_y(m_2-1,n-1,p).*obj.TFSF.fields1D.Hy(m_2);

%**** constant y faces (total-field nodes) ****
% nothing to correct on this face

%**** constant z faces (total-field nodes) ****
m = obj.TFSF.region(1,1):obj.TFSF.region(1,2)-1;
p_1 = obj.TFSF.region(3,1);
p_2 = obj.TFSF.region(3,2);
for n = obj.TFSF.region(2,1):obj.TFSF.region(2,2);
	% correct Ex at firstZ face by adding Hy_inc
	obj.Ex(m,n,p_1) = obj.Ex(m,n,p_1) + obj.coefEx_y(m,n-1,p_1-1).*obj.TFSF.fields1D.Hy(m);
	% correct Ex at lastZ face by subtracting Hy_inc
	obj.Ex(m,n,p_2) = obj.Ex(m,n,p_2) - obj.coefEx_y(m,n-1,p_2-1).*obj.TFSF.fields1D.Hy(m);
end

end
