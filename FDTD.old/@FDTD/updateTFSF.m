%% Apply 
function obj = updateTFSF(obj,time)
    if isempty(obj.TFSF) % check if initTFSF() has been called
        fprintf('\nUse initTFSF() method before updateTFSF()!...')
%                 exit(-1);
    end

    switch obj.Dimensionality
        case 1
            % correct Hy adjacent to TFSF boundary
            obj.Hy(obj.TFSF.region) = obj.Hy(obj.TFSF.region) ...
                - obj.coefHy_z(obj.TFSF.region)*SourceFunction(time, obj.TFSF.source,obj.CourantNumber(1));
            % correct Ez adjacent to TFSF boundary (time+0.5-(-0.5))
            obj.Ez(obj.TFSF.region+1) = obj.Ez(obj.TFSF.region+1) ...
                + SourceFunction(time+1,obj.TFSF.source,obj.CourantNumber(1));
        case 2
            if obj.isTE
                % correct H adjacent to TFSF boundary

                %****** constant x faces (scattered-field nodes) ******
                m_1 = obj.TFSF.region(1,1)-1;
                m_2 = obj.TFSF.region(1,2);
                n = obj.TFSF.region(2,1):obj.TFSF.region(2,2);
                % correct Hy at firstX-1/2 by subtracting Ez_inc
                obj.Hy(m_1,n) = obj.Hy(m_1,n) - obj.coefHy_z(m_1,n).*obj.TFSF.fields1D.Ez(m_1+1);
                % correct Hy at lastX + 1/2 by adding Ez_inc
                obj.Hy(m_2,n) = obj.Hy(m_2,n) + obj.coefHy_z(m_2,n).*obj.TFSF.fields1D.Ez(m_2);

                %**** constant y faces (scattered-field nodes) ****
                m = obj.TFSF.region(1,1):obj.TFSF.region(1,2);
                n_1 = obj.TFSF.region(2,1)-1;
                n_2 = obj.TFSF.region(2,2);
                % correct Hx at firstY-1/2 by adding Ez_inc
                obj.Hx(m,n_1) = obj.Hx(m,n_1) + obj.coefHx_z(m,n_1).*obj.TFSF.fields1D.Ez(m);
                % correct Hx at lastY+1/2 by subtracting Ez_inc
                obj.Hx(m,n_2) = obj.Hx(m,n_2) - obj.coefHx_z(m,n_2).*obj.TFSF.fields1D.Ez(m);

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
                n = obj.TFSF.region(2,1):obj.TFSF.region(2,2);
                % correct Ez at firstX face by subtracting Hy_inc
                obj.Ez(m_1,n) = obj.Ez(m_1,n) - obj.coefEz_x(m_1-1,n-1).*obj.TFSF.fields1D.Hy(m_1-1);
                % correct Ez at lastX face by adding Hy_inc
                obj.Ez(m_2,n) = obj.Ez(m_2,n) + obj.coefEz_x(m_2-1,n-1).*obj.TFSF.fields1D.Hy(m_2);

                % no need to correct Ez along top and bottom since
                % incident Hx is zero

            else
                % correct Hz along x-axis (start)
                m = obj.TFSF.region(1,1)-1;
                n = obj.TFSF.region(2,1):obj.TFSF.region(2,2)-1;
                obj.Hz(m,n) = obj.Hz(m,n) + obj.coefHz_y(m,n).*(-obj.TFSF.fields1D.Ez(m+1)); % Ez <- (-Ey)
                % correct Hz along x-axis (end)
                m = obj.TFSF.region(1,2);
                obj.Hz(m,n) = obj.Hz(m,n) - obj.coefHz_y(m,n).*(-obj.TFSF.fields1D.Ez(m)); % Ez <- (-Ey)

                % update 1D
                obj.TFSF.fields1D.updateH();         % update 1D magnetic field
                obj.TFSF.fields1D.updateTFSF(time);	% in 1D, apply TFSF
                obj.TFSF.fields1D.updateE();         % update 1D electric field
%                         obj.TFSF.fields1D.Ey(1) = SourceFunction(time,obj.TFSF.source,obj.CourantNumber); %<-Ey in 1D, set source node
                if isempty(obj.PML)
                    obj.TFSF.fields1D.applyABC();        % apply ABC
                end

                % correct Ex along y-axis (start)
                m = obj.TFSF.region(1,1):obj.TFSF.region(1,2)-1;
                n = obj.TFSF.region(2,1);
                obj.Ex(m,n) = obj.Ex(m,n) - obj.coefEx_z(m,n-1).*obj.TFSF.fields1D.Hy(m); % Hy <- Hz
                % correct Ex along y-axis (end)
                n = obj.TFSF.region(2,2);
                obj.Ex(m,n) = obj.Ex(m,n) + obj.coefEx_z(m,n-1).*obj.TFSF.fields1D.Hy(m); % Hy <- Hz

                % correct Ey field along x-axis (start)
                m = obj.TFSF.region(1,1);
                n = obj.TFSF.region(2,1):obj.TFSF.region(2,2)-1;
                obj.Ey(m,n) = obj.Ey(m,n) + obj.coefEy_z(m-1,n).*obj.TFSF.fields1D.Hy(m-1); % Hy <- Hz
                % correct Ey field along y-axis (end)
                m = obj.TFSF.region(1,2);
                obj.Ey(m,n) = obj.Ey(m,n) - obj.coefEy_z(m-1,n).*obj.TFSF.fields1D.Hy(m); % Hy <- Hz

                % no need to correct Ex along top and bottom since
                % incident Ex is zero

            end
        case 3
            % correct H adjacent to TFSF boundary

            %****** constant x faces (scattered-field nodes) ******
            m_1 = obj.TFSF.region(1,1)-1;
            m_2 = obj.TFSF.region(1,2);
            n = obj.TFSF.region(2,1):obj.TFSF.region(2,2);
            p = obj.TFSF.region(3,1):obj.TFSF.region(3,2)-1;
            % correct Hy at firstX-1/2 by subtracting Ez_inc
            obj.Hy(m_1,n,p) = obj.Hy(m_1,n,p) - obj.coefHy_x(m_1,n,p).*obj.TFSF.fields1D.Ez(m_1+1);
            % correct Hy at lastX + 1/2 by adding Ez_inc
            obj.Hy(m_2,n,p) = obj.Hy(m_2,n,p) + obj.coefHy_x(m_2,n,p).*obj.TFSF.fields1D.Ez(m_2);

            %**** constant y faces (scattered-field nodes) ****
            m = obj.TFSF.region(1,1):obj.TFSF.region(1,2);
            n_1 = obj.TFSF.region(2,1)-1;
            n_2 = obj.TFSF.region(2,2);
            for p = obj.TFSF.region(3,1):obj.TFSF.region(3,2)-1;
                % correct Hx at firstY-1/2 by adding Ez_inc
                obj.Hx(m,n_1,p) = obj.Hx(m,n_1,p) + obj.coefHx_y(m,n_1,p).*obj.TFSF.fields1D.Ez(m);
                % correct Hx at lastY+1/2 by subtracting Ez_inc
                obj.Hx(m,n_2,p) = obj.Hx(m,n_2,p) - obj.coefHx_y(m,n_2,p).*obj.TFSF.fields1D.Ez(m);
            end

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
            obj.Ez(m_1,n,p) = obj.Ez(m_1,n,p) - obj.coefEz_x(m_1-1,n-1,p).*obj.TFSF.fields1D.Hy(m_1-1);
            % correct Ez at lastX face by adding Hy_inc
            obj.Ez(m_2,n,p) = obj.Ez(m_2,n,p) + obj.coefEz_x(m_2-1,n-1,p).*obj.TFSF.fields1D.Hy(m_2);

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
end
