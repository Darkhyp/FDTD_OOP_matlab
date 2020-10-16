%% calculate source function at given time and location
function obj = add_Ez_inc(obj,Position,time,sourceStr)
%             if (sourceStr.ppw <= 0)
%                 fprintf(stderr,'E_z(inc): ezIncInit() must be called before ezInc.\n Points per wavelength must be positive.\n');
%                 exit(-1);
%             end

    source = SourceFunction(time,sourceStr,obj.CourantNumber(1));
    switch obj.Dimensionality
        case 1
            if sourceStr.isAdd
                obj.Ez(Position(1)) = obj.Ez(Position(1)) + source;
            else
                obj.Ez(Position(1)) = source;
            end
            fprintf('%d, %g\n',time,obj.Ez(Position(1)))
        case 2
            if sourceStr.isAdd
                obj.Ez(Position(1),Position(2)) = obj.Ez(Position(1),Position(2)) + source(1);
            else
                obj.Ez(Position(1),Position(2)) = source(1);
            end
            fprintf('%d, %g\n',time,obj.Ez(Position(1),Position(2)))
        case 3
            if sourceStr.isAdd
                obj.Ez(Position(1),Position(2),Position(3)) = obj.Ez(Position(1),Position(2),Position(3)) + source(1);
            else
                obj.Ez(Position(1),Position(2),Position(3)) = source(1);
            end
            fprintf('%d, %g\n',time,obj.Ez(Position(1),Position(2),Position(3)))
    end
end
