function obj = add_Hy_inc(obj,Position,time,sourceStr)
    global imp0
    source = SourceFunction(time,sourceStr,obj.CourantNumber(1)) / imp0;
    switch obj.Dimensionality
        case 1
            if sourceStr.isAdd
                obj.Hy(Position(1)) = obj.Hy(Position(1)) - source;
            else
                obj.Hy(Position(1)) = -source;
            end
            fprintf('%d, %g\n',time,obj.Hy(Position(1)))
        case 2
            if sourceStr.isAdd
                obj.Hy(Position(1),Position(2)) = obj.Hy(Position(1),Position(2)) - source;
            else
                obj.Hy(Position(1),Position(2)) = -source;
            end
            fprintf('%d, %g\n',time,obj.Hy(Position(1),Position(2)))
        case 3
            if sourceStr.isAdd
                obj.Hy(Position(1),Position(2),Position(3)) = obj.Hy(Position(1),Position(2),Position(3)) - source;
            else
                obj.Hy(Position(1),Position(2),Position(3)) = -source;
            end
            fprintf('%d, %g\n',time,obj.Hy(Position(1),Position(2),Position(3)))
    end
end
