function a = ChooseAction(x)
 r = rand;

    switch x
        case 0 
            if(r<0.5)
                a = 0;
            else
                a = 4;
            end
        case 1
            if(r<(1/3))
                a = 1;
            elseif(r<(2/3))
                a = 3;
            else
                a = 5;
            end
        case 2
            if(r<0.5)
                a = 2;
            else
                a = 3;
            end
        case 3
            if(r<0.25)
                a = 1;
            elseif(r<0.5)
                a = 2;
            elseif(r<0.75)
                a = 3;
            else
                a = 4;
            end
        case 4
            if(r<0.25)
                a = 0;
            elseif(r<0.5)
                a = 3;
            elseif(r<0.75)
                a = 4;
            else
                a = 5;
            end
        case 5
            if(r<(1/3))
                a = 1;
            elseif(r<(2/3))
                a = 4;
            else
                a = 5;
            end
    end