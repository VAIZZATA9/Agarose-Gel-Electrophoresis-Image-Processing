function pos = boxFinder(im)
    threshold = 0.51;
    imBW = im;
    pos = [1,1,1,1];

    for i = 1:size(im,1);
        for j = 1:size(im,2);
            if im(i,j) < threshold
                imBW(i,j) = 0;
            else
                imBW(i,j) = 1;
            end
        end
    end

    for i = 1:size(im,1);
        for j = 1:size(im,2);
            width = 1;
            height = 1;
            xpos = j;
            ypos = i;
            while imBW(i,j) == 1 && i < size(im,1)
                i = i + 1;
                height = height + 1;
            end
            i = ypos;
            while imBW(i,j) == 1 && j < size(im,2)
                j = j + 1;
                width = width + 1;
            end
            if width * height > pos(3) * pos(4);
                pos = [ypos,xpos,height,width];
            end
        end
    end
end