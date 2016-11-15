function [two_rect_horiz, two_rect_vert] = HaarF(I,y,x,length,b2)
% determine the Haar feature at the position (x,y) with a window length of
% length

    if(length == 1)
        %feature has to be determined only for one pixle

        %two rectangle feature, horizontal
        temp_1 = I((x-b2):x,(y-b2):(y+b2));
        temp_2 = I(x:(x+b2),(y-b2):(y+b2));
        two_rect_horiz = sum(temp_1(:)) - sum(temp_2(:));

        %two rectangle feature, vertical
        temp_1 = I((x-12):(x+b2),(y-b2):(y));
        temp_2 = I((x-12):(x+b2),(y):(y+b2));
        two_rect_vert = sum(temp_1(:)) - sum(temp_2(:));

    else
        %the case when the feature should determine for more than one pixle


    end
end