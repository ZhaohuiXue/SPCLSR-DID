function [ result_rgb ] = color_picture( result)
    result_rgb = zeros(size(result,1),size(result,2),3);
    color_rgb = [255 255 255;
        0 168 132;
        76 0 115;
        0 0 0;
        190 255 232;
        255 0 0;
        115 0 0;
        205 205 102;
        137 90 68;
        215 158 158;
        255 115 223;
        0 0 255;
        156 156 156;%
        115 223 255;%
        0 255 0;%
        255 255 0;%
        255 170 0];
    for n = 0:max(max(result))
        [hang,lie,] = find(result == n);
    
        for k = 1:size(hang)
            result_rgb(hang(k),lie(k),:) = color_rgb(n+1,:);
        end
    end
    result_rgb = uint8(result_rgb);
        
