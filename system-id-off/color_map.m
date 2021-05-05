function colorMap = color_map(color1, color2, nShades, showColors)

% color_map.m
%   Creates a colormap matrix to be used in plots.
%
%
%   Syntax:
%       colorMap = color_map(color1, color2, nShades)
%
%   Input:
%       - color1: RGB row vector of color for highest values.
%       - color2: RGB row vector of color for lowest values.
%       - nShades: number of shades inbetween color1 and color2
%       - showColors: boolean, true if want to plot colorMap
%
%   Output
%       - colorMap: matrix containing both colors and the shades inbetween.
%
%   Author:
%       - Luiz Felipe da S. Coelho -- luizfelipe.coelho@smt.ufrj.br
%


colorMap = zeros(nShades+2, 3);
colorMap(1, :) = color2;
colorMap(end, :) = color1;

differenceR = abs(color1(1)-color2(1));
differenceG = abs(color1(2)-color2(2));
differenceB = abs(color1(3)-color2(3));

aux = ones(1, 3);
if color1(1) < color2(1)
    aux(1) = -1*differenceR./(nShades+1);
else
    aux(1) = differenceR./(nShades+1);
end
if color1(2) < color2(2)
    aux(2) = -1*differenceG./(nShades+1);
else
    aux(2) = differenceG./(nShades+1);
end
if color1(3) < color2(3)
    aux(3) = -1*differenceB./(nShades+1);
else
    aux(3) = differenceB./(nShades+1);
end

for i = 2:nShades+1
    colorMap(i, :) = colorMap(i-1, :) + aux;
end

switch showColors
    case true
        nRows = 100;
        nColumns = 500;
        imgMatR = ones(nRows, nColumns);
        imgMatG = ones(nRows, nColumns);
        imgMatB = ones(nRows, nColumns);
        nColorBlocks = nColumns./(nShades+2);
        for i = 1:nColumns
            imgMatR(:, i) = imgMatR(:, i).*colorMap(ceil(i./nColorBlocks), 1);
            imgMatG(:, i) = imgMatG(:, i).*colorMap(ceil(i./nColorBlocks), 2);
            imgMatB(:, i) = imgMatB(:, i).*colorMap(ceil(i./nColorBlocks), 3);
        end
        image = zeros(nRows, nColumns, 3);
        image(:, :, 1) = imgMatR;
        image(:, :, 2) = imgMatG;
        image(:, :, 3) = imgMatB;
        figure('Name', 'ColorMap')
        imshow(image)
    otherwise
end


% EoF

