function draw_box(top_left, bottom_right, line_width)

% draw_box.m
%   Draws an empty black square in axes. From top left to bottom right.
%

top = top_left(2);
left = top_left(1);
bottom = bottom_right(2);
right = bottom_right(1);


plot([left right], [top top], 'k', 'linewidth', line_width,...
    'Handlevisibility', 'off'), hold on
plot([right right], [top bottom], 'k', 'linewidth', line_width,...
    'HandleVisibility', 'off')
plot([left right], [bottom bottom], 'k', 'linewidth', line_width,...
    'HandleVisibility', 'off')
plot([left left], [bottom top], 'k', 'linewidth', line_width,...
    'HandleVisibility', 'off')


