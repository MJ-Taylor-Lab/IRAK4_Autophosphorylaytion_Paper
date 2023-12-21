img  = zeros(10, 10);
x    = 6;   % Position of center
y    = 4;
r    = 3;   % Radius
[h, w] = size(img);
mask = ((1-x:h-x).' .^2 + (1-y:w-y) .^2) <= r^2;
disp(mask)
img(mask) = 2;
disp(img)
Result    = sum(img(mask), 'all')