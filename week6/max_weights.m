In = 10;
H = 36;

f = @(t) -(In*(t-1) + t*(H - t - 1) + H - t);

[x, xval] = fminunc(f, 0);
x
-xval