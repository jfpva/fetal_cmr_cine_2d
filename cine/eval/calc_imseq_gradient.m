function [g,gx,gy,gt] = calc_imseq_gradient( imseq, x, y, t, rrInterval )

x = reshape( x, [], 1, 1 );
y = reshape( y, 1, [], 1 );
t = reshape( t, 1, 1, [] );

if numel( x ) == 1,
    dx = x;
else
    dx = diff( x, 1, 1 ); 
    dx(end+1) = dx(end);
end

if numel( y ) == 1,
    dy = y;
else
    dy = diff( y, 1, 2 ); 
    dy(end+1) = dy(end);
end

if numel( t ) == 1,
    dt = t;
else
    dt = diff( cat( 3, t, t(1) + rrInterval ), 1, 3 ); 
end

gx = bsxfun( @rdivide, diff( padarray( imseq, [1,0,0], 'replicate', 'post'), 1, 1 ), dx );

gy = bsxfun( @rdivide, diff( padarray( imseq, [0,1,0], 'replicate', 'post'), 1, 2 ), dy );

gt = bsxfun( @rdivide, diff( padarray( imseq, [0,0,1], 'circular',  'post'), 1, 3 ), dt );

g = sqrt( abs(gx).^2 + abs(gy).^2 + abs(gt).^2 );


end