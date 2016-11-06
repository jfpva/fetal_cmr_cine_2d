function xtq = kspace_zeropad_interp( xt, xtqDim )

xt2kt = @( xt ) fftshift( fftshift( fft2( xt ), 2 ), 1 );
kt2xt = @( kt ) ifft2( ifftshift( ifftshift( kt, 2 ), 1 ) );

xtDim  = size( xt ); 

padDim  = round( ( xtqDim - xtDim ) / 2 );

kt    = xt2kt( xt );  

ktq = padarray( kt, padDim );

xtq = kt2xt( ktq );
   

end  % kspace_zeropad_interp(...)