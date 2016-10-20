function [ A, tform2param, param2tform ] = init_tform_struct( nFrame )
%INIT_TFORM_STRUCT  initialise image sequence transformation structure
%
%   [ A, tform2param, param2tform ] = INIT_TFORM_STRUCT( nFrame )
%
%   See also: rreg_imseq2imseq, transform_imseq

% jfpva (joshua.vanamerom@kcl.ac.uk) 

tform2param = @(A) [ A.T(3,1), A.T(3,2), asind(A.T(1,2)) ];  % NOTE: for reference, http://uk.mathworks.com/discovery/affine-transformation.html
param2tform = @(tx,ty,rz) affine2d( [cosd(rz), sind(rz), 0; -sind(rz), cosd(rz), 0; tx, ty, 1 ] );

A = struct([]);

for iF = 1:nFrame,
    A(iF).A  = affine2d( [ 1 0 0 ; ...     % null transform
                0 1 0 ; ...
                0 0 1 ] );
	p = tform2param( A(iF).A );
    A(iF).tx = p(1);
    A(iF).ty = p(2);
    A(iF).rz = p(3);
end

end  % init_tform_struct(...)