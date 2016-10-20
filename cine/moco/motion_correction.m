function [ A, dispMap, U, V, imTgt, imMvg, imFix ] = motion_correction( imSrc, trrSrc, dtSrc, rrInterval, mask, pixdim, varargin )
%MOTION_CORRECTION  rigid registration of fetal heart in 2d+t image sequence
%
%       imSrc           source image sequence
%       trrSrc          time since last cardiac trigger in milliseconds for 
%                       each frame in imSrc
%       dtSrc           temporal resolution of imSrc in milliseconds
%       rrInterval      R-R interval duration in milliseconds
%       mask            TODO: add description
%       pixdim          TODO: add description
%
%       A               TODO: add description
%       dispMap         TODO: add description
%       U               TODO: add description
%       V               TODO: add description
% 
%   Additional name-value input: 
%
%       transformations initial tform struct; default identity 
%       voxelprob       voxel-wise posterior probability maps
%       frameprob       frame-wise posterior probability vector
%       verbose         verbosity; default false
%
%   See also: rreg_imseq2imseq

% jfpva (joshua.vanamerom@kcl.ac.uk)  


%% Notes

% TOOD: allow specification of origin instead of using centroid of mask


%% Parse Input


default.A0          = struct([]);                 % initial transformations
default.voxProb     = ones( size( imSrc ) );      % voxel probabilities
default.frmProb     = ones( size( imSrc, 3), 1 ); % frame probabilities
default.isVerbose   = false;

p = inputParser;

addRequired(  p, 'imSrc', @(x) validateattributes( x, {'numeric'}, ...
        {'ndims',3}, mfilename ) );
    
addRequired(  p, 'trrSrc', @(x) validateattributes( x, {'numeric'}, ...
        {'real','vector','numel',size(imSrc,3)}, mfilename ) );

addRequired(  p, 'rrInterval', @(x) validateattributes( x, {'numeric'}, ...
        {'real','scalar','positive'}, mfilename ) );

addRequired(  p, 'dtRlt', @(x) validateattributes( x, {'numeric'}, ...
        {'real','scalar','positive'}, mfilename ) );

addRequired( p, 'mask', @(x) validateattributes( x, {'logical'}, ...
        {'size',[size(imSrc(:,:,1)),NaN]}, mfilename ) );    

addRequired( p, 'pixdim', @(x) validateattributes( x, {'numeric'}, ...
        {'positive','size',[1,2]}, mfilename ) );

addParameter( p, 'transformations', default.A0, ...
        @(x) validateattributes( x, {'struct'}, ...
        {'numel',size(imSrc,3)}, mfilename ) );

addParameter( p, 'voxelprob', default.voxProb, ...
        @(x) validateattributes( x, {'numeric'}, ...
        {'real','size',size(imSrc),'<=',1,'>=',0}, mfilename ) );

addParameter( p, 'frameprob', default.frmProb, ...
        @(x) validateattributes( x, {'numeric'}, ...
        {'real','vector','numel',size(imSrc,3),'<=',1,'>=',0}, mfilename ) );

addParameter( p, 'verbose',     default.isVerbose, ...
        @(x) validateattributes( x, {'logical'}, ...
        {}, mfilename ) );

parse( p, imSrc, trrSrc, dtSrc, rrInterval, mask, pixdim, varargin{:} );

A0          = p.Results.transformations;
voxProb     = p.Results.voxelprob;
frmProb     = p.Results.frameprob;
isVerbose   = p.Results.verbose;


%% Setup


nFrame = size(imSrc,3);


%% Initialise


% Initial Transformations

if isempty( A0 ),
    A0          = init_tform_struct( nFrame );
    imSrcA0     = imSrc;
    voxProbA0   = voxProb;
else
    imSrcA0     = transform_imseq( imSrc,   A0, mask, pixdim );
    voxProbA0   = transform_imseq( voxProb, A0, mask, pixdim, 'forward', 'linear' );
end


%% Register Source Image Sequence to Target


% Target Images

imTgt = imseq_kernel_smooth( trrSrc, imSrcA0, trrSrc, 'tPeriod', rrInterval, 'kSigma', dtSrc, 'kWidth', rrInterval, 'vWeight', voxProbA0, 'fWeight', frmProb, 'fExclude', num2cell(1:nFrame) );
 

% Rigid Registration
    
U = rreg_imseq2imseq( imSrc, imTgt, mask, pixdim );
    

%% Correct for 'Overfitting' (torision of heart across cardiac cycle)
    

% NOTE: correction may not work well using outlier rejection, check
% efficacy using outlier rejection for moco iter. 2+


% Source (Moving) Images

imSrcU = transform_imseq( imSrc, U, mask, pixdim );

if ( all( voxProb(:) == 1 ) && all( frmProb(:) == 1 ) ),
    
    imMvg = imseq_kernel_smooth( trrSrc, imSrcU, trrSrc, 'tPeriod', rrInterval, 'kSigma', dtSrc, 'kWidth', rrInterval );

else
    
    imMvg = zeros( size( imSrc ), 'like', imSrc );

    voxProbU = transform_imseq( voxProb, U, mask, pixdim, 'forward', 'linear' );

    for iF = 1:nFrame,  % ensure full inclusion of current frame
        voxProbUiF = voxProbU;
        voxProbUiF(:,:,iF) = 1;
        frmProbUiF = frmProb;
        frmProbUiF(iF) = 1;
        imMvg(:,:,iF) = imseq_kernel_smooth( trrSrc, imSrcU, trrSrc(iF), 'tPeriod', rrInterval, 'kSigma', dtSrc, 'kWidth', rrInterval, 'vWeight', voxProbUiF, 'fWeight', frmProbUiF );
    end
    
end


% Target (Fixed) Images

% if ( all( voxProb(:) == 1 ) && all( frmProb(:) == 1 ) ),

    imFix = imseq_kernel_smooth( trrSrc, imSrc,  trrSrc, 'tPeriod', rrInterval, 'kSigma', dtSrc, 'kWidth', rrInterval );  

% else
%     
%     imFix = zeros( size( imSrc ), 'like', imSrc );
%     
%     for iF = 1:nFrame,  % ensure full inclusion of current frame
%         voxProb0iF = voxProb;
%         voxProb0iF(:,:,iF) = 1;  
%         frmProb0iF = frmProb;
%         frmProb0iF(iF) = 1;
%         imFix(:,:,iF) = imseq_kernel_smooth( trrSrc, imSrc, trrSrc(iF), 'tPeriod', rrInterval, 'kSigma', dtSrc, 'kWidth', rrInterval, 'vWeight', voxProb0iF, 'fWeight', frmProb0iF );
%     end
%     
% end


% Rigid Registration

V = rreg_imseq2imseq( imMvg, imFix, mask, pixdim );


%% Combine Transformations

[ A, tform2param ] = init_tform_struct( nFrame );

for iF = 1:nFrame,
    A(iF).A.T = U(iF).A.T * V(iF).A.T;
    p = tform2param( A(iF).A );
    A(iF).tx = p(1);
    A(iF).ty = p(2);
    A(iF).rz = p(3);
end


%% Calculate Displacement

dispMap = zeros( size( imSrc ) );         

R = make_imref2d( imSrc(:,:,1), mask, pixdim ); 

for iF = 1:nFrame;
    
    dispMap(:,:,iF) = calc_dispmap_from_affine2d( A(iF).A, R );

end


%% Visualise

if ( isVerbose ),
    
    cMap = lines(3);
    
    k = 2/3;
    
    figure('Name',sprintf('motion_correction'))
    
    plot([A.tx],'-','Color',cMap(1,:),'LineWidth',2),
    hold on,
    plot([A.ty],'-','Color',cMap(2,:),'LineWidth',2),
    plot([A.rz],'-','Color',cMap(3,:),'LineWidth',2),
    plot([U.tx],'-.','Color',k*cMap(1,:),'LineWidth',1),
    plot([U.ty],'-.','Color',k*cMap(2,:),'LineWidth',1),
    plot([U.rz],'-.','Color',k*cMap(3,:),'LineWidth',1),
    plot([V.tx],'-','Color',k*cMap(1,:),'LineWidth',1),
    plot([V.ty],'-','Color',k*cMap(2,:),'LineWidth',1),
    plot([V.rz],'-','Color',k*cMap(3,:),'LineWidth',1),
    grid on,
    legend( 'A tx','A ty','A rz', 'U tx','U ty','U rz', 'V tx','V ty','V rz' )
    xlabel('real-time frame no.')
    ylabel('displacement [mm] / rotation [deg.]')
    set(gca,'YLim',max(abs(get(gca,'YLim')))*[-1,+1])
    
end


end  % motion_correction(...)
