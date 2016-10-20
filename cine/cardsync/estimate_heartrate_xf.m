function [ rrInterval, triggerTime ] = estimate_heartrate_xf( imSeq, frameDuration, varargin )
%ESTIMATE_HEARTRATE_XF   Estimate heartrate in x-f space.
%   [ rrInterval ] = ESTIMATE_HEARTRATE_XF( imSeq, frameDuration ) 
%   returns an estimate of R-R interval for 2D+time array of realtime 
%   images, imSeq, with temporal resolution, frameDuration.
%   Heartrate is estimated using peaks in x-f space at fundamental  
%   frequency and first harmonic.
%   All times are in units of seconds.
%
%   [ ..., triggerTime ] = ESTIMATE_HEARTRATE_XF( ... ) also returns
%   estimated cardiac trigger times.
%
%   ESTIMATE_HEARTRATE_XF( ..., 'roi', roi ) uses 2D logical array roi as a 
%   mask; default is full FOV ROI 
% 
%   ESTIMATE_HEARTRATE_XF( ..., 'useHarmonic', false ) look for peaks at
%   fundamental frequencies only; default is true
%
%   ESTIMATE_HEARTRATE_XF( ..., 'useUpsampleFreq', false ) upsample
%   frequency resolution;  default is true
%
%   ESTIMATE_HEARTRATE_XF( ..., 'hrRange', [minHR maxHR] ) look for peaks at
%   fundamental frequencies between min and max HR; default is [120 160]
% 
%   ESTIMATE_HEARTRATE_XF( ..., 'verbose', true ) shows verbose output; 
%   default is false 
%
%   If imSeq has a fourth dimension that is non-singleton, multi-channel
%   images are assumed and combined using root sum-of-squares.
% 

%   jfpva (joshua.vanamerom@kcl.ac.uk)


%% NOTES


% TODO: improve peak detection 


%% Parse Inputs


p = inputParser;

default.roi             = [];
default.useHarmonic     = true;
default.useUpsampleFreq = true;
default.hrRange         = [120 160];  % bpm
default.isVerbose       = false;

addRequired(  p, 'imSeq', ...
    @(x) validateattributes( x, {'numeric'}, {'size',[NaN NaN NaN NaN 1]},mfilename) );
addRequired(  p, 'frameDuration',  ...
    @(x) validateattributes( x, {'numeric'}, {'scalar','positive'}, mfilename) );
addParameter(  p, 'roi', default.roi, ...
    @(x) validateattributes( x, {'logical'}, {'2d'}, mfilename));
addParameter( p, 'useHarmonic',  default.useHarmonic, ...
    @(x) validateattributes( x, {'logical'}, {'scalar'}, mfilename) );
addParameter( p, 'useUpsampleFreq',  default.useUpsampleFreq, ...
    @(x) validateattributes( x, {'logical'}, {'scalar'}, mfilename) );
addParameter( p, 'hrRange',  default.hrRange, ...
    @(x) validateattributes( x, {'numeric'}, {'size',[1 2]}, mfilename) );
addParameter( p, 'verbose',  default.isVerbose, ...
    @(x) validateattributes( x, {'logical'}, {'scalar'}, mfilename) );

parse( p, imSeq, frameDuration, varargin{:} );

roi             = p.Results.roi;
useHarmonic     = p.Results.useHarmonic;
useUpsampleFreq = p.Results.useUpsampleFreq;
hrRange         = p.Results.hrRange;
isVerbose       = p.Results.verbose;


%% Setup


% Dimensions

[nX,nY,nDyn,nChan] = size(imSeq);


% Realtime Image Sequence

if nChan > 1,  % ensure single-channel images
    imSeq = sqrt( sum( imSeq.^2, 4 ) );
end


% Timing

minHR = min( hrRange );  % bpm
maxHR = max( hrRange );  % bpm

minRR = 60/maxHR;  % s
maxRR = 60/minHR;  % s

minFreq = 1/maxRR; % Hz
maxFreq = 1/minRR; % Hz

if useUpsampleFreq,
    dt = 0.0001;                   % target temporal resolution for estimation of heartrate
    df = 1/(maxRR) - 1/(maxRR+dt); % frequency resolution at min expected fundamental frequency
    nF = range(calc_freq( nDyn, frameDuration )) / df;  % number of x-f frames to reconstruct
else
    nF = nDyn;
end


% Check against ROI

if isempty( roi ),
    
    warning('No ROI specified, using full FOV'),
    
    roi = true(nX,nY);
    
    nF = nDyn;  % don't upscale frequency resolution for full-FOV ROI since
                % reconstructing many more frames in x-f space for full-FOV 
                % ROI is computationally demanding 

else
    
    nVox = sum(roi(:));
    
    while nF*nVox > 1e8,
        nF = nF/10;
    end
    
end


% x-t to x-f

xt2xf = @( xt ) fftshift( fft( xt, [], 3 ), 3 );



%% Get Image Signal in ROI


xtRoi = nan( sum(roi(:)), 1, nDyn );

for iF = 1:nDyn,
    imFrame = imSeq(:,:,iF);
    xtRoi(:,:,iF) = imFrame(roi);
end


%% Get x-f Magnitude in ROI


xfRoi = abs( xt2xf( padarray( xtRoi, [0,0,ceil((nF-nDyn)/2)] ) ) );

xfMeanSig = nan(1,size(xfRoi,3));

for iF = 1:length(xfMeanSig),
    xfMeanSig(iF) = mean( xfRoi(:,:,iF) );
end

xfMeanSig = smooth( xfMeanSig, nF/nDyn, 'moving' );  % remove ringing introdcued by zero-padding x-t space in time to upscale frequency resolution


%% Identify Fundamental Frequency


f = calc_freq( size(xfRoi,3), frameDuration );

minPeakDistance   = find(f>minFreq,1,'first') - find(f==0,1,'first');
indf0 = find(f<minFreq,1,'last'):find(f>maxFreq,1,'first');
minPeakProminence = prctile(xfMeanSig(indf0),100)-prctile(xfMeanSig(indf0),90);
                                            % WAS:  minPeakProminence = prctile(xfMeanSig(indf0),100)-prctile(xfMeanSig(indf0),67);
                                            % WAS:  minPeakProminence = prctile(xfMeanSig,50);  
                                            % NOTE: minPeakProminence = prctile(xfMeanSig,18); seems to work for HR estimation after rreg
                                            
[ ~, fPeak ] = findpeaks( xfMeanSig, 'MinPeakDistance', minPeakDistance, 'MinPeakProminence', minPeakProminence );

heartFreqPeak  = f( fPeak );


%% Identify Fundamental and Harmonic Frequencies


fundamental  = abs( heartFreqPeak( ( abs(heartFreqPeak) > minFreq ) & ( abs(heartFreqPeak) < maxFreq ) ) );
if isempty( fundamental ),
    error( '' )
elseif length( fundamental ) == 1,
    fundamental = [ fundamental, NaN ];
else
    fundamental = sort( fundamental );
    fundamental = fundamental(1:2);
end

harmonic1    = abs( heartFreqPeak( ( abs(heartFreqPeak) >= 2 * minFreq ) & ( abs(heartFreqPeak) < 2 * maxFreq ) ) / 2 );
if length( harmonic1 ) < 2,
    useHarmonic = false;
else
    harmonic1 = sort( harmonic1 );
    harmonic1 = harmonic1( 1:2 );
end


%% Calculate Heartrate


refFreq    = fundamental;

if ( useHarmonic ),
   refFreq    = [ refFreq, harmonic1 ];
end    

rrInterval = mean( 1 ./ refFreq( ~isnan( refFreq ) ) );


%% Calculate Cardiac Trigger Times


nTrigger = ceil( nDyn * frameDuration / rrInterval );

triggerTime = rrInterval * (0:nTrigger);


%% Verbose


if ( isVerbose ),
    
    figure( 'Name', 'heart_rate_estimate' );
    
    plot( f, xfMeanSig ),  % plot mean x-f magnitude in ROI v. frequency
    hold on
    plot( f( fPeak), xfMeanSig( fPeak ), 'o')
    title('Heart Rate Estimation');
    xlabel('frequency (Hz)');
    ylabel('mean-signal_{ROI} (a.u.)');
    legend('signal','peaks'),
    hAx = gca;
    hAx.YTick = 0;
    a = axis;
    text(a(1)+0.05*(a(2)-a(1)),a(4)*0.5, sprintf('mean R-R interval = %.1f ms', rrInterval * 1000 ) ),
    if ( ~isempty( fundamental ) ),
        text(a(1)+0.05*(a(2)-a(1)),a(4)*0.4, sprintf('f_{0} = %.3f, %.3f Hz => %.1f, %.1f ms', fundamental(1), fundamental(2), 1000/fundamental(1), 1000/fundamental(2) ) ),
    end
    if ( ( useHarmonic ) && ~isempty( harmonic1 ) ),
        text(a(1)+0.05*(a(2)-a(1)),a(4)*0.3, sprintf('f_{1} = %.3f, %.3f Hz => %.1f, %.1f ms', 2*harmonic1(1), 2*harmonic1(2), 1000/harmonic1(1), 1000/harmonic1(2) ) )
    end
    
    if any(~roi(:)),  % show ROI
        
        axes('Position',[.55 .3 .4 .4])
        im = abs( mean( imSeq, 3 ) );
        B = bwboundaries( roi );
        imshow( im, [ prctile(im(:),1), prctile(im(:),99)] )
        hold on
        line(B{1}(:,2),B{1}(:,1),'LineWidth',1,'Color','c')
    
    end
    
end


end  % estimate_heartrate_xf(...)