function [count,p,coords] = seqdotplotsilent(seq1,seq2,varargin)
% SEQDOTPLOTSILENT generates a dotplot of sequence matches.
%
%   SEQDOTPLOTSILENT(S,T) plots the sequence matches of sequences S and T.
%
%   SEQDOTPLOTSILENT(S,T,WINDOW,NUM) plots sequence matches when there are at
%   least NUM matches in a window of size WINDOW. For nucleotide sequences
%   the literature recommends a WINDOW of 11 and NUM of 7.
%
%   MATCHES = SEQDOTPLOTSILENT(...) returns the number of dots in the dotplot
%   matrix.
%
%   [MATCHES, MATRIX] = SEQDOTPLOTSILENT(...) returns the dotplot as a sparse
%   matrix.
%
%   SEQDOTPLOTSILENT(...,'DISPLAY',true) turns the verbose mode on. Useful for
%   large queries where the algorithm takes large amounts of memory and
%   time.
%
%   Example:
%
%       moufflon = getgenbank('AB060288','sequence',true)
%       takin = getgenbank('AB060290','sequence',true)
%       seqdotplot(moufflon,takin,11,7)
%
%   This shows the similarities between prion protein (PrP) nucleotide
%   sequences of two ruminants, the moufflon and the golden takin.
%
%   See also NWALIGN, SWALIGN.

% 	Example reference:
% 	Comparative analysis of the prion protein open reading frame nucleotide
% 	sequences of two wild ruminants, the moufflon and golden takin.
%
% 	Seo SW, Hara K, Kubosaki A, Nasu Y, Nishimura T, Saeki K, Matsumoto Y,
% 	Endo H, Onodera T.

% Copyright 2003-2004 The MathWorks, Inc.
% $Revision: 1.1 $   $Date: 2006/05/29 14:27:44 $

path
%%% setting some constants
memoryBlock    = 2^26;
limitForSparse = 2^28;  % (in bytes, i.e. can contain limitForSparse/4 matches)

%%% setting defaults for parameters
verbosity = false;
window = 1;
stringency =1;

%%% get space limit for image in pixels
rootUnits = get(0, 'Units');
set(0, 'Units', 'pixels');
screnSize = get(0,'ScreenSize')*[0 0 1 0;0 0 0 1]';
windowLimits = screnSize - [20 110];
windowLimits = round( windowLimits * 0.85);
set(0, 'Units', rootUnits);

%%% processing inputs
switch nargin
    case 0,
        seq1 = randseq(15);
        seq2 = seq1;
        fprintf('Creating some random data... \n')
    case 1,
        seq2 = seq1;
    otherwise
        if nargin > 2
            val = varargin{1};
            if isreal(val) && numel(val)==1 && isnumeric(val)  
                window = val;
                varargin(1) = [];
            end
        end
        if nargin > 3
            val = varargin{1};
            if isreal(val) && numel(val)==1 && isnumeric(val)   
                stringency = val;
                varargin(1) = [];
            else
                stringency = window;
            end    
        end
end

% processing remaining varargins
if  numel(varargin)
    if rem(numel(varargin),2) == 1
        error('Bioinfo:IncorrectNumberOfArguments',...
            'Incorrect number of arguments to %s.',mfilename);
    end
    okargs = {'display'};
    for j=1:2:numel(varargin)
        pname = varargin{j};
        pval = varargin{j+1};
        k = strmatch(lower(pname), okargs); %#ok
        if isempty(k)
            error('Bioinfo:UnknownParameterName',...
                'Unknown parameter name: %s.',pname);
        elseif length(k)>1
            error('Bioinfo:AmbiguousParameterName',...
                'Ambiguous parameter name: %s.',pname);
        else
            switch(k)
                case 1  % verbosity on/off
                    verbosity = pval == true;
            end
        end
    end
end

% If the input is a structure then extract the Sequence data.
if isstruct(seq1)
    try
        seq1 = seqfromstruct(seq1);
    catch
        rethrow(lasterror);
    end
end
if isstruct(seq2)
    try
        seq2 = seqfromstruct(seq2);
    catch
        rethrow(lasterror);
    end
end

% testing valid conditions of input parameters
if stringency > window
  error('Bioinfo:DotplotTooStringent',...
  'The number of matches per window cannot be more than the window size.');
end

% convert to integers and check that both sequences are of the same type
if (ischar(seq1) && ischar(seq2))
    if (isnt(seq1) && isnt(seq2))
        iSeq1 = nt2int(seq1);
        iSeq2 = nt2int(seq2);
    elseif (isaa(seq1) && isaa(seq2))
        iSeq1 = aa2int(seq1);
        iSeq2 = aa2int(seq2);
    else
      error('Bioinfo:DotplotMismatch',...
      'Both sequences must be valid sequences of the same type (AA or NT).');
    end
elseif (isnumeric(seq1) && isnumeric(seq2))
    if (isnt(seq1) && isnt(seq2))
        iSeq1 = seq1;
        iSeq2 = seq2;
    elseif (isaa(seq1) && isaa(seq2))
        iSeq1 = seq1;
        iSeq2 = seq2;
    else
      error('Bioinfo:DotplotMismatch',...
      'Both sequences must be valid sequences of the same type (AA or NT).');
    end
else
   error('Bioinfo:DotplotMismatch',...
   'Both sequences must be the same type (alpha symbol or numeric array).');
end

% length of input sequences
le1=length(iSeq1); 
le2=length(iSeq2);

if le1>le2 &&  le1 < 0
    switchedSequences = true;
    tem = iSeq1;
    iSeq1 = iSeq2;
    iSeq2 = tem;
    le1=length(iSeq1); 
    le2=length(iSeq2);
else
    switchedSequences = false;
end    


% pad sequences to contain window at the end
padSeq = zeros(1,window-1);
iSeq1 = [iSeq1 padSeq];
iSeq2 = [iSeq2 padSeq];

mapSize = double(max([iSeq1(:);iSeq2(:)]));

% we can not index 0's so we re-map to the last position of map
iSeq1(iSeq1 == 0) = mapSize+1;
iSeq2(iSeq2 == 0) = mapSize+1;
 
% create map, we do not want zeros to match zeros (last position in map)
map = diag(uint8([true(1,mapSize),false]));
%memoryBlock = le1*le2;
columnBlock = floor(memoryBlock/le1); %including window columns

p=uint32([]);
count = 0;

s = zeros(le1,columnBlock-window+1,'uint8');
firstCols =  1:columnBlock-window+1:le2;
if (numel(firstCols) > 1) && verbosity 
   t0 = clock;
   fprintf('Search is divided in %d memory blocks \n', numel(firstCols))
end

for i = 1:numel(firstCols)
    j = firstCols(i);
    % compute some constants for this block
    lastCol = min( le2+window-1 , j+columnBlock-1);
    numCols = lastCol - j + 1 - window + 1;
    colInd  = j:lastCol;

    mask = map(iSeq1,iSeq2(colInd));

    if (numCols+window-1 == columnBlock)
        s(:)=uint8(0);
        for k = 1:window
            s = s + mask(k:le1+k-1,k:numCols+k-1);
        end
        s = s .* mask(1:le1,1:numCols);
        h = find(s >= stringency);
        count = count + numel(h);
        if (4*count) > limitForSparse
          error('Bioinfo:LimitExceeded',['Match matrix exceeded allowed'...
                ' limit.\n    Try to reduce the sequence lengths or'...
                ' increase\n    the window stringency parameters']);
        end
        p = [p;uint32((j-1)*le1+h)];
    else % when last block is incomplete
        s = zeros(le1,numCols,'uint8');
        
        for k = 1:window
            s = s + mask(k:le1+k-1,k:numCols+k-1);
        end
        s = s .* mask(1:le1,1:numCols);
%	fprintf('size of s: \n');
%	size(s)
        h = find(s >= stringency);
	[hr, hc] = find(s >= stringency);
        count = count + numel(h);
        if (4*count) > limitForSparse
          error('Bioinfo:LimitExceeded',['Match matrix exceeded allowed'...
                ' limit.\n    Try to reduce the sequence lengths or'...
                ' increase\n    the window stringency parameters']);
        end
        p = [p;uint32((j-1)*le1+h)];
%	fprintf('hsize: \n');
    end
    if (numel(firstCols) > 1) && verbosity 
      fprintf('Memory block: %d            Elapsed time: %f seconds \n',...
               i,etime(clock,t0))  
      fprintf('Match matrix size: %d  Memory used: %d bytes\n',...
               count,count*4)  
    end
end

% match matrix info should be in sparse, as it is needed fot the output,
% but also for easy downsampling computation
p=reshape(sparse(double(p),1,true,le1*le2,1),le1,le2);
nc = size(hc,1);
coords(1:nc,1) = hr;
coords(1:nc,2) = hc;
coords(1:nc,3) = hr + window;
coords(1:nc,4) = hc + window;

return;
%hFig = figure('tag','seqDotPlot',...
%              'nextplot','replacechildren');

%set(hFig,'Visible','on')


imagesc(full(p));
colormap(1-gray);


% computing the scalign factors in the case the match matrix does not fit in
% the screen
downSampleX = ceil( le2 / windowLimits(1) );
downSampleY = ceil( le2 / windowLimits(2) );
downSampleX = 1;
downSampleY = 1;

if (downSampleX > 1) || (downSampleY > 1)
%if (1 == 0)
    warning('Bioinfo:imageTooBigForScreen',...
        ['Match matrix has more points than available screen pixels.'...
        '\n         Scaling image by factors of %d in X and %d in Y'],...
        downSampleX,downSampleY)

    
    imageSizePixelsX = ceil(size(p,2)/downSampleX);
    imageSizePixelsY = ceil(size(p,1)/downSampleY);
    
    % setting the size to a multiple of downSample for easy downsampling,
    % (like padding with zeros a full matrix)
    padCordinatesX = imageSizePixelsX*downSampleX;
    padCordinatesY = imageSizePixelsY*downSampleY;
    
    if ~isequal(size(p),[padCordinatesY,padCordinatesX]) % do not pad if p 
                                                     % is already that size
        p(padCordinatesY,padCordinatesX)=false;
    end
    
    % setting the size and class of the x-downsampled matrix
    tp = sparse([],[],false,padCordinatesY,imageSizePixelsX);
    
    % x-downsampling in the sparse domain
    for k = 1:downSampleX
        tp = tp | p(:,k:downSampleX:padCordinatesX);
    end
    
    % setting the size and class of the output Image
    I = false(imageSizePixelsY,imageSizePixelsX);
    
    % y-downsampling 
    I(ceil(find(tp(:))/downSampleY))=true;
    
    % resize the padded p so the output has the sequence lengths
    p = p(1:le1,1:le2);
    fprintf('assigning p to tp\n');
    p = tp;
else % the image fits in the screen, no need to downsample
  fprintf('not downsampling\n');
    I =full(p);
end
%I = flipdim(I,1);
if switchedSequences
    xLabel = 'Sequence 1';
    yLabel = 'Sequence 2';
else
    xLabel = 'Sequence 2';
    yLabel = 'Sequence 1';
end

hFig = figure('tag','seqDotPlot',...
              'nextplot','replacechildren');


%xData = (0:size(I,2)-1)*downSampleX+1;
%yData = (0:size(I,1)-1)*downSampleY+1;       
%
%xDataMax = max(xData);
%yDataMax = max(yData);
%imagesc(xData, yDataMax-yData, I,'xData',xData,'yData',yData);


hAxis = gca;
axis xy;
axis(hAxis,'image'); 
set(hFig,'Visible','on')


if (downSampleX==1 && downSampleY==1) 
    axis(hAxis,'image'); 
end

%set(hAxis, 'YDir',  'reverse');
%set(hAxis, 'XDir',  'reverse');
set(hAxis,'XAxisLocation','bottom','Units','pixels' );
xlabel(xLabel,'fontsize',12);
ylabel(yLabel,'fontsize',12);

windowSize = ceil(size(I)/.85);

if any(windowSize>get(hAxis,'Position')*[0 0 0 1;0 0 1 0]')
    figPos([1 3]) = [ceil((screnSize(1)-windowSize(2))/2) windowSize(2)];
    figPos([2 4]) = [screnSize(2) - 90 - windowSize(1) windowSize(1)];
    set(hFig,'Position',figPos)

    axesPos([2 1]) = ceil(windowSize.*[0.05 0.1])+.5 ;
    axesPos([4 3]) = size(I);
    axesPos(2) = axesPos(2) + max(0,floor((windowSize(1)*0.85 - size(I,1))/2));
    set(hAxis,'Position',axesPos)
    set(hAxis,'Units','Normalized')
else % image fits do not need to resize
    set(hAxis,'Units','Normalized')
    pos = get(hAxis,'position');
    set(hAxis,'position',pos - [0 .05 0 0]);
    % for short sequences, force the tick marks to be integer values.
    if le1<6  
        set(hAxis,'ytick',1:le1); 
    end
    if le2<6  
        set(hAxis,'xtick',1:le2); 
    end
end

set(hFig,'Visible','on')

% with output args force everything to be doubles, clear not used outputs
switch nargout
    case 2; 
        if switchedSequences 
            p = double(p');
        else
            p = double(p);
        end
    case 1; clear p;
    otherwise; clear p; clear count;
end

p = I;
