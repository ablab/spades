function [tr,tree,names]  = adjlphytreeread(filename)
%PHYTREEREAD reads a NEWICK tree formatted file.
%
%  TREE = PHYTREEREAD(FILENAME) reads a NEWICK tree formatted file
%  FILENAME, returning the data in the file as a PHYTREE object. FILENAME
%  can also be a URL or MATLAB character array that contains the text of a
%  NEWICK format file. 
%
%  The NEWICK tree format is found at: 
%        http://evolution.genetics.washington.edu/phylip/newicktree.html
%  
%  Note: This implementation only allows binary trees, non-binary trees
%  will be translated into a binary tree with extra branches of length 0.
%
%   Example:
%
%      tr = phytreeread('pf00002.tree')
%
%   See also GETHMMTREE, PHYTREE, PHYTREETOOL, PHYTREEWRITE.

% Copyright 2003-2004 The MathWorks, Inc.
% $Revision: 1.1 $ $Author: mchaisso $ $Date: 2006/06/28 23:57:09 $


if nargin==0
     [filename, pathname] = uigetfile({'*.tree';'*.dnd'},'Select Phylogenetic Tree File');
     if ~filename
        disp('Canceled, file not read.');
        tr=[];
        return;
    end
    filename = [pathname, filename];
end
    
% check input is char
% in a future version we may accept also cells
if ~ischar(filename)
    error('Bioinfo:InvalidInput','Input must be a character array')
end

if size(filename,1)>1  % is padded string ?
    strin = cellstr(filename);
    strin = [strin{:}];
elseif (strfind(filename(1:min(10,end)), '://')) % is an url ?
    if (~usejava('jvm'))
        error('Bioinfo:NoJava','Reading from a URL requires Java.')
    end
    try
        strin = urlread(filename);
    catch
        error('Bioinfo:CannotReadURL','Cannot read URL "%s".', filename);
    end
    strin = strread(strin,'%c','delimiter','\n')';
elseif  (exist(filename,'file') || ...
         exist(fullfile(cd,filename),'file') )    %  is a valid filename ?
    strin = textread(filename,'%c','delimiter','\n')';
else  % must be single a string with '\n'
    strin = strread(filename,'%c','delimiter','\n')';
end

% characterizing the string
numBranches = sum(strin==',');
numLeaves   = numBranches + 1;
numLabels   = numBranches + numLeaves;

if (numBranches == 0)
    error('Bioinfo:NoCommaInInputString', ...
          ['There is not any comma in the data,\ninput string may not '...
           'be in Newick style or is not a valid filename.'])
end

% find the string features: open and close parentheses and leaves
leafPositions = regexp(strin,'[(,][^(,)]')+1;
parenthesisPositions = regexp(strin,'[()]');
strFeatures = strin(sort([leafPositions parenthesisPositions]));

% some consistency checking on the parenthesis
temp = cumsum((strFeatures=='(') - (strFeatures==')'));
if any(temp(1:end-1)<1) || (temp(end)~=0)
    error('Bioinfo:InconsistentParentheses','The parentheses structure is inconsistent,\ninput string may not be in Newick style or is not a valid filename.')
end

dist = zeros(numLabels,1);             % allocating space for distances
tree = zeros(numBranches,2);           % allocating space for tree pointers
names = cell(numLabels,1);             % allocating space for tree labels

try

% extract label information for the leaves
leafData = regexp(strin,'[(,][^(,);\[\]]+','match');
for j=1:numel(leafData)
    coi = find(leafData{j}==':',1,'last');
    if isempty(coi) % if no colon no length, the whole label is the name
        dist(j) = 0;
        names{j} = leafData{j}(2:end);
    else % if there is colon, get name and length
        dist(j) = strread(leafData{j}(coi+1:end),'%f');
        names{j} = leafData{j}(2:coi-1);
    end
end 
% uniformizing empty cells, value inside the brackets can never be empty
% because branch names will always be empty
[names{cellfun('isempty',names)}] = deal('');
        
% extract label information for the parenthesis
parenthesisData = regexp(strin,')[^(,);\[\]]*','match');
parenthesisDist = zeros(numel(parenthesisData),1);
for j=1:numel(parenthesisData)
    coi = find(parenthesisData{j}==':',1,'last');
    if isempty(coi) % if no colon no length, the whole label is the name
        parenthesisDist(j) = 0;
        parenthesisData{j} = parenthesisData{j}(2:end);
    else % if there is colon, get name and length
        parenthesisDist(j) = strread(parenthesisData{j}(coi+1:end),'%f');
        parenthesisData{j} = parenthesisData{j}(2:coi-1);
    end
end 
% uniformizing empty cells, value inside brackes may be empty
if any(cellfun('isempty',parenthesisData))
    [parenthesisData{cellfun('isempty',parenthesisData)}] = deal('');
end

li = 1; bi = 1; pi = 1;          % indexes for leaf, branch and parentheses
queue = zeros(1,2*numLeaves); qp = 0; % setting the queue (worst case size)

j = 1;

while j <= numel(strFeatures)
    switch strFeatures(j)
        case ')' % close parenthesis, pull values from the queue to create
                 % a new branch and push the new branch # into the queue                 
            lastOpenPar = find(queue(1:qp)==0,1,'last');     
            numElemInPar = min(3,qp-lastOpenPar);
            switch numElemInPar
                case 2  % 99% of the cases, two elements in the parenthesis
                    bp = bi + numLeaves;
                    names{bp} = parenthesisData{pi};      % set name
                    dist(bp) = parenthesisDist(pi);       % set length 
                    tree(bi,:) = queue(qp-1:qp); 
                    qp = qp - 2; % writes over the open par mark
                    queue(qp) = bp;
                    bi = bi + 1;
                    pi = pi + 1;
                case 3  % find in non-binary trees, create a phantom branch
                    bp = bi + numLeaves;
                    names{bp} = '';      % set name
                    dist(bp) = 0;        % set length 
                    tree(bi,:) = queue(qp-1:qp); 
                    qp = qp - 1; % writes over the left element
                    queue(qp) = bp;
                    bi = bi + 1;
                    j = j - 1; %repeat this closing branch to get the rest
                case 1  % parenthesis with no meaning (holds one element)
                    qp = qp - 1;
                    queue(qp) = queue(qp+1);
                    pi = pi + 1;
                case 0  % an empty parenthesis pair
                    error('Bioinfo:ParenthesisPairWithNoData', ...
                          ['Found parenthesis pair with no data,\n', ...
                           'input string may not be in Newick style or',...
                           'is not a valid filename.'])
            end % switch numElemInPar
            
        case '(' % an open parenthesis marker (0) pushed into the queue
            qp = qp + 1;
            queue(qp) = 0;
            
        otherwise % a new leaf pushed into the queue
            qp = qp + 1;
            queue(qp) = li;
            li = li + 1;
    end % switch strFeatures
    j = j + 1;
end % while j ...

catch 
   le = lasterror;
   if strcmp(le.identifier,'Bioinfo:ParenthesisPairWithNoData')
       rethrow(le)
   else
       error('Bioinfo:IncorrectString',...
             ['An error occurred while trying to interpret the data,\n'...
              'input string may not be in Newick style or is not a '...
              'valid filename.'])
   end
end    

% make sure all dists are greater than 0
dist = max(0,dist);

if sum(dist) == 0  % there was no distance information so force to an unitary ultrametric tree
    tr = phytree(tree,names);
elseif sum(dist(1:numLeaves)) == 0 % no dist infor for leaves, so force an ultrametric tree    
    tr = phytree(tree,dist(numLeaves+1:end),names);
else % put all info into output object
    tr = phytree(tree,dist,names);
end
