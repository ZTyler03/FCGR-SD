%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
% This program uses the FCGR-SD method to calculate the 
% Neighbor-Joining phylogenetic tree of DNA sequences.
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
K = 7;
% Controls the fineness of the fcgr chunking, which essentially 
% corresponds to the k-value of the K-mer in the sequence
    
fastaname = '41mito.fasta';
% Fasta file path string, the fasta file should contain all the sequences
% [Please change to your own file path]
    
seqtype = "41mito";
outhead = strcat('AGCT_',seqtype,'_K');
% Filename tagging of resulting phylonegetic tree files
    
File = fastaread(fastaname);
% Read fasta file
    
n = size(File,1);
% Get the number of sequence entries
    
N = 2^K;
% Obtain the fineness of the fcgr delineation
   
fcgr = cell(n,1);
% Store the FCGR matrices (in frequency form), 
% each fcgr{i} is an N*N matrix, here N=2^k
    
A_trunk = cell(n,1);
% Storing truncated feature matrices
    
inner_product = cell(1);
% Storing the product matrix of truncated feature matrices
      
theta = cell(1);
% Store the collection of protagonists
    
threshold = cos(10^(-5));
% The singularity threshold that determines whether the principal angle 
% computation applies acos or asin
    
distance_grassman = zeros(n,n);
% Storing the generalized Grassman distance matrix
    
fprintf('Calculation begins!\n');
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 1: Calculating the FCGR Matrix
    
    timer = 0;
    tic;
        
    for i = 1 : n
        fcgr{i} = kmer_fcgr(File(i).Sequence,K,3);
    end
    timer = timer + toc;
        
    fprintf('Step %d finished in %f seconds!\n',1,timer);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 2: The first stage of singular value decomposition to find the 
% set of truncated feature matrices {A_i}
    
    tic;
    for j = 1 : n
        [U,diag1,~] = svd(fcgr{j},"vector");
        tmp = sum(diag1);
        for r = 1 : N
            if sum(diag11(1:r)) >= P * tmp
                A_trunk{j} = U(:,1:r);
                break;
            end
        end
        % Intercepting the left eigenmatrix based on the cumulative 
        % proportion of singular values
    end
    timer = timer + toc;
        
    fprintf('Step 1 to %d finished in %f seconds!\n',2,timer);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 3: Stage 2 singular value decomposition to compute the 
% principal angles and distance matrices
        
    tic;
    for s = 1 : n - 1
        for t = s + 1 : n

            inner_product = A_trunk{s}' * A_trunk{t};
            % First compute the product of the two truncated matrices

            [~,diag2,VV] = svd(inner_product,"vector");
            % Singular value decomposition of the product matrix

            teser = size(diag2,1);
            theta = zeros(teser,1);
            % Since the length and value of theta changes with each loop,
            % initialize theta before calculating it.

            if diag2(1) < threshold
                theta = acos(diag2);
                % If the largest singular value is below the threshold 
                % (the smallest principal angle is greater than the 
                % threshold), then directly use the acos method

            else
                kapa = find(diag2 < threshold,1);
                % Otherwise, we first determine the location of the cutoff
                % where the singular value is below the threshold (the 
                % principal angle is greater than the threshold)

                if isempty(kapa)
                    kapa = teser + 1;
                    % Prevent subscript overflow when all singular values
                    % are higher than the threshold (all principal angles 
                    % are less than the threshold)
                else
                    theta(kapa:teser) = acos(diag2(kapa:teser));
                    % For the portion of the singular values below the 
                    % threshold, the acos method is taken
                end

                Vcos = A_trunk{t} * VV;
                Bmat = Vcos(:,1:kapa - 1) - A_trunk{s} * (A_trunk{s}' * Vcos(:,1:kapa - 1));
                % The singular values of Bmat are the sine values 
                % corresponding to those small principal angles

                [~,SS1,~] = svd(Bmat,"vector");
                diag2tmp = SS1;
                theta(1:kapa - 1) = asin(diag2tmp);
                % For the portion of the singular values greater than 
                % or equal to the threshold (small angles), the asin 
                % method is used
            end

            distance_grassman(s,t) = sqrt(theta(:)' * theta(:));
            distance_grassman(t,s) = distance_grassman(s,t);
            % Calculate the generalized Grassmannian distance
        end 
    end

    timer = timer + toc;
        
    fprintf('Step 1 to %d finished in %f seconds!\n',3,timer);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Export the phylogenetic tree

tree = seqneighjoin(distance_grassman,'equivar',File);
% where the parameter 'equivar' indicates that the Neighbor-Joining 
% method was used to build the tree, if it is changed to 'firstorder' 
% then the reslting tree is built by BIONJ method.

phytreewrite(strcat(outhead,num2str(K),'_testcutnum.nwk'),tree);
% Exporting phylogenetic tree files
    


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%