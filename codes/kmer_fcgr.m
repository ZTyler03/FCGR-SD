function [fcgr]=kmer_fcgr(seq,n,ka)
% seq: DNA sequence used to calculate FCGR
% n: length of kmer fragments corresponding to FCGR k
% ka: Alphabetical parameter of CGR, with values from 1 to 6, 
% corresponding to ACGT,ACTG,AGCT,AGTC,ATCG,ATGC, respectively.

    tmp = chaos(seq,ka);
    CGR_coordinates = tmp(:,n:end);

    N = 2 ^ n;
    fcgr = zeros(N,N);
    t = size(CGR_coordinates,2);
    for j = 1 : t
        x = ceil(N * CGR_coordinates(1,j));
        y = ceil(N * CGR_coordinates(2,j));
        fcgr(x,y) = fcgr(x,y) + 1;
    end
end
