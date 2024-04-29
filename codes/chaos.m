function x=chaos(seq,kap)
% seq: DNA sequence used to convert to CGR image coordinates
% kap: Alphabetical parameter of CGR, with values from 1 to 6, 
% corresponding to ACGT,ACTG,AGCT,AGTC,ATCG,ATGC respectively.
    
    seq=upper(seq);
    % Uniform capitalization of sequence letters

    if kap == 1
        p_A=[0;0];
        p_C=[0;1];
        p_G=[1;1];
        p_T=[1;0];
    elseif kap == 2
        p_A=[0;0];
        p_C=[0;1];
        p_T=[1;1];
        p_G=[1;0];
    elseif kap == 3
        p_A=[0;0];
        p_G=[0;1];
        p_C=[1;1];
        p_T=[1;0];
    elseif kap == 4
        p_A=[0;0];
        p_G=[0;1];
        p_T=[1;1];
        p_C=[1;0];
    elseif kap == 5
        p_A=[0;0];
        p_T=[0;1];
        p_C=[1;1];
        p_G=[1;0];
    elseif kap == 6
        p_A=[0;0];
        p_T=[0;1];
        p_G=[1;1];
        p_C=[1;0];
    end
    % Determine which bases are represented in each of the four corners 
    % of the unit square based on the value of kap

    m=size(seq,2);
    x=zeros(2,m+1);
    x(1,1)=1/2;
    x(2,1)=1/2;
    % The initial point of the CGR is (1/2,1/2)

    for j=1:m

        if j >= 50
        % After the 50th nucleotide, it is possible to begin to be 
        % affected by the precision of the floating point number.

            if x(1,j) <= 2^(-500)
                x(1,j) = x(1,j) + 2^(-31);

            elseif abs(x(1,j)-0.5) <= 2^(-50)
                if x(1,j-1) < 0.5
                    x(1,j) = x(1,j) + 2^(-31);
                else
                    x(1,j) = x(1,j) - 2^(-31);
                end
            elseif 1-x(1,j) <= 2^(-50)
                x(1,j) = x(1,j) - 2^(-31);
            end
            if x(2,j) <= 2^(-500)
                x(2,j) = x(2,j) + 2^(-31);
            elseif abs(x(2,j)-0.5) <= 2^(-50)
                if x(2,j-1) < 0.5
                    x(2,j) = x(2,j) + 2^(-31);
                else
                    x(2,j) = x(2,j) - 2^(-31);
                end
            elseif 1-x(2,j) <= 2^(-50)
                x(2,j) = x(2,j) - 2^(-31);
            end
        end

        switch seq(j)
            case 'A'
                x(:,j+1)=(x(:,j)+p_A)/2;
            case 'C'
                x(:,j+1)=(x(:,j)+p_C)/2;
            case 'G'
                x(:,j+1)=(x(:,j)+p_G)/2;
            case 'T'
                x(:,j+1)=(x(:,j)+p_T)/2;
        end
    end
    
    x=x(:,2:end);
    % Remove the initial point
end