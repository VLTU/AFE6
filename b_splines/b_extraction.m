function C = b_extraction(U,p)
    nel = length(unique(U))-1;
    m = length(U);
    a = p+1;
    b = a+1;
    nb = 1;
    C{1}=eye(p+1);
    while b<m
        C{nb+1}= eye(p+1); 
        %/ Initialize the next extraction operator.
        i=b;
        %Count multiplicity of the knot at location b.
        while b < m && U(b+1) == U(b)
            b=b+1;
        end
        mult = b-i+1;
        if mult < p
            %Use (10) to compute the alphas.
            numer = U(b)-U(a);
            for j = p:-1:mult+1
                alphas(j-mult) = numer / (U(a+j)-U(a));
            end
            r = p-mult;
            %Update the matrix coefficients for r new knots
            for j=1:r
                save = r-j+1;
                s = mult+j;
                for k=p+1:-1:s+1
                    alpha = alphas(k-s);
                    %The following line corresponds to (9).
                    C{nb}(:,k) = alpha*C{nb}(:,k) + (1.0-alpha)*C{nb}(:,k-1);
                end
                if b<m 
                    %Update overlapping coefficients of the next operator.
                    C{nb+1}(save:j+save,save) = C{nb}(p-j+1:p+1,p+1);
                end
            end
        end
        nb = nb + 1;
        %// Finished with the current operator.
        if b<m 
            %Update indices for the next operator.
            a=b;b = b+1;
        end
    end
    C = C(1:nel);
end 