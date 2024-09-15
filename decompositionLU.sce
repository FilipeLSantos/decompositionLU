function [Inv]=LU(A)
    function [A, Pivot]=decomposicaoLU(A)
        n = size(A, 1);
    
        for i = 1 : n
            Pivot(i) = i;
        end
        
        for j = 1 : n -1
            p = j;
            Amax = abs(A(j,j));
            
            for k = j + 1 : n
                if abs(A(k,j)) > Amax then
                    Amax = abs(A(k, j));
                    p = k;
                end
            end
            
            if p ~= j then
                for k = 1 : n
                    t = A(j, k);
                    A(j, k) = A(p, k);
                    A(p, k) = t;
                end
                
                m = Pivot(j);
                Pivot(j) = Pivot(p);
                Pivot(p) = m;
            end
                
            if abs(A(j, j)) ~= 0 then
                r =  1 / A(j, j);
                
                for i = j + 1 : n
                    mult = A(i, j) * r;
                    A(i, j) = mult;
                    
                    for k = j + 1 : n
                        A(i, k) = A(i, k) - mult * A(j, k);
                    end
                end
            end
        end    
    endfunction
    function [y]=subsSucessivas(L, b, Pivot)
        n = size(L, 1);
        
        k = Pivot(1);
        y(1) = b(k);
        
        for i = 2 : n
            Soma = 0;
            
            for j = 1 : i - 1
                Soma = Soma + L(i, j) * y(j);
            end
            
            k = Pivot(i);
            y(i) = b(k) - Soma;
        end
        
    endfunction

    function [x]=subsRetroativas(U, d)
        n = size(U, 1);
        
        x(n) = d(n)/U(n, n);
        
        for i = n - 1 : -1 : 1
            Soma = 0;
            
            for j = i + 1 : n
                Soma = Soma + U(i, j) * x(j);
            end
            
            x(i) = (d(i) - Soma)/U(i, i);
        end
    endfunction
    
    n = size(A, 1);
    I = eye(n,n);
    [LU, Pivot] = decomposicaoLU(A);
    
    for i = 1 : size(A, 1)
        y = subsSucessivas(LU, I(:,i), Pivot);
        Inv(:,i) = subsRetroativas(LU, y);       
    end
endfunction
