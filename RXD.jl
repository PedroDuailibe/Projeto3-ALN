using LinearAlgebra

function R_Xd_manual(X, d, b)
    m = length(X)
    V = zeros(m, d+1)
    for i in 1:m
        for j in 0:d
            V[i, j+1] = X[i]^j
        end
    end
    c = V \ b   # resolve min ||Vc - b||_2
    return c
    # return x -> sum(c[j+1] * x^j for j in 0:d)
end

function teste_linearidade()

    println("#=== Variando o grau dos polinômios ===#")
    println("# teste R_Xd(a + b) = R_Xd(a) + R_Xd(b):")
    for d in 5:10
        println("# Teste com polinômios de grau $d")
        x = rand(1:10, (d + 2))
        a = rand(1:10, (d + 2))
        b = rand(1:10, (d + 2))
        
        r1 = R_Xd_manual(x, d, a)
        r2 = R_Xd_manual(x, d, b)
        r3 = R_Xd_manual(x, d, (a + b))

        if (norm((r3 - r2 - r1)) < 1e-10)
            println("True")
        else
            println("False")
        end

    end

    sleep(3)

    println("\n\n#=== teste R_Xd(ka) = k R_Xd(a) ===#")
    for d in 5:10
        println("# Teste com polinômios de grau $d")
        x = rand(1:10, (d + 2))
        a = rand(1:10, (d + 2))
        c = rand(1:10)
        
        r1 = R_Xd_manual(x, d, a)
        r2 = R_Xd_manual(x, d, (c*a))
        
        if (norm((r2 - c * r1)) < 1e-10)
            println("True")
        else
            println("False")
        end

    end

    sleep(3)

    println("\n#=== Variando o total de pontos ===#\n")
    for pts in 10:15
        println("# Teste com $pts pontos")
        x = rand(1:10, pts)
        a = rand(1:10, pts)
        b = rand(1:10, pts)
        a1 = rand(1:10)
        b1 = rand(1:10)

        r1 = R_Xd_manual(x, 8, a)
        r2 = R_Xd_manual(x, 8, b)
        r3 = R_Xd_manual(x, 8, (a1 * a + b1 * b))

        if(norm((r3 - a1 * r1 - b1 * r2)) < 1e-10)
            println("True")
        else
            println("False")
        end

    end
end

teste_linearidade()