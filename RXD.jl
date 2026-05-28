using LinearAlgebra
using Plots

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

function cond_Lx(x, m, d) 
    V = zeros(m, d + 1)
    for i in 1:m
        for j in 0:d
            V[i, j+1] = x[i]^j
        end
    end

    return cond(V)
end

function teste_cond_d()
    grau_polinomio = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    x = [1, 2, 3, 4, 5, 6, 7, 8, 9]
    d = zeros(10)

    for i in 1:10
        d[i] = log(cond_Lx(x, 6, i))
    end

    plot(grau_polinomio, d, marker=:circle)

end

function teste_cond_m()
    numero_pontos = [6, 7, 8, 9, 10]
    d = zeros(5)

    for i in 6:10

        x = zeros(i)
        for j in 1:i
            x[j] = j
        end

        d[i - 5] = log(cond(x, i, 5))
    end

    plot(numero_pontos, d, marker=:circle)
end
teste_cond()

#=
# Letra d)

Lema: O condicionamento da matriz de R_{Xd} é igual ao condicionamento da matriz de L_{X}

L_{X} na base canônica de Pd é a matriz V de vandermonde, onde V_{ij} = x_{i}^{j-1}. Também,
R_{Xd} é a pseudo-inversa da matriz de L_{X}, ou seja, (V^{T}V)^{-1}V^{T} = V†. Seja s_{1}, ..., s_{n}
os valores singulares de V, com s_{1} ≥ s_{2} ≥ ... ≥ s_{n}. Então existem listas ortonormais e_{1}, ..., e_{n} e f_{1}, ..., f_{n}
tais que:
Vu = s_{1}\left\langle u, e_{1} \right\rangle f_{1} + ... + s_{n} \left\langle u, e_{n} \right\rangle f_{n}
e
V†u = \frac{\left\langle u, f_{1} \right\rangle}{s_{1}} e_{1} + ... + \frac{\left\langle u, f_{n} \right\rangle}{s_{n}} e_{n}

Sabemos que o condicionamento de uma matriz é o maior valor singular dividido pelo menor. Então
k(V) = \frac{s_{1}}{s_{n}}
Enquanto isso, o maior valor singular de V† é \frac{1}{s_{n}}, enquanto o menor é \frac{1}{s_{1}}. Logo,
k(V†) = \frac{\frac{1}{s_{n}}} {\frac{1}{s_{1}}} = \frac{s_{1}}{s_{n}} = k(V)

Logo, os condicionamentos são iguais
Desse modo, basta analisar o comportamento do condicionamento da matriz de L_{X} quando d e m variam, 
pois este será igual ao da matriz de R_{Xd}.
=#