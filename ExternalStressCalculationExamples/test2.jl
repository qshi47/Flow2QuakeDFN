function main()
    P = 28
    nu = 0.15
    Phi = 70
    Sigma_H = P*nu/(1-nu)
    t = sqrt( (Sigma_H*sind(Phi))^2+(P*cosd(70))^2 )
    Sigma_N = Sigma_H*sind(Phi)^2 + P*cosd(Phi)^2
    println("Sigma_N is ", Sigma_N)
    Sigma_S = sqrt(t^2 - Sigma_N^2)
    println("Sigma_S is ", Sigma_S)
end

main()