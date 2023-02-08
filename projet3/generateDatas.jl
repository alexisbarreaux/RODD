
function generateM(nbInd::Int64=8, nbLocus::Int64=5, nbAllele::Int64=2)
    M = zeros(nbInd,nbLocus, nbAllele)
    for p in 1:nbInd
        for i in 1:nbLocus
            allelesLeft = nbAllele
            for j in 1:nbAllele
                if allelesLeft == 0
                    break
                else
                    if j == nbAllele
                        allelesOfTypeJ = allelesLeft
                    else
                        allelesOfTypeJ = rand(1:allelesLeft)
                    end
                    M[p,i,j] = allelesOfTypeJ / nbAllele
                    allelesLeft -= allelesOfTypeJ
                end
            end
        end
    end
    return M
end