
function generateM(nbInd::Int64=8, nbLocus::Int64=5, nbAllele::Int64=2)
    M = zeros(nbInd,nbLocus, nbAllele)
    for i in 1:nbInd
        for j in 1:nbLocus
            for k in 1:nbAllele
                M[i,j,k] = rand([0, 0.5, 1])
            end
        end
    end
    return M
end