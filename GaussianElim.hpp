//  GaussianElim.hpp
//  A function that solves a linear system by Gaussian Elimination
//
//  Project 2
//
//  Created by Yaojia Huang on 2018/11/14.
//  Initially run in XCode(Version 10.1 (10B61))

#ifndef GaussianElim_hpp
#define GaussianElim_hpp
#include <string>


////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The function solves Ax = b, where A is a matrix and b is a vector                                      //
// Input: a N*N matrix, a vector of length N                                                              //
// Output: the vector solution                                                                            //
////////////////////////////////////////////////////////////////////////////////////////////////////////////
std::vector<double> GaussianElim(std::vector<std::vector<double>> mat, std::vector<double> b)
{
    // Dimension of the matrix
    long nrow = mat.size();

    // Normalise the first row
    // Make A[0][0] = 1
    b[0] /= mat[0][0];
    std::for_each(mat[0].begin(), mat[0].end(), [mat](double &el){el /= mat[0][0]; });
    
    // Loop over rows
    for (int i = 0; i < nrow; i++)
    {
        // Loop over rows preceding the current row
        for (int j = 0; j < i; j++)
        {
            // Subtract a multiple of preceding rows
            // to make 0s at certain position of the current row
            b[i] -= (mat[i][j] * b[j]);
            std::vector<double> replace(nrow);
            for (int k = 0; k < nrow; k++)
                replace[k] = mat[i][k] - mat[i][j] * mat[j][k];
                mat[i] = replace;
        }
        
        // Normalise the diagonal element
        b[i] /= mat[i][i];
        std::vector<double> replace(nrow);
        for (int k = 0; k < nrow; k++)
            replace[k] = mat[i][k] / mat[i][i];
            mat[i] = replace;
    }
    
    // Solve the remaining equation by substitution
    for (int i = 1; i < nrow+1; i++)
    {
        for (int j = 1; j < i; j++)
            b[nrow-i] -= mat[nrow-i][nrow-j] * b[nrow-j];
    }
    
    return b;
}
#endif
