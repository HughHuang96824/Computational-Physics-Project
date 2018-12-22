//  function.hpp
//
//  A function class that uses various metheds(e.g. Newton-Coates
//  and Monte Carlo methods) to do integrations
//
//  Project 2
//
//  Created by Yaojia Huang on 2018/11/14.
//  Initially run in XCode(Version 10.1 (10B61))

#ifndef function_hpp
#define function_hpp
#include <functional>
#include <iostream>
#include <random>


class Function {
private:
    // Default constructor
    // Not allowed to be accessed
    Function();
    
    // Attributes:
    // The function
    std::function<double(double)> f;
    // The starting and ending points, and the range
    double start, end, range;
    // The maximum of the function within the range
    double fmax;
    // Minimum samples that the program has to obtain to give the estimation
    int min_sample;
    
    // Structure that serves for Monte Carlo with adaptive importance sampling
    struct division;
    
    // Get the maximum of the function within the range
    double GetMax(const std::function<double(double)> &p) const;
    
    // Return kernel functions of randomly sampled points.
    // N is the number of the points want to be sampled.
    std::function<double(double)> Get_Kernel(const int &N) const;
    
    // Fit the inputted kernel function with a polynomail of a specified power
    std::vector<double> Fitting(const std::function<double(double)> &Kernel, const int &power) const;
    
public:
    
    // Constructors
    Function(const std::function<double(double)> &func);
    Function(const Function &func);
    
    // Destructor
    ~Function();
    
    
    // Making the function callable
    double operator()(const double &x) const;
    
    // Set integration interval range
    void Range(const double &a, const double &b);
    
    // Set the minimum number of samples
    void MinSample(const int &n);
    
    // Set the function
    void Func(const std::function<double(double)> &func);
    
    // Use trapezoidal rule to estimate the area
    // Take epsilon as input and return the area
    double Trapezoidal(const double &epsilon) const;
    
    // Use Simpson rule to estimate the area
    // Take epsilon as input and return the area
    double Simpson(const double &epsilon) const;
    
    // Use naive Monte Carlo to estimate the area
    // Take epsilon as input and return the area
    double MC(const double &err) const;
    
    // Use Monte Carlo with importance sampling to estimate the area
    // Take a sampling PDF, a inverted CDF of the PDF, and a specified err
    // Return the area
    double MC(const std::function<double(double)> &q, const std::function<double(double)> &invert_CDF, const double &err) const;
    
    // Use Monte Carlo with importance sampling to estimate the area
    // Take a sampling PDF and a specified err
    // Return the area
    double MC(const std::function<double(double)> &q, const double &err) const;
    
    // Use Monte Carlo with non-parametric importance sampling to estimate the area
    // Input: the numer of random samples for kernel density estimation, power of the fitting polynomial,
    //        a specified error, the number of times that the user want to update the PDF,
    //        the path of the txt file that stores the information of the PDF
    // Output: the result of integration
    double NPIS(const int &N, const int &degree, const double &err, int Nupdate = 0, std::string path = " ") const;
    
    // Use Monte Carlo with adaptive importance sampling to estimate the area
    // The method uses VEGAS algorithm that updates the stratified PDF based on last estimated numerics.
    // Input: Initial number of bins for stratified sampling, a specified error
    // Output: the result of integration
    double AIS(int N, const double &err) const;
};

#endif
