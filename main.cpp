//  main.cpp
//  Test program for the Function class
//
//  Project 2
//
//  Created by Yaojia Huang on 2018/11/14.
//  Initially run in XCode(Version 10.1 (10B61))

#define _USE_MATH_DEFINES
#include <iostream>
#include <functional>
#include <cmath>
#include <vector>
#include <sstream>
#include <iomanip>
#include "function.hpp"
using namespace std;


int main() {

    // Function that needs to be integrated
    Function prb([](double x) {return pow(M_PI, -0.5) * exp(-x*x);});
    // Integration range
    prb.Range(0, 2);
    // Minimum samples required
    prb.MinSample(10000);


    // Use the trapezoidal rule
    double area2 = prb.Trapezoidal(pow(10, -6));
    std::cout << "Trapezoidal: " << std::setprecision(8) << area2 << "\n" << endl;


    // Use Simpson's rule
    double area1 = prb.Simpson(pow(10, -6));
    std::cout << "Simpson: " << area1 << "\n" << std::endl;


//    double naive = prb.MC(0.001);
//    cout << "Naive Monte Carlo: " << naive << "\n" << endl;
//
//
//    // Uniform sampling PDF function
//    Function pdf1([](double x) {return 1.0/2.0;});
//    // Inverted CDF of the PDF
//    Function invert_CDF([](double x) {return 2.0 * x;});
//
//    // Monte Carlo integration with uniform sampling
//    double result = prb.MC(pdf1, invert_CDF, 0.0001);
//    cout << "Uniform Sampling: " << result << "\n" <<  endl;


//    // Weighted sampling PDF function
//    Function pdf2([](double x) {return -0.48 * x + 0.98;});
//
//    // Monte Carlo integration with weighted sampling
//    double result2 = prb.MC(pdf2, 0.0001);
//    cout << "Linear Weighted Sampling: " <<  result2 << "\n" <<  endl;


//    // Set a minimum sample size
//    prb.MinSample(100000);
//
//    // Monte Carlo integration with non-parametric importance sampling
//    double NPIS_RESULT = prb.NPIS(500000, 6, 0.000001, 5);
//    cout << "Non-parametric Sampling: " << std::setprecision(7) << NPIS_RESULT<< endl;
    
    
//    // Monte Carlo integration with adaptive importance sampling
//    double x = prb.AIS(100, 0.0000001);
//    cout << std::setprecision(8) << x << endl;
}
