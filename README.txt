This program was coded using C++ in XCode 10.1
It is a class for 1D function. 
The main functionality of this class are evaluating integrals using numerical methods.

CAUTION: Since this program is initially designed on Mac system, some parts of the code, such as file path, may (or may not) need adjustment.




Files:
function.hpp
function.cpp
GaussianElim.hpp
main.cpp
coefficient.txt (may or may not exist, does not matter)





Recommendation for reading comments:

When reading the comments in function.cpp, it would be easier to read the comments of the
public member functions before reading those of private member functions because private
functions are used to serve the public functions. So public functions would be a proper
Start.





Instruction for running the program:

The functions open to public are already coded in main.cpp. The user only needs to uncomment the code to run. It could take a longer time to run if all lines are uncommented. All functions will show the number of samples/divisions used and the epsilon of the result.

The order I wrote the arguments of the function below is also the order should be followed in the program.


Example for creating a function object:
Function prb([](double x) {return pow(M_PI, -0.5) * exp(-x*x);});


The integration range must be set to do integration.
Such as prb.Range(0, 2);


The user can set a minimum sample size. The default value is 5000.
prb.MinSample(10000);


Trapezoidal() and Simpson() use the trapezoidal and Simpson's rule to calculate the integral. They both takes a relative accuracy as input.



MC() is a overloaded function.

If user enters only a relative accuracy, the integration is evaluated using naive Monte Carlo method -- throwing dots.

If the user enters a PDF and a relative accuracy, the integration is evaluated using importance sampling, but samples are picked using rejection method.

If the user provides a PDF, an inverted CDF of the PDF and a relative accuracy, the integration is evaluated using importance sampling. The samples are generated from the inverted CDF.



NPIS() is an abbreviation for non-parametric importance sampling. It uses kernel density estimation to get a rough PDF, and refines it by fitting a polynomial function.

The inputs for this function are: 
1. the number of kernel functions the user want to use to estimate the PDF (larger the better, but slower);
2. the power of the fitting polynomial function;
3. a relative accuracy;
4. the number of polynomial PDFs the user wants to generate (only one of them is adopted and stored). If this value is not entered, the program uses previous generated PDF.

This function would show the relative accuracy of the accepted PDF during execution. After the required number of PDF are generated, the program will ask the user whether he/she wants to produce more PDFs or calculate the integral straight away. Please follow the instruction given by the program. After the program terminates, there will be a text file that records the coefficients of the polynomial and the relative accuracy in the same folder as the program files.



AIS() uses VEGAS algorithm. It asks the user for the number of initial bins and a relative accuracy.

