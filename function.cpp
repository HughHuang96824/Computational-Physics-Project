//  function.cpp
//  A function class that uses various metheds(e.g. Newton-Coates
//  and Monte Carlo methods) to do integrations
//
//  Project 2
//
//  Created by Yaojia Huang on 2018/11/14.
//  Initially run in XCode(Version 10.1 (10B61))

#define _USE_MATH_DEFINES
#include <set>
#include <cmath>
#include <cfloat>
#include <fstream>
#include <sstream>
#include <iterator>
#include <iomanip>
#include <map>
#include "function.hpp"
#include "GaussianElim.hpp"


////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////Private Function////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////


// This structure is designed for adaptive importance sampling(AIS) method.  Not for private usage
// See member function AIS() for detail
struct Function::division
{
    double average = 0;
    int NSubsub = 0;
};


////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The function estimates the maximum value of a given function within the integration region.            //
// Input: a function                                                                                      //
// Output: maximum value of the function within the integration region                                    //
////////////////////////////////////////////////////////////////////////////////////////////////////////////
double Function::GetMax(const std::function<double(double)> &p) const
{
    std::random_device rd;
    std::mt19937 mt(rd());
    std::uniform_real_distribution<double> x_dist(start, end);
    
    // Step constraint
    double step = range;
    std::uniform_real_distribution<double> move(-step, step);
    
    // Initial temperature and "Boltzann constant"
    const double k = 1e-6 ;
    double T = 100 * range / k;
    
    // Boltzmann Distribution
    const std::function<double(double, double, double)> boltzmann =
    [k](double T, double E1, double E2) {return exp((E2 - E1)/k/T);};
    
    double x0 = (start + end) / 2; // Initial loacation, in the middle
    double x_new;                  // Next location
    double E1 = p(x0);             // Initial function value
    double E2;                     // Next function value
    
    // Variable that precludes deadloop
    int no_deadloop = 0;

    // Iterating until temperature reaches near 0 or samples are trapped.
    while (T > 200 && no_deadloop < 1000)
    {   // Get next location if it is within the boundaries
        do x_new = x0 + move(mt);
        while(x_new < start || x_new > end);
        
        // Get the value of the function at the next location
        E2 = p(x_new);
        
        // If the value becomes larger or is accepted by chance
        std::uniform_real_distribution<double> accept(0, 1);
        if (E2 > E1 || boltzmann(T, E1, E2) > accept(mt))
        {   // Update the x-coordinate
            x0 = x_new;
            // Update the maximum value of the function
            E1 = E2;
    
            T -= 200;        // Annealing
            no_deadloop = 0; // Reset no_deadloop
        }
        else
            no_deadloop++;   // Increment no_deadloop if rejected
    }
    
    // Return the maximum value of function within the boundaries
    // Coefficient 1.1 is used to prevent possible error
    return 1.1 * p(x0);
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The function estimates the PDF that has the same shape as the original function by sampling points     //
// Input: number of sampling points/kernel functions                                                      //
// Output: Estimated PDF                                                                                  //
////////////////////////////////////////////////////////////////////////////////////////////////////////////
std::function<double(double)> Function::Get_Kernel(const int &N) const
{
    std::random_device rd;
    std::mt19937 mt(rd());
    std::uniform_real_distribution<double> x_dist(start, end);
    std::uniform_real_distribution<double> y_dist(0, fmax);
    int count = 0;
    
    // x and y samples
    double _x, _y;

    // A single kernel function
    std::function<double(double)> K;
    
    // Width of each gaussian kernel function, determined by the rule of thumb
    double b = 1.06 * pow(N, -0.2);
    
    // Vector of N kernel functions
    std::vector<std::function<double(double)>> Kernel_list; Kernel_list.reserve(N);
    while (count < N)
    {   // Get randomly sampled x and y
        _x = x_dist(mt);
        _y = y_dist(mt);
        
        // Accept if following conditions are satisfied
        if (f(_x) >= _y && _x != start && _x != end)
        {   // Set single kernel function
            K = [_x, b, this](double x)
            {return 1 / sqrt(2*M_PI) *
                    exp(-pow(log( (x-start) * (end-_x) / (_x-start) / (end-x) ) / b, 2) / 2);};
            
            // Add individual kernel function to the kernel list
            Kernel_list.push_back(K);
            count++;
        }
    }
    
    // This is the nomalising coefficient of kernel estimation
    // Add it to the end of the kernel list
    K = [N, b, this](const double& x)
    {return range/(double)N/b/(x-start)/(end - x);};
    Kernel_list.push_back(K);
    
    // Get the PDF by adding individual kernel functions and normalising
    std::function<double(double)> Kernel = [Kernel_list](double x)
    {   double sum = 0;
        for (int i = 0; i < Kernel_list.size() - 1; i++)
            sum += Kernel_list[i](x);
        sum *= (*std::prev(Kernel_list.end()))(x); // Normalising
        return sum;
    };
    
    // Return kernel estimation
    return Kernel;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The function calcultes the polynomial function that fits the kernel estimation function best.          //
// Input: the Kernel function, the power of the fitting polynomial                                        //
// Output: a vector of coefficients of the polynomial terms, ranked from the highest power to the lowest  //
////////////////////////////////////////////////////////////////////////////////////////////////////////////
std::vector<double> Function::Fitting(const std::function<double(double)> &Kernel, const int &power) const
{
    std::random_device rd;
    std::mt19937 mt(rd());
    std::uniform_real_distribution<double> dist(-0.005, 0.005);
    
    // Select 25 * range samples from the kernel density estimation
    int Npoints = (int)(25 * range);
    double step = (end - start) / (Npoints - 1);
    
    // Vector for values of the kernel density estimation at uniformly sampled x coordinates
    std::vector<double> values; values.reserve(Npoints);
    
    // Iterate over sampled points
    // Excluding the value at x = start or end
    // because the kernel function has unpredictable behaviours there
    for (double i = start + step; i < end - step; i += step)
    {
        // Store the value of kernel density estimation at x_i
        values.push_back(Kernel(i));
    }
    
    // To get the coefficients of the best fitting function analytical,
    // laborious calculations are involved. The following number lists,
    // number1, number2 and number3, contain the elements that are required
    // to compute this system.
    
    // Number1
    std::vector<double> numbers1; numbers1.reserve(2*power);
    for (int i = 2 * power; i >= 1; i--)
    {
        double tmp1 = 0;
        std::vector<double>::const_iterator iter = values.begin();
        for (double j = start + step; j < end - step; j += step)
        {
            tmp1 += pow(j, i) / pow((*iter), 2);
            iter++;
        }
        numbers1.push_back(tmp1);
    }
    
    // Number2
    std::vector<double> numbers2; numbers2.reserve(power);
    for (int i = power; i >= 1; i--)
    {
        double tmp2 = 0;
        std::vector<double>::const_iterator iter = values.begin();
        for (double j = start + step; j < end - step; j += step)
        {
            tmp2 += pow(j, i) / (*iter);
            iter++;
        }
        numbers2.push_back(tmp2);
    }
    
    // Number3
    std::vector<double> numbers3; numbers3.reserve(power);
    for (int i = power + 1; i >= 2; i--)
        numbers3.push_back((pow(end, i) - pow(start, i))/range/i);
    
    // Construct the matrix using above number lists
    std::vector<std::vector<double>> matrix(power);
    for (int i = 0; i < power; i++)
    {
        std::vector<double> row(power);
        for (int j = 0; j < power; j++)
            row[j] = numbers1[i+j] - numbers3[j] * numbers1[i+power];
        matrix[i] = row;
    }
    
    // Construct the b vector (Ax = b)
    std::vector<double> b(power);
    for (int i = 0; i < power; i++)
    {
        b[i] = numbers2[i] - numbers1[i+power] / range;
    }
    
    // Get the coefficients of the fitting function excluding the constant term
    // Coefficient of highest power term is in the most front of the vector
    std::vector<double> coefficient = GaussianElim(matrix, b);
    
    // Get the constant term, which ensures that the area of the fitting function
    // within the integration region is 1 so that it is a PDF.
    double const_term = 0;
    for (int i = power + 1; i > 1; i--)
        const_term += (pow(end, i) - pow(start, i)) * coefficient[power+1-i] / i;
    const_term = (1 - const_term) / range;
    
    // Add the constant term to the coefficient vector
    coefficient.push_back(const_term);
    
    // Return the coefficient vector
    return coefficient;
}



////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////Public Function////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////


/////////////////////////////////////////Constructors & Destructor//////////////////////////////////////////

// Constructor that sets the function
Function::Function(const std::function<double(double)> &func):
f(func), start(0), end(0), range(-1), min_sample(5000) {};

// Copy Constructor
Function::Function(const Function &func):
f(func.f), start(func.start), end(func.end), range(func.range), fmax(func.fmax), min_sample(5000) {};

// Destructor
Function::~Function() {}

//////////////////////////////////////////////////Getters///////////////////////////////////////////////////

// Calculate the value of the function given the value of x
double Function::operator()(const double &x) const
{
    return f(x);
}

/////////////////////////////////////////////////Modifiers//////////////////////////////////////////////////

// Set integration interval
void Function::Range(const double &a, const double &b)
{   // The smaller number is start and the other end
    start = a <= b? a:b;
    end   = a <= b? b:a;
    range = end - start;
    // Get the maximum
    fmax = GetMax(f);
}


// Set the minimum number of samples
void Function::MinSample(const int &n)
{
    min_sample = n;
}


// Set the function
void Function::Func(const std::function<double(double)> &func)
{
    f = func;
    
    // If the range has be set(default value is -1), get the maximum
    if (range >= 0)
        fmax = GetMax(f);
}


//////////////////////////////////////////Newton-Coates integrations////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Use trapezoidal rule to estimate the area                                                              //
// Input: a specifed error                                                                                //
// Output: the result of integration                                                                      //
////////////////////////////////////////////////////////////////////////////////////////////////////////////
double Function::Trapezoidal(const double &epsilon) const
{
    // Count number of interval divisions
    int counter = 0;
    
    // Initial step size
    double step = end - start;
    
    // Set that stores the sampling points in order
    std::set<double> x_points;
    x_points.insert(start);
    x_points.insert(end);
    
    // Calculte the initial area using only one interval
    double old_area = step * ((*this)(start) / 2 + (*this)(end) / 2);
    
    // Variable for updated area
    double new_area;
    
    // Variable that stores old area temporarily
    double tmp;
    
    new_area = old_area;
    do
    {
        
        // iterator over previous sampled points
        for (std::set<double>::iterator iter = std::next(x_points.begin()); iter != x_points.end(); iter++)
        {
            // Get a new sampling point in the middle of two old ones
            double new_x = ( *iter + *(std::prev(iter)) ) / 2;
            
            // Get the value of the function at the new point and increase the area
            new_area += step * (*this)(new_x);
            
            // Put newly sampled point in the set
            x_points.insert(new_x);
        }
        // Divide the step and the area by 2
        step /= 2;
        new_area /= 2;
        
        // Increase count
        counter++;
        
        tmp = old_area;
        old_area = new_area;
        
        // If the percentage change is larger than epsilon, repeat the above process.
    }   while (abs(new_area - tmp)/tmp > epsilon);

    std::cout << "Epsilon:   " << abs(new_area - tmp)/tmp << std::endl;
    std::cout << "Divisions: " << pow(2, counter) << std::endl;
    
    return new_area;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Use Simpson rule to estimate the area                                                                  //
// It is basically a little modification of the trapezoidal rule                                          //
// Input: a specifed error                                                                                //
// Output: the result of integration                                                                      //
////////////////////////////////////////////////////////////////////////////////////////////////////////////
double Function::Simpson(const double &epsilon) const
{
    // Count number of interval divisions
    int counter = 0;
    
    // Initial step size
    double step = end - start;
    
    // Set that stores the sampling points in order
    std::set<double> x_points;
    x_points.insert(start);
    x_points.insert(end);
    
    // Calculte the initial trapezoidal area using only one interval
    double old_T = step * ((*this)(start) / 2 + (*this)(end) / 2);
    
    // Variable for areas
    double new_T;          // Updated trapezoidal area
    double new_S = old_T;  // Updated Simpson area
    double old_S;          // Last Simpson area
    new_T = old_T;
    
    do
    {
        old_S = new_S;
        
        // iterator over previous sampled points
        for (std::set<double>::iterator iter = std::next(x_points.begin()); iter != x_points.end(); iter++)
        {
            // Get a new sampling point in the middle of two old ones
            double new_x = ( *iter + *(std::prev(iter)) ) / 2;
            
            // Get the value of the function at the new point and increase the area
            new_T += step * (*this)(new_x);
            
            // Put newly sampled point in the set
            x_points.insert(new_x);
        }
        // Divide the step and the area by 2
        step /= 2;
        new_T /= 2;
        
        // Increase count
        counter++;
        
        // The i-th Simpson rule can be obtained by weighing the (i+1)-th and
        // the i-th tratrapezoidal areas differently, and adding them together
        new_S = 4.0 * new_T / 3.0 - old_T / 3.0;
        old_T = new_T;
        
        // If the percentage change is larger than epsilon, repeat the above process.
    }   while (abs(new_S - old_S)/old_S > epsilon);

    std::cout << "Epsilon:   " << abs(new_S - old_S)/old_S << std::endl;
    std::cout << "Divisions: " << pow(2, counter) << std::endl;
    
    return new_S;
}


//////////////////////////////////////////Monte Carlo Integration///////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////////////////////////////////
// This function is the implementation of the naive Monte Carlo Method                                    //
// Input: a specified error                                                                               //
// Output: the result of integration                                                                      //
////////////////////////////////////////////////////////////////////////////////////////////////////////////
double Function::MC(const double &err) const
{
    // Get the maximum of the function within the range
    std::random_device rd;
    std::mt19937 mt(rd());
    std::uniform_real_distribution<double> dist_x(start, end);
    std::uniform_real_distribution<double> dist_y(0, fmax);
    
    // Number of accepted points and total points
    long accepted = 0;
    long count = 0;
    
    // Area of the bounding region
    double Outer_area = fmax * range;
    // Error
    double Epsilon;
    double error_smaller = 0;
    
    // Get samples
    do
    {   // Random sampled coordinates
        double x = dist_x(mt);
        double y = dist_y(mt);
        
        // If the coordinate is below f(x), accept the sample
        if (y <= f(x))
            accepted++;
        
        // Increase total count
        count++;
        
        // Accept probability
        double p = accepted / (double)count;
        // Get error
        Epsilon = sqrt((p * (1-p)) * count) / accepted;

        if (Epsilon < err && count > min_sample)
            error_smaller++;
        else
            error_smaller = 0;
        
      // Continue if the the Epsilon is not successively less than err for 100 times
      // Prevent random fluctuation
    } while (error_smaller < 100);
    
    std::cout << "Epsilon: " << Epsilon << "\n" << "NSample: " << count << std::endl;
    
    // Return the estimated area
    return Outer_area * (accepted / (double)count);
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////
// This function is for Monte Carlo integration with importance sampling                                  //
// Use this function if the CDF of the sampling distribution is invertible                                //
// Input: the sampling PDF, the inverted samling CDF, a specifed error                                    //
// Output: the result of integration                                                                      //
////////////////////////////////////////////////////////////////////////////////////////////////////////////
double Function::MC
(const std::function<double(double)> &q,
 const std::function<double(double)> &invert_CDF, const double &err) const
{
    std::random_device rd;
    std::mt19937 mt(rd());
    std::uniform_real_distribution<double> dist(0, 1);
    double area = 0;
    double squared_sum = 0; // Sum of the square of each sample value
    double Epsilon = 0;
    long count = 0;         // Number of samples
    int error_smaller = 0;  // Number of successive times that the Epsilon is less than err
    
    // Continue if the the Epsilon is not successively less than err for 100 times
    // Prevent random fluctuation
    while (error_smaller < 100)
    {
        // Sampling x coordinate within the intergration range
        double x = invert_CDF(dist(mt)) + start;
        
        // tmp is the value of each sample
        double tmp = f(x) / q(x);
        
        // Increase area and the sum of tmp^2
        area += tmp;
        squared_sum += pow(tmp, 2);
        // Increase count
        count++;
        // Calculate the Epsilon
        Epsilon = sqrt((squared_sum/count - pow(area/count, 2)) * count) / area;

        // Increment error_smaller or make it 0
        if (Epsilon < err && count > min_sample)
            error_smaller++;
        else
            error_smaller = 0;
    }

    // Take the average
    area /= count;
    std::cout << "Epsilon: " << Epsilon << "\n" << "NSample: " << count << std::endl;
    
    return area;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////
// This function is for Monte Carlo integration with importance sampling                                  //
// Use this function if the CDF of the sampling distribution is NOT invertible                            //
// Input: the sampling PDF, a specifed error                                                              //
// Output: the result of integration                                                                      //
////////////////////////////////////////////////////////////////////////////////////////////////////////////
double Function::MC(const std::function<double(double)> &q, const double &err) const
{
    double max = GetMax(q);
    std::random_device rd;
    std::mt19937 mt(rd());
    std::uniform_real_distribution<double> dist_x(start, end);
    std::uniform_real_distribution<double> dist_y(0, max);
    double area = 0;
    double squared_sum = 0; // Sum of the square of each sample value
    double Epsilon = 0;
    long count = 0;         // Number of samples
    int error_smaller = 0;  // Number of successive times that the Epsilon is less than err
    
    // Continue if the the Epsilon is not successively less than err for 100 times
    // Prevent random fluctuation
    while (error_smaller < 100)
    {
        double x;
        
        do x = dist_x(mt);
        while(dist_y(mt) > q(x));
        
        // tmp is the value of each sample
        double tmp = f(x) / q(x);
        
        // Increase area and the sum of tmp^2
        area += tmp;
        squared_sum += pow(tmp, 2);
        // Increase count
        count++;
        // Calculate the Epsilon
        Epsilon = sqrt((squared_sum/count - pow(area/count, 2)) * count) / area;
        
        // Increment error_smaller or make it 0
        if (Epsilon < err && count > min_sample)
            error_smaller++;
        else
            error_smaller = 0;
    }

    // Take the average
    area /= count;
    std::cout << "Epsilon: " << Epsilon << "\n" << "NSample: " << count << std::endl;
    
    return area;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////
// This function is for Monte Carlo integration with non-parametric importance sampling                   //
// The sampling PDF is obtained by fitting a polynomial to the kernel density estimation function.        //
// The most fitting PDF is used for sampling.                                                             //
// Input: the numer of random samples for kernel density estimation, power of the fitting polynomial,     //
//        a specified error, the number of times that the user want to update the PDF,                    //
//        the path of the txt file that stores the information of the PDF.                                //
// Output: the result of integration                                                                      //
////////////////////////////////////////////////////////////////////////////////////////////////////////////
double Function::NPIS(const int &N, const int &power, const double &err,
                      int NUpdate, std::string path) const
{
    // Set path in which the best PDF is stored
    // The default place is the directory of main.cpp
    // The name of the created document is coefficient.txt
    if (path == " ")
    {
        std::string file_path = __FILE__;
        path = file_path.substr(0, file_path.rfind("\\"));
        path.erase(path.find_last_of("/\\"), std::string::npos);
    
        std::stringstream path_sstream;
        path_sstream << path << "/coefficient.txt";
        path = path_sstream.str();
    }
    
    // If the txt file is empty, have to get a PDF even if NUpdate is 0
    std::ifstream file(path);
    if (file.peek() == std::ifstream::traits_type::eof() && !NUpdate)
        NUpdate = 1;
    
    // Numer of samples uniformly selected from the range
    int NStep = (int)(1000 * range);
    double step = range/NStep;
    
    // AS long as there are more updates need to be done
    while (NUpdate > 0)
    {
        // Get kernel density estimation
        std::function<double(double)> Kernel = Get_Kernel(N);
        
        // Get the coefficients of the fitting polynomial PDF
        std::vector<double> coefficient = Fitting(Kernel, power);
        
        // If the coefficients cannot be solved(rare cases), do the next loop
        if (isnan(coefficient[0]))
        {
            NUpdate--;
            std::cout << "Rejected!" << std::endl;
            // If the file is still empty during the last loop, abort the execution.
            if (file.peek() == std::ifstream::traits_type::eof() && !NUpdate)
                abort();
            continue;
        }
        
        // Construct a polynomial PDF using the obtained coefficients and the given power
        std::function<double(double)> fit = [power, coefficient](double x)
        {   double sum = 0;
            for (int i = power; i >= 0; i--)
                sum += coefficient[power-i] * pow(x, i);
            return sum;
        };
        
        double area = 0;
        double squared_sum = 0;
        int count = 0;
        
        // Iterate over the x intervals
        for (double i = start; i <= end; i += step)
        {   // If fit(i) is not 0 (denominator cannot be 0)
            if (fit(i))
            {   // Get the value of the sample and add relavent numbers to other variables
                double tmp = f(i) / fit(i);
                area += tmp;
                squared_sum += pow(tmp, 2);
                count++;
            }
        }
        // Calculate Epsilon for the selected points
        double Epsilon = sqrt((squared_sum/count - pow(area/count, 2)) / count);
        
        
        // Take last stored coefficients from the txt file
        // Different from the coefficients above!
        std::vector<double> coe;
        std::ifstream stream(path);
        std::copy(std::istream_iterator<double>(stream),
                  std::istream_iterator<double>(),
                  std::back_inserter(coe));
        
        // Stored std
        double pre_Epsilon = 0;
        // If the txt document is not empty, get the value of the last epsilon
        if (coe.size()) pre_Epsilon = coe.back();
        
        // Rewrite the previous file if following conditions are satisfied
        // 1. The file has no data; 2. Newly obtained epsilon is less than pre_std
        // 3. The value of pre_std is nan
        if (!coe.size() || Epsilon < pre_Epsilon || isnan(pre_Epsilon))
        {
            std::ofstream myfile;
            myfile.open(path);
            if (myfile.is_open())
            {   // Write in the coefficients of the polynomial PDF
                for (int i = 0; i < coefficient.size(); i++)
                    myfile << std::fixed << std::setprecision(20) << coefficient[i] << " ";
                // Write in the Epsilon in the end
                myfile << Epsilon;
                myfile.close();
                std::cout << "Accepted! Epsilon: " << Epsilon << std::endl;
            }
            else std::cout << "Unable to open file";
        }
        else std::cout << "Rejected!" << std::endl;
        
        NUpdate--;
        
        // Ask user wether keep updating PDF or not
        if (!NUpdate)
        {
            std::cout << "Do you want to continue optimising the PDF?\n(1 for Yes, 0 for No): ";
            bool YN;
            std::cin >> YN;
            if (YN)
            {   // Ask for number of trials
                std::cout << "How many more times do you want to try?\nEnter a positive integer: ";
                std::cin >> NUpdate;
            }
        }
    }
    std::cout << "Update Completed!" << std::endl;
    
    // Get the coefficients from the txt file
    // This is the best polynomial obtained so far!
    std::vector<double> coefficient;
    std::ifstream stream(path);
    std::copy(std::istream_iterator<double>(stream),
              std::istream_iterator<double>(),
              std::back_inserter(coefficient));
    
    // Get the power of the polynomial PDF(Excluding the constant term and the epsilon)
    long stored_power = coefficient.size() - 2;
    
    // Construct the PDF
    std::function<double(double)> fit = [stored_power, coefficient](double x)
    {   double sum = 0;
        for (long i = stored_power; i >= 0; i--)
            sum += coefficient[stored_power - i] * pow(x, i);
        return sum;
    };
    
    // Call importance sampling with the above PDF and return the result
    return MC(fit, err);
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////
// This function is for Monte Carlo integration with adaptive importance sampling                         //
// The method uses VEGAS algorithm that updates the stratified PDF based on last estimated numerics.      //
// Input: Initial number of bins for stratified sampling, a specified error                               //
// Output: the result of integration                                                                      //
////////////////////////////////////////////////////////////////////////////////////////////////////////////
double Function::AIS(int N, const double &err) const
{
    std::random_device rd;
    std::mt19937 mt(rd());
    
    const int K = 10;        // Parameter that determines the number of sub-sub-divisions
    const int NSample = 10;  // Number of samples need to be drawn in a sub-division
    
    // Sub-divisions
    std::map<double, division> sub;
    
    // Initial width of sub-divisions
    double step = range / N;
    
    // Varible used to determine how many sub-sub-divisions should be created in each sub-division
    double sum_average = 0;
    
    // Create initial sub-divisions
    for (int i = 0; i < N; i++)
    {
        // Partition location
        double loc = end - i * step;
        division _sub;
        
        // Use to take samples in the sub-division
        std::uniform_real_distribution<double> dist(loc - step, loc);
        
        // Take NSamples from the sub-division
        for (int j = 0; j < NSample; j++)
        {
            double x = dist(mt);
            _sub.average += pow(f(x) * range, 2);
        }
        
        _sub.average = sqrt(_sub.average);
        sum_average += _sub.average * step;
        
        // The key for the sub-division is its location
        sub[loc] = _sub;
    }
    
    // Get the partition at the boundary
    division _sub;
    sub[start] = _sub;
    
    double epsilon;
    double area;

    do
    {
        // Number of sub-sub-divisions
        int Ninterval = 0;
        for (std::map<double, division>::iterator iter = std::next(sub.begin()); iter != sub.end(); iter++)
        {
            // Number of sub-divisions within the current sub-division
            int m = nearbyint(K * (iter->second).average / sum_average);
            // Minimum number sub-sub-divisions within the sub-division is 1(sub-division itself)
            (iter->second).NSubsub = m > 0? m:1;
            // Get total number of sub-sub-divisions
            Ninterval += (iter->second).NSubsub;
        }
        
        // Add sub-sub-divisions to sub-divisions
        for (std::map<double, division>::iterator iter = std::next(sub.begin()); iter != sub.end(); iter++)
        {
            double prev_loc = std::prev(iter)->first;        // Previous location
            double delta = iter->first - prev_loc;           // Width of the current sub-division
            double sub_step = delta / (iter->second).NSubsub;// Width for sub-sub-divisions
            
            // Add sub-sub-divisions to sub-divisions
            for (int i = 1; i < (iter->second).NSubsub; i++)
                // Recycle _sub bacause its attributes does not matter(only location matters)
                sub[prev_loc + i * sub_step] = _sub;
        }
        
        // Sub-sub-divisions will be allocated evenly to new sub-divisions
        int residual = Ninterval % N; // Leftovers after allocation
        Ninterval = Ninterval / N;    // Average number of sub-sub-divisions in each sub-division
        
        // New sub-divisions(by regrouping sub-sub-divisions)
        std::map<double, division> new_sub;
        
        // Get new sub-divisions by regrouping sub-sub-divisions
        for (std::map<double, division>::iterator iter = sub.begin();
             iter != std::prev(sub.end(), 1); advance(iter, Ninterval))
        {   // last several elements, allocat residual sub-sub-divisions
            if (new_sub.size() >= N + 1 - residual)
            {   // Jump over one more element
                iter++;
                if (iter == std::prev(sub.end(), 1)) break;
            }
            // Get new locations and reset attributes
            new_sub[iter->first] = iter->second;
            (iter->second).average = 0;
            (iter->second).NSubsub = 0;
        }
        
        // Let the new divisions replace the old one
        sub = new_sub;
        // Add the last location which was not counted in the last loop
        sub[end] = _sub;
        
        
        area = 0;
        sum_average = 0;
        
        // Middle-way variable to calculate the Epsilon
        double condition = 0;
        
        for (std::map<double, division>::iterator iter = std::next(sub.begin()); iter != sub.end(); iter++)
        {
            std::uniform_real_distribution<double> dist(std::prev(iter)->first, iter->first);
            
            // Find the current width of the sub-division
            step = iter->first - std::prev(iter)->first;
            
            double sub_area = 0;
            double squared = 0;
            
            // Take samples and collect essential information to
            // update the widths of sub-divisions/calculate Epsilon/calculate area
            for (int j = 0; j < NSample; j++)
            {
                double x = dist(mt);
                double tmp = f(x) * step * N;
                sub_area += tmp;
                area += tmp;
                squared += pow(tmp, 2);
                (iter->second).average += squared;
            }
            sub_area /= NSample;
            condition += squared / NSample - pow(sub_area, 2);
            (iter->second).average = sqrt((iter->second).average);
            sum_average += (iter->second).average * step;
        }
        // Calculate area
        area = area / N / NSample;
    
        // Calculate epsilon
        epsilon = sqrt(condition / NSample) / N / area;
        
        // Increase bins (The multiplicator is set so that the error is halved for every update)
        N = (int)(N * 1.59);

    } while (epsilon > err);
    
    std::cout << "Epsilon: " << epsilon << "\n" << "NSample: " << N * NSample << std::endl;
    
    return area;
}
