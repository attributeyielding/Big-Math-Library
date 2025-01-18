using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MathBigLibrary
{
    using System;
    using System.Collections.Generic;
    using System.Numerics;

    /// <summary>
    /// Provides a collection of mathematical utility functions and constants.
    /// </summary>
    public static class MathUtils
    {
        // Constants
        public const double Pi = 3.14159265358979323846;
        public const double E = 2.71828182845904523536;

        // Precision for iterative calculations
        public static double Precision { get; set; } = 1e-15;


        public const int MaxIterations = 1000;

        // Custom random number generator
        private static Random _random = new Random();

        // Helper functions
        private static double Abs(double x) => x < 0 ? -x : x;
        private static double Mod(double x, double y) => x - y * Floor(x / y);
        private static double Floor(double x) => x >= 0 ? (long)x : (long)x - 1;

        // =============================================
        // Arithmetic Operations
        // =============================================
        public static class Arithmetic
        {
            /// <summary>
            /// Adds two numbers.
            /// </summary>
            /// <param name="a">The first number.</param>
            /// <param name="b">The second number.</param>
            /// <returns>The sum of a and b.</returns>
            public static double Add(double a, double b) => a + b;

            /// <summary>
            /// Subtracts the second number from the first.
            /// </summary>
            /// <param name="a">The first number.</param>
            /// <param name="b">The second number.</param>
            /// <returns>The difference of a and b.</returns>
            public static double Subtract(double a, double b) => a - b;

            /// <summary>
            /// Multiplies two numbers.
            /// </summary>
            /// <param name="a">The first number.</param>
            /// <param name="b">The second number.</param>
            /// <returns>The product of a and b.</returns>
            public static double Multiply(double a, double b) => a * b;

            /// <summary>
            /// Divides the first number by the second.
            /// </summary>
            /// <param name="a">The numerator.</param>
            /// <param name="b">The denominator.</param>
            /// <returns>The quotient of a and b.</returns>
            /// <exception cref="DivideByZeroException">Thrown when b is zero.</exception>
            public static double Divide(double a, double b)
            {
                if (b == 0) throw new DivideByZeroException("Division by zero is not allowed.");
                return a / b;
            }
        }

        // =============================================
        // Exponential and Logarithmic Functions
        // =============================================
        public static class Exponential
        {
            /// <summary>
            /// Computes the exponential function e^x using a Taylor series expansion.
            /// </summary>
            /// <param name="x">The input value.</param>
            /// <returns>e raised to the Exponential.Power of x.</returns>
            public static double Exp(double x)
            {
                // Handle very large negative values to avoid stack overflow
                if (x < -100) return 0; // e^x for very large negative x is approximately 0

                // Handle negative values iteratively
                if (x < 0) return 1.0 / Exp(-x);

                // Constants
                const double ln2 = 0.6931471805599453; // Natural logarithm of 2

                // Reduce x to a smaller range using the property e^x = e^(n * ln2 + r)
                double n = Floor(x / ln2); // Integer part
                double r = x - n * ln2;    // Remainder

                // Taylor series expansion for e^r
                double term = 1.0, sum = term;
                for (int i = 1; i <= MaxIterations; i++)
                {
                    term *= r / i; // Update term: r^i / i!
                    sum += term;   // Add term to the sum
                    if (Abs(term) < Precision) break; // Stop if the term is small enough
                }

                // If the sum is zero, throw an exception (should not happen)
                if (sum == 0) throw new Exception("Exp function did not converge.");

                // Multiply by 2^n to account for the reduction step
                return sum * Exponential.Pow(2, n);
            }

            /// <summary>
            /// Computes the natural logarithm of x using Newton-Raphson iteration.
            /// </summary>
            /// <param name="x">The input value.</param>
            /// <returns>The natural logarithm of x.</returns>
            /// <exception cref="ArgumentException">Thrown when x <= 0.</exception>
            public static double Log(double x)
            {
                // Handle invalid inputs
                if (x <= 0)
                    throw new ArgumentException("Logarithm of non-positive number is undefined.");

                // Initial guess for y (better approximation for faster convergence)
                double y = x - 1; // Initial guess based on the Taylor series approximation of log(1 + x)

                // Newton-Raphson iteration
                for (int i = 0; i < MaxIterations; i++)
                {
                    double expY = Exp(y); // Compute e^y
                    double nextY = y - (expY - x) / expY; // Update y

                    // Check for convergence
                    if (Abs(nextY - y) < Precision)
                        return nextY;

                    y = nextY; // Update y for the next iteration
                }

                // If the method did not converge, throw an exception
                throw new Exception("Log function did not converge.");
            }

            /// <summary>
            /// Computes baseValue raised to the Exponential.Power of exponent.
            /// </summary>
            /// <param name="baseValue">The base value.</param>
            /// <param name="exponent">The exponent.</param>
            /// <returns>The result of baseValue raised to the Exponential.Power of exponent.</returns>
            /// <exception cref="ArgumentException">Thrown when baseValue < 0 and exponent is not an integer.</exception>
            public static double Pow(double baseValue, double exponent)
            {
                if (baseValue == 0)
                {
                    if (exponent == 0)
                        return double.NaN; // 0^0 is undefined
                    else if (exponent > 0)
                        return 0; // 0^exponent = 0 for exponent > 0
                    else
                        throw new ArgumentException("Negative exponent with zero base is not supported.");
                }

                if (baseValue < 0 && exponent % 1 != 0)
                    throw new ArgumentException("Negative base with non-integer exponent is not supported.");

                double result = 1.0;
                for (int i = 0; i < exponent; i++)
                {
                    result *= baseValue;
                }
                return result;
            }
        }

        // =============================================
        // Trigonometric Functions
        // =============================================
        public static class Trigonometry
        {
            /// <summary>
            /// Computes the sine of x using a Taylor series expansion.
            /// </summary>
            /// <param name="x">The input value in radians.</param>
            /// <returns>The sine of x.</returns>
            public static double Sin(double x)
            {
                x = Mod(x, 2 * Pi);
                if (x < -Pi) x += 2 * Pi;
                else if (x > Pi) x -= 2 * Pi;
                double term = x, sum = term;
                const int MaxIterations = 100;
                for (int n = 3; n <= MaxIterations; n += 2)
                {
                    term *= -x * x / (n * (n - 1));
                    sum += term;
                    if (Abs(term) < Precision) break;
                }
                if (Abs(term) >= Precision) throw new Exception("Sin function did not converge.");
                return sum;
            }

            /// <summary>
            /// Computes the cosine of x using a Taylor series expansion.
            /// </summary>
            /// <param name="x">The input value in radians.</param>
            /// <returns>The cosine of x.</returns>
            public static double Cos(double x)
            {
                x = Mod(x, 2 * Pi);
                if (x < -Pi) x += 2 * Pi;
                else if (x > Pi) x -= 2 * Pi;
                double term = 1.0, sum = term;
                const int MaxIterations = 100;
                for (int n = 2; n <= MaxIterations; n += 2)
                {
                    term *= -x * x / (n * (n - 1));
                    sum += term;
                    if (Abs(term) < Precision) break;
                }
                if (Abs(term) >= Precision) throw new Exception("Cos function did not converge.");
                return sum;
            }

            /// <summary>
            /// Computes the tangent of x.
            /// </summary>
            /// <param name="x">The input value in radians.</param>
            /// <returns>The tangent of x.</returns>
            /// <exception cref="DivideByZeroException">Thrown when cosine of x is zero.</exception>
            public static double Tan(double x)
            {
                double cosValue = Cos(x);
                if (cosValue == 0) throw new DivideByZeroException("Tangent is undefined at this point.");
                return Sin(x) / cosValue;
            }

            /// <summary>
            /// Computes the arcsine of x.
            /// </summary>
            /// <param name="x">The input value.</param>
            /// <returns>The arcsine of x in radians.</returns>
            /// <exception cref="ArgumentException">Thrown when |x| > 1.</exception>
            public static double ASin(double x)
            {
                if (x < -1 || x > 1) throw new ArgumentException("ASin is undefined for |x| > 1.");
                return ATan(x / Roots.Sqrt(1 - x * x));
            }

            /// <summary>
            /// Computes the arccosine of x.
            /// </summary>
            /// <param name="x">The input value.</param>
            /// <returns>The arccosine of x in radians.</returns>
            /// <exception cref="ArgumentException">Thrown when |x| > 1.</exception>
            public static double ACos(double x)
            {
                if (x < -1 || x > 1) throw new ArgumentException("ACos is undefined for |x| > 1.");
                return Pi / 2 - ASin(x);
            }

            /// <summary>
            /// Computes the arctangent of x.
            /// </summary>
            /// <param name="x">The input value.</param>
            /// <returns>The arctangent of x in radians.</returns>
            public static double ATan(double x)
            {
                if (x < -1 || x > 1)
                {
                    if (x > 0) return Pi / 2 - ATan(1 / x);
                    else return -Pi / 2 - ATan(1 / x);
                }
                double t = x, sum = t, term = t;
                const int MaxIterations = 100;
                for (int n = 1; n <= MaxIterations; n++)
                {
                    term *= -x * x * x * x / ((2 * n + 1) * (2 * n + 1));
                    sum += term / (2 * n + 1);
                    if (Abs(term) < Precision) break;
                }
                if (Abs(term) >= Precision) throw new Exception("ATan function did not converge.");
                return sum;
            }
        }

        // =============================================
        // Hyperbolic Functions
        // =============================================
        public static class Hyperbolic
        {
            /// <summary>
            /// Computes the hyperbolic sine of x.
            /// </summary>
            /// <param name="x">The input value.</param>
            /// <returns>The hyperbolic sine of x.</returns>
            public static double Sinh(double x) => (Exponential.Exp(x) - Exponential.Exp(-x)) / 2;

            /// <summary>
            /// Computes the hyperbolic cosine of x.
            /// </summary>
            /// <param name="x">The input value.</param>
            /// <returns>The hyperbolic cosine of x.</returns>
            public static double Cosh(double x) => (Exponential.Exp(x) + Exponential.Exp(-x)) / 2;

            /// <summary>
            /// Computes the hyperbolic tangent of x.
            /// </summary>
            /// <param name="x">The input value.</param>
            /// <returns>The hyperbolic tangent of x.</returns>
            /// <exception cref="DivideByZeroException">Thrown when cosh(x) is zero.</exception>
            public static double Tanh(double x)
            {
                double coshValue = Cosh(x);
                if (coshValue == 0) throw new DivideByZeroException("Tanh is undefined at this point.");
                return Sinh(x) / coshValue;
            }

            /// <summary>
            /// Computes the inverse hyperbolic sine of x.
            /// </summary>
            /// <param name="x">The input value.</param>
            /// <returns>The inverse hyperbolic sine of x.</returns>
            public static double ASinh(double x) => Exponential.Log(x + Roots.Sqrt(x * x + 1));

            /// <summary>
            /// Computes the inverse hyperbolic cosine of x.
            /// </summary>
            /// <param name="x">The input value.</param>
            /// <returns>The inverse hyperbolic cosine of x.</returns>
            /// <exception cref="ArgumentException">Thrown when x < 1.</exception>
            public static double ACosh(double x)
            {
                if (x < 1) throw new ArgumentException("ACosh is undefined for x < 1.");
                return Exponential.Log(x + Roots.Sqrt(x * x - 1));
            }

            /// <summary>
            /// Computes the inverse hyperbolic tangent of x.
            /// </summary>
            /// <param name="x">The input value.</param>
            /// <returns>The inverse hyperbolic tangent of x.</returns>
            /// <exception cref="ArgumentException">Thrown when |x| >= 1.</exception>
            public static double ATanh(double x)
            {
                if (Abs(x) >= 1) throw new ArgumentException("ATanh is undefined for |x| >= 1.");
                return 0.5 * Exponential.Log((1 + x) / (1 - x));
            }
        }

        // =============================================
        // Root and Exponential.Power Functions
        // =============================================
        public static class Roots
        {
            /// <summary>
            /// Computes the square root of x using Newton-Raphson iteration.
            /// </summary>
            /// <param name="x">The input value.</param>
            /// <returns>The square root of x.</returns>
            /// <exception cref="ArgumentException">Thrown when x < 0.</exception>
            public static double Sqrt(double x)
            {
                if (x < 0) throw new ArgumentException("Square root of negative number is undefined.");
                if (x == 0) return 0;
                double guess = x / 2;
                const int MaxIterations = 100;
                for (int i = 0; i < MaxIterations; i++)
                {
                    double nextGuess = (guess + x / guess) / 2;
                    if (Abs(nextGuess - guess) < Precision) return nextGuess;
                    guess = nextGuess;
                }
                throw new Exception("Sqrt function did not converge.");
            }

            /// <summary>
            /// Computes the cube root of x using Newton-Raphson iteration.
            /// </summary>
            /// <param name="x">The input value.</param>
            /// <returns>The cube root of x.</returns>
            public static double Cbrt(double x)
            {
                if (x < 0) return -Cbrt(-x);
                double guess = x / 3;
                const int MaxIterations = 100;
                for (int i = 0; i < MaxIterations; i++)
                {
                    double nextGuess = (2 * guess + x / (guess * guess)) / 3;
                    if (Abs(nextGuess - guess) < Precision) return nextGuess;
                    guess = nextGuess;
                }
                throw new Exception("Cbrt function did not converge.");
            }
        }

        // =============================================
        // Probability and Statistics
        // =============================================
        public static class Statistics
        {
            /// <summary>
            /// Computes the mean of an array of values.
            /// </summary>
            /// <param name="values">The array of values.</param>
            /// <returns>The mean of the values.</returns>
            /// <exception cref="ArgumentException">Thrown when the array is null or empty.</exception>
            public static double Mean(double[] values)
            {
                if (values == null || values.Length == 0) throw new ArgumentException("Array cannot be null or empty.");
                double sum = 0;
                foreach (double value in values) sum += value;
                return sum / values.Length;
            }

            /// <summary>
            /// Computes the variance of an array of values.
            /// </summary>
            /// <param name="values">The array of values.</param>
            /// <returns>The variance of the values.</returns>
            /// <exception cref="ArgumentException">Thrown when the array is null or empty.</exception>
            public static double Variance(double[] values)
            {
                if (values == null || values.Length == 0) throw new ArgumentException("Array cannot be null or empty.");
                double mean = Mean(values);
                double sum = 0;
                foreach (double value in values) sum += Exponential.Pow(value - mean, 2);
                return sum / values.Length;
            }

            /// <summary>
            /// Computes the standard deviation of an array of values.
            /// </summary>
            /// <param name="values">The array of values.</param>
            /// <returns>The standard deviation of the values.</returns>
            /// <exception cref="ArgumentException">Thrown when the array is null or empty.</exception>
            public static double StandardDeviation(double[] values) => Roots.Sqrt(Variance(values));

            /// <summary>
            /// Computes the probability density function (PDF) of the normal distribution.
            /// </summary>
            /// <param name="x">The input value.</param>
            /// <param name="mean">The mean of the distribution.</param>
            /// <param name="stdDev">The standard deviation of the distribution.</param>
            /// <returns>The PDF value at x.</returns>
            public static double NormalPDF(double x, double mean = 0, double stdDev = 1) =>
                Exponential.Exp(-0.5 * Exponential.Pow((x - mean) / stdDev, 2)) / (stdDev * Roots.Sqrt(2 * Pi));

            /// <summary>
            /// Computes the cumulative distribution function (CDF) of the normal distribution.
            /// </summary>
            /// <param name="x">The input value.</param>
            /// <param name="mean">The mean of the distribution.</param>
            /// <param name="stdDev">The standard deviation of the distribution.</param>
            /// <returns>The CDF value at x.</returns>
            public static double NormalCDF(double x, double mean = 0, double stdDev = 1) =>
                0.5 * (1 + SpecialFunctions.Erf((x - mean) / (stdDev * Roots.Sqrt(2))));
        }

        // =============================================
        // Special Mathematical Functions
        // =============================================
        public static class SpecialFunctions
        {
            private static Dictionary<int, double> factorialCache = new Dictionary<int, double>();

            /// <summary>
            /// Computes the factorial of n iteratively and caches results.
            /// </summary>
            /// <param name="n">The input value.</param>
            /// <returns>n!.</returns>
            /// <exception cref="ArgumentException">Thrown when n < 0.</exception>
            public static double Factorial(int n)
            {
                if (n < 0) throw new ArgumentException("Factorial is not defined for negative numbers.");
                if (n == 0 || n == 1) return 1;
                if (factorialCache.ContainsKey(n)) return factorialCache[n];
                double result = 1;
                for (int i = 1; i <= n; i++)
                {
                    result *= i;
                    if (!factorialCache.ContainsKey(i)) factorialCache.Add(i, result);
                }
                return result;
            }

            /// <summary>
            /// Computes the Gamma function using the Lanczos approximation.
            /// </summary>
            /// <param name="x">The input value.</param>
            /// <returns>The Gamma function value at x.</returns>
            /// <exception cref="ArgumentException">Thrown when x is a non-positive integer.</exception>
            public static double Gamma(double x)
            {
                if (x <= 0 && x % 1 == 0) throw new ArgumentException("Gamma function is undefined for non-positive integers.");
                if (x > 20) return Roots.Sqrt(2 * Pi / x) * Exponential.Pow(x / E, x);
                double[] p = { 676.5203681218851, -1259.1392167224028, 771.32342877765313, -176.61502916214059, 12.507343278686905, -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7 };
                double g = 7;
                if (x < 0.5) return Pi / (Trigonometry.Sin(Pi * x) * Gamma(1 - x));
                x -= 1;
                double a = 0.99999999999980993;
                double t = x + g + 0.5;
                for (int i = 0; i < p.Length; i++) a += p[i] / (x + i + 1);
                return Roots.Sqrt(2 * Pi) * Exponential.Pow(t, x + 0.5) * Exponential.Exp(-t) * a;
            }

            /// <summary>
            /// Computes the Beta function in terms of the Gamma function.
            /// </summary>
            /// <param name="x">The first parameter.</param>
            /// <param name="y">The second parameter.</param>
            /// <returns>The Beta function value at (x, y).</returns>
            public static double Beta(double x, double y) => Gamma(x) * Gamma(y) / Gamma(x + y);

            /// <summary>
            /// Computes the error function using a rational approximation.
            /// </summary>
            /// <param name="x">The input value.</param>
            /// <returns>The error function value at x.</returns>
            public static double Erf(double x)
            {
                double t = 1.0 / (1.0 + 0.3275911 * Abs(x));
                double tau = t * (0.254829592 + t * (-0.284496736 + t * (1.421413741 + t * (-1.453152027 + t * 1.061405429))));
                double erf = 1.0 - tau * Exponential.Exp(-x * x);
                return x >= 0 ? erf : -erf;
            }

            /// <summary>
            /// Computes the Gaussian hypergeometric function 2F1 using a series expansion.
            /// </summary>
            /// <param name="a">The first parameter.</param>
            /// <param name="b">The second parameter.</param>
            /// <param name="c">The third parameter.</param>
            /// <param name="x">The input value.</param>
            /// <returns>The hypergeometric function value at (a, b, c, x).</returns>
            /// <exception cref="ArgumentException">Thrown when |x| >= 1.</exception>
            public static double Hypergeometric2F1(double a, double b, double c, double x)
            {
                if (Abs(x) >= 1) throw new ArgumentException("Hypergeometric function is only defined for |x| < 1.");
                double term = 1.0, sum = term;
                const int MaxIterations = 1000;
                for (int n = 1; n <= MaxIterations; n++)
                {
                    term *= (a + n - 1) * (b + n - 1) * x / (c + n - 1) / n;
                    sum += term;
                    if (Abs(term) < Precision) break;
                }
                if (Abs(term) >= Precision) throw new Exception("Hypergeometric2F1 function did not converge.");
                return sum;
            }

            /// <summary>
            /// Computes the Bessel function of the first kind J_n(x) using a series expansion.
            /// </summary>
            /// <param name="n">The order of the Bessel function.</param>
            /// <param name="x">The input value.</param>
            /// <returns>The Bessel function value J_n(x).</returns>
            public static double BesselJ(int n, double x)
            {
                double sum = 0;
                double term = 0;

                for (int k = 0; k <= MaxIterations; k++)
                {
                    term = Math.Pow(-1, k) * Math.Pow(x / 2, 2 * k + n) / (Factorial(k) * Factorial(k + n));
                    sum += term;

                    if (Math.Abs(term) < Precision)
                        break;
                }

                if (Math.Abs(term) >= Precision)
                    throw new Exception("BesselJ function did not converge.");

                return sum;
            }

            /// <summary>
            /// Computes the Bessel function of the second kind Y_n(x) using a series expansion.
            /// </summary>
            /// <param name="n">The order of the Bessel function.</param>
            /// <param name="x">The input value.</param>
            /// <returns>The Bessel function value Y_n(x).</returns>
            /// <exception cref="ArgumentException">Thrown when x <= 0.</exception>
            public static double BesselY(int n, double x)
            {
                // Handle invalid input (x <= 0)
                if (x <= 0)
                    throw new ArgumentException("BesselY is undefined for non-positive x.");

                // Constants
                const double EulerMascheroni = 0.5772156649015329; // Euler-Mascheroni constant

                // Compute J_n(x) using the existing BesselJ method
                double Jn = BesselJ(n, x);

                // Compute the first term: (2/π) * (ln(x/2) + γ) * J_n(x)
                double term1 = (2 / MathUtils.Pi) * (Math.Log(x / 2) + EulerMascheroni) * Jn;

                // Compute the second term: sum over k
                double sum = 0;
                for (int k = 0; k <= 100; k++) // Limit iterations to 100 for practical purposes
                {
                    // Compute the numerator: (-1)^k * (x/2)^(n + 2k)
                    double numerator = Math.Pow(-1, k) * Math.Pow(x / 2, n + 2 * k);

                    // Compute the denominator: k! * (n + k)!
                    double denominator = Factorial(k) * Factorial(n + k);

                    // Compute the digamma terms: ψ(n + k + 1) + ψ(k + 1)
                    double digammaNk = Digamma(n + k + 1);
                    double digammaK = Digamma(k + 1);

                    // Add the term to the sum
                    sum += numerator / denominator * (digammaNk + digammaK);

                    // Stop if the term is small enough
                    if (Math.Abs(numerator / denominator) < MathUtils.Precision)
                        break;
                }

                // Compute the second term: (-1/π) * sum
                double term2 = (-1 / MathUtils.Pi) * sum;

                // Return the result: term1 + term2
                return term1 + term2;
            }
            /// <summary>
            /// Computes the digamma function (ψ) of a given value x.
            /// The digamma function is the logarithmic derivative of the gamma function.
            /// </summary>
            /// <param name="x">The input value. Must be a positive number.</param>
            /// <returns>The value of the digamma function at x.</returns>
            /// <exception cref="ArgumentException">Thrown when x is less than or equal to 0, as the digamma function is undefined for non-positive values.</exception>
            /// <remarks>
            /// The digamma function is computed using the following approximation for large x:
            /// <code>
            /// ψ(x) ≈ ln(x) - 1/(2x) - 1/(12x^2) + 1/(120x^4) - 1/(252x^6) + ...
            /// </code>
            /// For small x, the recurrence relation ψ(x) = ψ(x + 1) - 1/x is used.
            /// </remarks>
            /// <example>
            /// <code>
            /// double result = Digamma(5.0); // Computes ψ(5)
            /// </code>
            /// </example>
            public static double Digamma(double x)
            {
                if (x <= 0)
                    throw new ArgumentException("Digamma is undefined for non-positive x.");

                // Use an approximation for large x
                if (x > 6)
                {
                    double result = Math.Log(x) - 1 / (2 * x) - 1 / (12 * x * x) + 1 / (120 * x * x * x * x);
                    return result;
                }

                // Use the recurrence relation ψ(x) = ψ(x + 1) - 1/x for small x
                return Digamma(x + 1) - 1 / x;
            }
        }

        // =============================================
        // Numerical Methods
        // =============================================
        public static class NumericalMethods
        {
            /// <summary>
            /// Computes the integral of f from a to b using Simpson's rule.
            /// </summary>
            /// <param name="f">The integrand function.</param>
            /// <param name="a">The lower limit of integration.</param>
            /// <param name="b">The upper limit of integration.</param>
            /// <param name="n">The number of intervals (must be even).</param>
            /// <returns>The approximate integral of f from a to b.</returns>
            public static double Integrate(Func<double, double> f, double a, double b, int n = 1000)
            {
                if (n % 2 != 0) n += 1; // Ensure n is even
                double h = (b - a) / n;
                double sum = f(a) + f(b);
                for (int i = 1; i < n; i++)
                {
                    double x = a + i * h;
                    sum += (i % 2 == 0 ? 2 : 4) * f(x);
                }
                return sum * h / 3;
            }

            /// <summary>
            /// Finds the root of a function using the Newton-Raphson method.
            /// </summary>
            /// <param name="f">The function.</param>
            /// <param name="df">The derivative of the function.</param>
            /// <param name="initialGuess">The initial guess for the root.</param>
            /// <returns>The approximate root of the function.</returns>
            /// <exception cref="InvalidOperationException">Thrown when the derivative is too small.</exception>
            public static double FindRoot(Func<double, double> f, Func<double, double> df, double initialGuess)
            {
                double x = initialGuess;
                const int MaxIterations = 100;
                for (int i = 0; i < MaxIterations; i++)
                {
                    double fx = f(x);
                    double dfx = df(x);
                    if (Math.Abs(dfx) < Precision)
                        throw new InvalidOperationException("Derivative is too small; cannot continue.");
                    double nextX = x - fx / dfx;
                    if (Math.Abs(nextX - x) < Precision) return nextX;
                    x = nextX;
                }
                throw new Exception("FindRoot function did not converge.");
            }
        }

        // =============================================
        // Number Theory Functions
        // =============================================
        public static class NumberTheory
        {
            /// <summary>
            /// Computes the greatest common divisor (GCD) of two numbers using the Euclidean algorithm.
            /// </summary>
            /// <param name="a">The first number.</param>
            /// <param name="b">The second number.</param>
            /// <returns>The GCD of a and b.</returns>
            public static int GCD(int a, int b)
            {
                while (b != 0)
                {
                    int temp = b;
                    b = a % b;
                    a = temp;
                }
                return a;
            }

            /// <summary>
            /// Computes the least common multiple (LCM) of two numbers.
            /// </summary>
            /// <param name="a">The first number.</param>
            /// <param name="b">The second number.</param>
            /// <returns>The LCM of a and b.</returns>
            public static int LCM(int a, int b) => a / GCD(a, b) * b;

            /// <summary>
            /// Checks if a number is prime using trial division.
            /// </summary>
            /// <param name="n">The number to check.</param>
            /// <returns>True if the number is prime; otherwise, false.</returns>
            public static bool IsPrime(int n)
            {
                if (n <= 1) return false;
                if (n <= 3) return true;
                if (n % 2 == 0 || n % 3 == 0) return false;
                for (int i = 5; i * i <= n; i += 6)
                    if (n % i == 0 || n % (i + 2) == 0) return false;
                return true;
            }

            /// <summary>
            /// Computes the nth Fibonacci number iteratively.
            /// </summary>
            /// <param name="n">The index of the Fibonacci number.</param>
            /// <returns>The nth Fibonacci number.</returns>
            public static int Fibonacci(int n)
            {
                if (n <= 0) return 0;
                if (n == 1) return 1;
                int a = 0, b = 1;
                for (int i = 2; i <= n; i++)
                {
                    int temp = a + b;
                    a = b;
                    b = temp;
                }
                return b;
            }
        }

        // =============================================
        // Sorting and Searching Algorithms
        // =============================================
        public static class Algorithms
        {
            /// <summary>
            /// Sorts an array using the bubble sort algorithm.
            /// </summary>
            /// <param name="arr">The array to sort.</param>
            public static void BubbleSort(int[] arr)
            {
                int n = arr.Length;
                for (int i = 0; i < n - 1; i++)
                    for (int j = 0; j < n - i - 1; j++)
                        if (arr[j] > arr[j + 1])
                        {
                            int temp = arr[j];
                            arr[j] = arr[j + 1];
                            arr[j + 1] = temp;
                        }
            }

            /// <summary>
            /// Sorts an array using the quicksort algorithm.
            /// </summary>
            /// <param name="arr">The array to sort.</param>
            /// <param name="low">The starting index.</param>
            /// <param name="high">The ending index.</param>
            public static void QuickSort(int[] arr, int low, int high)
            {
                if (low < high)
                {
                    int pi = Partition(arr, low, high);
                    QuickSort(arr, low, pi - 1);
                    QuickSort(arr, pi + 1, high);
                }
            }

            private static int Partition(int[] arr, int low, int high)
            {
                int pivot = arr[high];
                int i = low - 1;
                for (int j = low; j < high; j++)
                    if (arr[j] < pivot)
                    {
                        i++;
                        int temp = arr[i];
                        arr[i] = arr[j];
                        arr[j] = temp;
                    }
                int temp1 = arr[i + 1];
                arr[i + 1] = arr[high];
                arr[high] = temp1;
                return i + 1;
            }

            /// <summary>
            /// Searches for a target value in an array using linear search.
            /// </summary>
            /// <param name="arr">The array to search.</param>
            /// <param name="target">The target value.</param>
            /// <returns>The index of the target value, or -1 if not found.</returns>
            public static int LinearSearch(int[] arr, int target)
            {
                for (int i = 0; i < arr.Length; i++)
                    if (arr[i] == target) return i;
                return -1;
            }

            /// <summary>
            /// Searches for a target value in a sorted array using binary search.
            /// </summary>
            /// <param name="arr">The sorted array to search.</param>
            /// <param name="target">The target value.</param>
            /// <returns>The index of the target value, or -1 if not found.</returns>
            public static int BinarySearch(int[] arr, int target)
            {
                int low = 0, high = arr.Length - 1;
                while (low <= high)
                {
                    int mid = low + (high - low) / 2;
                    if (arr[mid] == target) return mid;
                    if (arr[mid] < target) low = mid + 1;
                    else high = mid - 1;
                }
                return -1;
            }
        }

        // =============================================
        // Matrix Functions
        // =============================================
        public static class Matrix
        {
            /// <summary>
            /// Adds two matrices.
            /// </summary>
            /// <param name="A">The first matrix.</param>
            /// <param name="B">The second matrix.</param>
            /// <returns>The sum of A and B.</returns>
            /// <exception cref="ArgumentException">Thrown when matrices have different dimensions.</exception>
            public static double[,] Add(double[,] A, double[,] B)
            {
                int rows = A.GetLength(0), cols = A.GetLength(1);
                if (rows != B.GetLength(0) || cols != B.GetLength(1))
                    throw new ArgumentException("Matrices must have the same dimensions.");
                double[,] result = new double[rows, cols];
                for (int i = 0; i < rows; i++)
                    for (int j = 0; j < cols; j++)
                        result[i, j] = A[i, j] + B[i, j];
                return result;
            }

            /// <summary>
            /// Multiplies two matrices.
            /// </summary>
            /// <param name="A">The first matrix.</param>
            /// <param name="B">The second matrix.</param>
            /// <returns>The product of A and B.</returns>
            /// <exception cref="ArgumentException">Thrown when the number of columns in A does not equal the number of rows in B.</exception>
            public static double[,] Multiply(double[,] A, double[,] B)
            {
                int rowsA = A.GetLength(0), colsA = A.GetLength(1);
                int rowsB = B.GetLength(0), colsB = B.GetLength(1);
                if (colsA != rowsB)
                    throw new ArgumentException("Number of columns in A must equal number of rows in B.");
                double[,] result = new double[rowsA, colsB];
                for (int i = 0; i < rowsA; i++)
                    for (int j = 0; j < colsB; j++)
                        for (int k = 0; k < colsA; k++)
                            result[i, j] += A[i, k] * B[k, j];
                return result;
            }

            /// <summary>
            /// Computes the determinant of a 2x2 matrix.
            /// </summary>
            /// <param name="matrix">The 2x2 matrix.</param>
            /// <returns>The determinant of the matrix.</returns>
            /// <exception cref="ArgumentException">Thrown when the matrix is not 2x2.</exception>
            public static double Determinant2x2(double[,] matrix)
            {
                if (matrix.GetLength(0) != 2 || matrix.GetLength(1) != 2)
                    throw new ArgumentException("Matrix must be 2x2.");
                return matrix[0, 0] * matrix[1, 1] - matrix[0, 1] * matrix[1, 0];
            }

            /// <summary>
            /// Computes the determinant of a square matrix using recursion.
            /// </summary>
            /// <param name="matrix">The square matrix.</param>
            /// <returns>The determinant of the matrix.</returns>
            /// <exception cref="ArgumentException">Thrown when the matrix is not square.</exception>
            public static double Determinant(double[,] matrix)
            {
                int n = matrix.GetLength(0);
                if (n != matrix.GetLength(1))
                    throw new ArgumentException("Matrix must be square.");
                if (n == 2) return Determinant2x2(matrix);
                double det = 0;
                for (int col = 0; col < n; col++)
                {
                    double[,] submatrix = new double[n - 1, n - 1];
                    for (int i = 1; i < n; i++)
                    {
                        int subcol = 0;
                        for (int j = 0; j < n; j++)
                        {
                            if (j == col) continue;
                            submatrix[i - 1, subcol] = matrix[i, j];
                            subcol++;
                        }
                    }
                    det += matrix[0, col] * (col % 2 == 0 ? 1 : -1) * Determinant(submatrix);
                }
                return det;
            }

            /// <summary>
            /// Computes the inverse of a 2x2 matrix.
            /// </summary>
            /// <param name="matrix">The 2x2 matrix.</param>
            /// <returns>The inverse of the matrix.</returns>
            /// <exception cref="InvalidOperationException">Thrown when the matrix is singular.</exception>
            public static double[,] Inverse2x2(double[,] matrix)
            {
                double det = Determinant2x2(matrix);
                if (det == 0)
                    throw new InvalidOperationException("Matrix is singular and cannot be inverted.");
                double[,] inverse = new double[2, 2];
                inverse[0, 0] = matrix[1, 1] / det;
                inverse[0, 1] = -matrix[0, 1] / det;
                inverse[1, 0] = -matrix[1, 0] / det;
                inverse[1, 1] = matrix[0, 0] / det;
                return inverse;
            }
        }

        // =============================================
        // Geometry Functions
        // =============================================
        public static class Geometry
        {
            /// <summary>
            /// Computes the Euclidean distance between two points in 2D space.
            /// </summary>
            /// <param name="x1">The x-coordinate of the first point.</param>
            /// <param name="y1">The y-coordinate of the first point.</param>
            /// <param name="x2">The x-coordinate of the second point.</param>
            /// <param name="y2">The y-coordinate of the second point.</param>
            /// <returns>The Euclidean distance between the two points.</returns>
            public static double Distance2D(double x1, double y1, double x2, double y2)
            {
                return Roots.Sqrt(Exponential.Pow(x2 - x1, 2) + Exponential.Pow(y2 - y1, 2));
            }

            /// <summary>
            /// Computes the Euclidean distance between two points in 3D space.
            /// </summary>
            /// <param name="x1">The x-coordinate of the first point.</param>
            /// <param name="y1">The y-coordinate of the first point.</param>
            /// <param name="z1">The z-coordinate of the first point.</param>
            /// <param name="x2">The x-coordinate of the second point.</param>
            /// <param name="y2">The y-coordinate of the second point.</param>
            /// <param name="z2">The z-coordinate of the second point.</param>
            /// <returns>The Euclidean distance between the two points.</returns>
            public static double Distance3D(double x1, double y1, double z1, double x2, double y2, double z2)
            {
                return Math.Sqrt(Exponential.Pow(x2 - x1, 2) + Exponential.Pow(y2 - y1, 2) + Exponential.Pow(z2 - z1, 2));
            }
        }

        // =============================================
        // Complex Numbers
        // =============================================
        public struct Complex
        {
            public double Real { get; }
            public double Imaginary { get; }

            public Complex(double real, double imaginary)
            {
                Real = real;
                Imaginary = imaginary;
            }

            public static Complex operator +(Complex a, Complex b)
            {
                return new Complex(a.Real + b.Real, a.Imaginary + b.Imaginary);
            }

            public static Complex operator *(Complex a, Complex b)
            {
                return new Complex(a.Real * b.Real - a.Imaginary * b.Imaginary, a.Real * b.Imaginary + a.Imaginary * b.Real);
            }

            public override string ToString()
            {
                return $"{Real} + {Imaginary}i";
            }
        }

        // =============================================
        // Random Number Generation
        // =============================================
        /// <summary>
        /// Generates a random double between 0 and 1.
        /// </summary>
        /// <returns>A random double between 0 and 1.</returns>
        public static double Random() => _random.NextDouble();

        /// <summary>
        /// Generates a random integer within the specified range.
        /// </summary>
        /// <param name="minValue">The minimum value.</param>
        /// <param name="maxValue">The maximum value.</param>
        /// <returns>A random integer between minValue and maxValue.</returns>
        public static int RandomInt(int minValue, int maxValue) => _random.Next(minValue, maxValue);

        /// <summary>
        /// Reseeds the random number generator.
        /// </summary>
        /// <param name="seed">The new seed value.</param>
        public static void ReseedRandom(long seed) => _random = new Random((int)seed);
    }
}
