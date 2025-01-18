using System;
using MathBigLibrary;

namespace MathBigLibraryTests
{
    class Program
    {
        static void Main(string[] args)
        {
            Console.WriteLine("Running tests for MathBigLibrary...\n");

            // Test Arithmetic Operations
            TestArithmetic();

            // Test Exponential and Logarithmic Functions
            TestExponential();

            // Test Trigonometric Functions
            TestTrigonometry();

            // Test Hyperbolic Functions
            TestHyperbolic();

            // Test Root and Power Functions
            TestRoots();

            // Test Probability and Statistics
            TestStatistics();

            // Test Special Mathematical Functions
            TestSpecialFunctions();

            // Test Numerical Methods
            TestNumericalMethods();

            // Test Number Theory Functions
            TestNumberTheory();

            // Test Sorting and Searching Algorithms
            TestAlgorithms();

            // Test Matrix Functions
            TestMatrix();

            // Test Geometry Functions
            TestGeometry();

            // Test Complex Numbers
            TestComplexNumbers();

            // Test Random Number Generation
            TestRandom();

            Console.WriteLine("\nAll tests completed.");
        }

        static void TestArithmetic()
        {
            Console.WriteLine("Testing Arithmetic Operations...");
            double a = 10, b = 4;

            Console.WriteLine($"Add({a}, {b}) = {MathUtils.Arithmetic.Add(a, b)}");
            Console.WriteLine($"Subtract({a}, {b}) = {MathUtils.Arithmetic.Subtract(a, b)}");
            Console.WriteLine($"Multiply({a}, {b}) = {MathUtils.Arithmetic.Multiply(a, b)}");
            Console.WriteLine($"Divide({a}, {b}) = {MathUtils.Arithmetic.Divide(a, b)}");

            try
            {
                Console.WriteLine($"Divide({a}, 0) = {MathUtils.Arithmetic.Divide(a, 0)}");
            }
            catch (DivideByZeroException ex)
            {
                Console.WriteLine($"Divide({a}, 0) threw exception: {ex.Message}");
            }
            Console.WriteLine();
        }

        static void TestExponential()
        {
            Console.WriteLine("Testing Exponential and Logarithmic Functions...");
            double x = 2;

            Console.WriteLine($"Exp({x}) = {MathUtils.Exponential.Exp(x)}");
            Console.WriteLine($"Log({x}) = {MathUtils.Exponential.Log(x)}");
            Console.WriteLine($"Pow({x}, 3) = {MathUtils.Exponential.Pow(x, 3)}");

            try
            {
                Console.WriteLine($"Log(-1) = {MathUtils.Exponential.Log(-1)}");
            }
            catch (ArgumentException ex)
            {
                Console.WriteLine($"Log(-1) threw exception: {ex.Message}");
            }
            Console.WriteLine();
        }

        static void TestTrigonometry()
        {
            Console.WriteLine("Testing Trigonometric Functions...");
            double angle = MathUtils.Pi / 4; // 45 degrees

            Console.WriteLine($"Sin({angle}) = {MathUtils.Trigonometry.Sin(angle)}");
            Console.WriteLine($"Cos({angle}) = {MathUtils.Trigonometry.Cos(angle)}");
            Console.WriteLine($"Tan({angle}) = {MathUtils.Trigonometry.Tan(angle)}");
            Console.WriteLine($"ASin(0.7071) = {MathUtils.Trigonometry.ASin(0.7071)}");
            Console.WriteLine($"ACos(0.7071) = {MathUtils.Trigonometry.ACos(0.7071)}");
            Console.WriteLine($"ATan(1) = {MathUtils.Trigonometry.ATan(1)}");
            Console.WriteLine();
        }

        static void TestHyperbolic()
        {
            Console.WriteLine("Testing Hyperbolic Functions...");
            double x = 1;

            Console.WriteLine($"Sinh({x}) = {MathUtils.Hyperbolic.Sinh(x)}");
            Console.WriteLine($"Cosh({x}) = {MathUtils.Hyperbolic.Cosh(x)}");
            Console.WriteLine($"Tanh({x}) = {MathUtils.Hyperbolic.Tanh(x)}");
            Console.WriteLine($"ASinh({x}) = {MathUtils.Hyperbolic.ASinh(x)}");
            Console.WriteLine($"ACosh({x + 1}) = {MathUtils.Hyperbolic.ACosh(x + 1)}");
            Console.WriteLine($"ATanh(0.5) = {MathUtils.Hyperbolic.ATanh(0.5)}");
            Console.WriteLine();
        }

        static void TestRoots()
        {
            Console.WriteLine("Testing Root and Power Functions...");
            double x = 16;

            Console.WriteLine($"Sqrt({x}) = {MathUtils.Roots.Sqrt(x)}");
            Console.WriteLine($"Cbrt({x}) = {MathUtils.Roots.Cbrt(x)}");
            Console.WriteLine();
        }

        static void TestStatistics()
        {
            Console.WriteLine("Testing Probability and Statistics...");
            double[] values = { 1, 2, 3, 4, 5 };

            Console.WriteLine($"Mean = {MathUtils.Statistics.Mean(values)}");
            Console.WriteLine($"Variance = {MathUtils.Statistics.Variance(values)}");
            Console.WriteLine($"Standard Deviation = {MathUtils.Statistics.StandardDeviation(values)}");
            Console.WriteLine($"NormalPDF(0) = {MathUtils.Statistics.NormalPDF(0)}");
            Console.WriteLine($"NormalCDF(0) = {MathUtils.Statistics.NormalCDF(0)}");
            Console.WriteLine();
        }

        static void TestSpecialFunctions()
        {
            Console.WriteLine("Testing Special Mathematical Functions...");
            int n = 5;
            double x = 2.5;

            Console.WriteLine($"Factorial({n}) = {MathUtils.SpecialFunctions.Factorial(n)}");
            Console.WriteLine($"Gamma({x}) = {MathUtils.SpecialFunctions.Gamma(x)}");
            Console.WriteLine($"Beta({x}, {x}) = {MathUtils.SpecialFunctions.Beta(x, x)}");
            Console.WriteLine($"Erf({x}) = {MathUtils.SpecialFunctions.Erf(x)}");
            Console.WriteLine($"Hypergeometric2F1(1, 1, 1, 0.5) = {MathUtils.SpecialFunctions.Hypergeometric2F1(1, 1, 1, 0.5)}");
            Console.WriteLine($"BesselJ(1, {x}) = {MathUtils.SpecialFunctions.BesselJ(1, x)}");
            Console.WriteLine($"BesselY(1, {x}) = {MathUtils.SpecialFunctions.BesselY(1, x)}");
            Console.WriteLine();
        }

        static void TestNumericalMethods()
        {
            Console.WriteLine("Testing Numerical Methods...");
            Func<double, double> f = x => x * x;
            double a = 0, b = 1;

            Console.WriteLine($"Integrate(x^2, {a}, {b}) = {MathUtils.NumericalMethods.Integrate(f, a, b)}");
            Console.WriteLine($"FindRoot(x^2 - 4, 2) = {MathUtils.NumericalMethods.FindRoot(x => x * x - 4, x => 2 * x, 2)}");
            Console.WriteLine();
        }

        static void TestNumberTheory()
        {
            Console.WriteLine("Testing Number Theory Functions...");
            int a = 12, b = 18;

            Console.WriteLine($"GCD({a}, {b}) = {MathUtils.NumberTheory.GCD(a, b)}");
            Console.WriteLine($"LCM({a}, {b}) = {MathUtils.NumberTheory.LCM(a, b)}");
            Console.WriteLine($"IsPrime(17) = {MathUtils.NumberTheory.IsPrime(17)}");
            Console.WriteLine($"Fibonacci(10) = {MathUtils.NumberTheory.Fibonacci(10)}");
            Console.WriteLine();
        }

        static void TestAlgorithms()
        {
            Console.WriteLine("Testing Sorting and Searching Algorithms...");
            int[] arr = { 5, 3, 8, 1, 2 };

            Console.WriteLine("Original Array: " + string.Join(", ", arr));
            MathUtils.Algorithms.BubbleSort(arr);
            Console.WriteLine("BubbleSort: " + string.Join(", ", arr));

            int[] arr2 = { 5, 3, 8, 1, 2 };
            MathUtils.Algorithms.QuickSort(arr2, 0, arr2.Length - 1);
            Console.WriteLine("QuickSort: " + string.Join(", ", arr2));

            Console.WriteLine($"LinearSearch(8) = {MathUtils.Algorithms.LinearSearch(arr, 8)}");
            Console.WriteLine($"BinarySearch(8) = {MathUtils.Algorithms.BinarySearch(arr, 8)}");
            Console.WriteLine();
        }

        static void TestMatrix()
        {
            Console.WriteLine("Testing Matrix Functions...");
            double[,] A = { { 1, 2 }, { 3, 4 } };
            double[,] B = { { 5, 6 }, { 7, 8 } };

            Console.WriteLine("Matrix A + B:");
            PrintMatrix(MathUtils.Matrix.Add(A, B));

            Console.WriteLine("Matrix A * B:");
            PrintMatrix(MathUtils.Matrix.Multiply(A, B));

            Console.WriteLine($"Determinant of A = {MathUtils.Matrix.Determinant(A)}");
            Console.WriteLine("Inverse of A:");
            PrintMatrix(MathUtils.Matrix.Inverse2x2(A));
            Console.WriteLine();
        }

        static void TestGeometry()
        {
            Console.WriteLine("Testing Geometry Functions...");
            double x1 = 0, y1 = 0, x2 = 3, y2 = 4;

            Console.WriteLine($"Distance2D({x1}, {y1}, {x2}, {y2}) = {MathUtils.Geometry.Distance2D(x1, y1, x2, y2)}");
            Console.WriteLine($"Distance3D({x1}, {y1}, 0, {x2}, {y2}, 0) = {MathUtils.Geometry.Distance3D(x1, y1, 0, x2, y2, 0)}");
            Console.WriteLine();
        }

        static void TestComplexNumbers()
        {
            Console.WriteLine("Testing Complex Numbers...");
            MathUtils.Complex c1 = new MathUtils.Complex(1, 2);
            MathUtils.Complex c2 = new MathUtils.Complex(3, 4);

            Console.WriteLine($"c1 + c2 = {c1 + c2}");
            Console.WriteLine($"c1 * c2 = {c1 * c2}");
            Console.WriteLine();
        }

        static void TestRandom()
        {
            Console.WriteLine("Testing Random Number Generation...");
            Console.WriteLine($"Random() = {MathUtils.Random()}");
            Console.WriteLine($"RandomInt(1, 100) = {MathUtils.RandomInt(1, 100)}");
            Console.WriteLine();
        }

        static void PrintMatrix(double[,] matrix)
        {
            for (int i = 0; i < matrix.GetLength(0); i++)
            {
                for (int j = 0; j < matrix.GetLength(1); j++)
                {
                    Console.Write(matrix[i, j] + "\t");
                }
                Console.WriteLine();
            }
        }
    }
}