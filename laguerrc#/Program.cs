using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace Laguerr
{
    public class Integral
    {
        private double _a;
        private double _b;
        private double _e;

        public Integral(double a, double b, double e)
        {
            this.A = a;
            this.B = b;
            this.E = e;
        }

        public double A
        {
            get { return _a; }
            set { _a = value; }
        }

        public double B
        {
            get { return _b; }
            set
            {
                if (value < this.A)
                    throw new ArgumentException("b must be greater than a");
                _b = value;
            }
        }

        public double E
        {
            get { return _e; }
            set
            {
                if (value < 0)
                    throw new ArgumentException("e must be greater than 0");
                _e = value;
            }
        }

        public double RectangleIntegral(Func<double, double> f, int steps = 1000)
        {
            double res1 = 0;
            double res2 = 0;

            for (int i = 0; i < steps; i++)
            {
                res1 += f(this.A + (this.B - this.A) / steps * i);
            }
            res1 *= (this.B - this.A) / steps;

            steps *= 2;

            for (int i = 0; i < steps; i++)
            {
                res2 += f(this.A + (this.B - this.A) / steps * i);
            }
            res2 *= (this.B - this.A) / steps;

            while (Math.Abs(res1 - res2) > this.E)
            {
                res1 = res2;
                steps *= 2;
                res2 = 0;

                for (int i = 0; i < steps; i++)
                {
                    res2 += f(this.A + (this.B - this.A) / steps * i);
                }
                res2 *= (this.B - this.A) / steps;
            }

            return Math.Round(res2, (int)Math.Log10(1 / this.E));
        }
    }

    public class Laguerre
    {
        private double _beta;
        private double _sigma;

        public double Beta
        {
            get { return _beta; }
            set
            {
                if (value < 0)
                    throw new ArgumentException("beta must be greater than 0");
                _beta = value;
            }
        }

        public double Sigma
        {
            get { return _sigma; }
            set
            {
                if (value < this.Beta)
                    throw new ArgumentException("sigma must be greater than beta");
                _sigma = value;
            }
        }

        public Laguerre(double beta, double sigma)
        {
            this.Beta = beta;
            this.Sigma = sigma;
        }

        public double LaguerreFunction(double t, int n)
        {
            double l0 = Math.Sqrt(this.Sigma) * Math.Exp(-this.Beta * t / 2);
            double l1 = Math.Sqrt(this.Sigma) * (1 - this.Sigma * t) * Math.Exp(-this.Beta * t / 2);

            if (n == 0)
                return l0;
            if (n == 1)
                return l1;
            if (n >= 2)
            {
                double lNext = ((2 * n - 1 - t * this.Sigma) / (double)n) * l1 - ((double)(n - 1) / n) * l0;
                for (int j = 3; j <= n; j++)
                {
                    double lTemp = l1;
                    l1 = lNext;
                    l0 = lTemp;
                    lNext = ((2 * j - 1 - t * this.Sigma) / j) * l1 - ((double)(j - 1) / j) * l0;
                }
                return lNext;
            }
            return 0;
        }

        public List<double> TabulateLaguerre(double T, int n, double step = 0.1)
        {
            List<double> values = Enumerable.Range(0, (int)(T / step)).Select(x => x * step).ToList();
            List<double> results = new List<double>();
            foreach (var i in values)
            {
                results.Add(LaguerreFunction(i, n));
            }

            return results;
        }

        /* public List<double> TransformLaguerre(Func<double, double> f, double T, int N)
         {
             List<double> ns = Enumerable.Range(0, N + 1).Select(x => (double)x).ToList();
             List<double> results = new List<double>();

             foreach (var i in ns)
             {
                 Func<double, double> func = (double x) => f(x) * LaguerreFunction(x, (int)i) * Math.Exp(-(this.Sigma - this.Beta) * x);
                 results.Add(new Integral(0, T, 0.0001).RectangleIntegral(func));
             }

             return results;
         }*/
        public List<double> TransformLaguerre(Func<double, double> f, double T, int N)
        {
            List<double> results = new List<double>();

            for (int i = 0; i <= N; i++)
            {
                Func<double, double> func = x => f(x) * LaguerreFunction(x, i) * Math.Exp(-(Sigma - Beta) * x);
                results.Add(new Integral(0, T, 0.0001).RectangleIntegral(func));
            }

            return results;
        }

        /*public double ReversedTransformLaguerre(List<double> seq, double t)
        {
            double sumRes = 0;

            for (int i = 0; i < seq.Count; i++)
            {
                sumRes += seq[i] * LaguerreFunction(t, i);
            }

            return sumRes;
        }*/
        public double ReversedTransformLaguerre(List<double> hList, double t)
        {
            double resultSum = 0;
            List<double> hListNew = hList.FindAll(x => x != 0);
            for (int i = 0; i < hListNew.Count; i++)
            {
                resultSum += hListNew[i] * LaguerreFunction(t, i);
            }

            return resultSum;
        }
        public List<double> ReversedTransformLaguerre(List<double> hList, double t1, double t2)
        {
            List<double> result = new List<double>() { };
            int size = 1000;
            double step = (t2 - t1) / size;
            for (int i = 0; i < size + 1; i++)
            {
                double t = t1 + i * step;

                result.Add(ReversedTransformLaguerre(hList, t));
            }

            return result;
        }
    }


    public class Experiment
    {
        private Laguerre _laguerre;

        public Experiment(Laguerre laguerre)
        {
            this.Laguerre = laguerre;
        }

        public Laguerre Laguerre
        {
            get { return _laguerre; }
            set
            {
                if (value == null)
                    throw new ArgumentNullException("laguerre cannot be null");
                _laguerre = value;
            }
        }

        public Tuple<List<double>, double> RunExperiment(double T, int N = 20, double eps = 0.001)
        {
            List<double> t = Enumerable.Range(0, 1001).Select(x => T * x / 1000).ToList();
            foreach (var i in t)
            {
                bool check = true;
                for (int n = 0; n <= N; n++)
                {
                    if (Math.Abs(Laguerre.LaguerreFunction(i, n)) >= eps)
                    {
                        check = false;
                        break;
                    }
                }
                if (check)
                {
                    List<double> ns = Enumerable.Range(0, N + 1).Select(x => (double)x).ToList();
                    return Tuple.Create(ns, i);
                }
            }
            return null;
        }
    }

    public class Functions
    {
        public static double F(double t)
        {
            if (t >= 0 && t <= 2 * Math.PI)
            {
                return Math.Sin(t - Math.PI / 2) + 1;
            }
            else
            {
                return 0;
            }
        }

        public static Func<double, double> Gauss(double mu, double lambda)
        {
            if (lambda <= 0)
                throw new ArgumentException("lambda must be greater than 0");

            return t =>
            {
                double exp = -Math.Pow((t - mu), 2) / (2 * Math.Pow(lambda, 2));
                double denom = lambda * Math.Sqrt(2 * Math.PI);

                return Math.Exp(exp) / denom;
            };
        }
        class Program
        {

            static void Main(string[] args)
            {
                double beta = 2;
                double sigma = 4;
                double T = 100;
                int N = 20;
                double eps = 0.001;

                Laguerre laguerre = new Laguerre(beta, sigma);
                Experiment exp = new Experiment(laguerre);

                var result = exp.RunExperiment(T, N, eps);

                if (result != null)
                {
                    List<double> ns = result.Item1;
                    double ans = result.Item2;

                    var transformed = laguerre.TransformLaguerre(Functions.F, ans, N);
                    WriteToCsv(transformed, "transformed_laguerre.csv");

                    List<double> tabulated = laguerre.TabulateLaguerre(ans, 5);
                    WriteToCsv(tabulated, "tabulated_laguerre.csv");

                    // List<double> reversed = new List<double>();
                    int t = 2;
                   
                    var rtransformed = laguerre.ReversedTransformLaguerre(transformed, t);
                    WriteToCsv(new List<double> { rtransformed }, "rtransformed_laguerre.csv");


                    /*List<double> reversed1 = new List<double>();

                    for (int i = 0; i < 20; i++)
                    {
                        reversed1.Add(i);
                    }

                    double t = 2;
                    List<double> reversedTransform = laguerre.ReversedTransformLaguerrens(, N);
                    WriteToCsv(reversedTransform, "reversedlaguerre.csv");*/

                    Func<double, double> gaussFunc = Functions.Gauss(7, 9);
                    List<double> gaussTabulated = laguerre.TransformLaguerre(gaussFunc, T, N);
                    WriteToCsv(gaussTabulated, "gauss_tabulated.csv");

                    Func<double, double> gaussFunc3 = Functions.Gauss(3, 3);
                    List<double> gaussTabulated3 = laguerre.TransformLaguerre(gaussFunc3, T, N);
                    WriteToCsv(gaussTabulated3, "gauss_tabulated_3_lambda.csv");

                    Func<double, double> gaussFunc5 = Functions.Gauss(6, 5);
                    List<double> gaussTabulated5 = laguerre.TransformLaguerre(gaussFunc5, T, N);
                    WriteToCsv(gaussTabulated5, "gauss_tabulated_5_lambda.csv");


                    //
                    Console.WriteLine("\nExperiment mu = 3*lambda");
                    double lambda_1 = 5;
                    var Gaus_transform_tabulate_exp_3 = laguerre.TransformLaguerre(Functions.Gauss(3 * lambda_1, lambda_1), T, N);
                    WriteToCsv(Gaus_transform_tabulate_exp_3, "transformedgausstabulated_3.csv");

                    Console.WriteLine("\nExperiment mu = 6*lambda");
                    double lambda_2 = 3;
                    var Gaus_transform_tabulate_exp_5 = laguerre.TransformLaguerre(Functions.Gauss(3 * lambda_2, lambda_2), T, N);
                    WriteToCsv(Gaus_transform_tabulate_exp_5, "transformedgausstabulated_5.csv");


                }
                else
                {
                    Console.WriteLine("Experiment did not converge.");
                }
            }

            static void WriteToCsv(IEnumerable<double> data, string filePath)
            {
                using (var writer = new StreamWriter(filePath))
                {
                    foreach (var item in data)
                    {
                        writer.WriteLine(item);
                    }
                }
            }
        }
    }
}
