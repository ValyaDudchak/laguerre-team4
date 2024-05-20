using System;
using System.Collections;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Runtime.Intrinsics.X86;

namespace Laguerr
{
    public class Integral
    {
        private double _a;
        private double _b;
        private double _e;

        public Integral(double a, double b, double e)
        {
            A = a;
            B = b;
            E = e;
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

        public double RectangleIntegral(Func<double, double> f,double a,double b, int steps = 1000)
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
        private double? experiment;

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
            this.experiment = null;

        }

        public double LaguerreFunction(double t, int n)
        {

            double l0 = Math.Sqrt(this.Sigma) * Math.Exp(-this.Beta * t / 2);
            double l1 = Math.Sqrt(this.Sigma) * (1 - this.Sigma * t) * Math.Exp(-this.Beta * t / 2);

            if (n == 0)
                return l0;
            if (n == 1)
                return l1;
           
            for (int i = 2; i <= n; i++)
            {
                double temp = l1;
                l1 = ((2 * i - 1 - Sigma * t) * l1 / i) - ((i - 1) * l0 / i);
                l0 = temp;
            }

            return l1;
        }

        public List<Tuple<double, double>> TabulateLaguerre(double T, int n, int s = 100)
        {
            List<Tuple<double, double>> values = new List<Tuple<double, double>>();
            double[] t = new double[s];
            for (int i = 0; i < s; i++)
            {
                t[i] = i * T / s;
                double l = LaguerreFunction(t[i], n);
                values.Add(Tuple.Create(t[i], l));
            }

            return values;
        }

        public double LaguerreTransform(Func<double, double> f, int n)
        {
            Func<double, double> integrand = t => f(t) * LaguerreFunction(t, n) * Math.Exp(-t * (Sigma - Beta));
            double b = Experiment(100).Item1 ?? double.PositiveInfinity;
            return new Integral(0, b, 1e-5).RectangleIntegral(integrand,0,b,1000);
        }

        public List<double> TabulateTransform(Func<double, double> f, int n)
        {
            List<double> res = new List<double>();
            for (int i = 0; i < n; i++)
            {
                double transform = LaguerreTransform(f, n);
                res.Add(transform);
            }
            return res;
        }

        public double ReversedTransformLaguerre(List<double> seq, double t)
        {
            double sumRes = 0;
            List<double> hListNew = seq.FindAll(x => x != 0);

            for (int i = 0; i < seq.Count; i++)
            {
                sumRes += seq[i] * LaguerreFunction(t, i);
            }

            return sumRes;
        }
     

        public List<double> ReversedLaguerreTransformationTabulate(List<double> seq, double t1, double t2,double step=0.1)
        {
            List<double> result = new List<double>() { };
            //int size = 1000;
            int size = (int)((t2 - t1) / step);

            //double step = (t2 - t1) / size;
            for (int i = 0; i < size + 1; i++)
            {
                double t = t1 + i * step;

                result.Add(ReversedTransformLaguerre(seq, t));
            }

            return result;
        }


        public Tuple<double?, List<Tuple<double, List<double>>>> Experiment(double T, double epsilon = 1e-3, int N = 20)
        {
            List<double> t = new List<double>();
            for (int i = 0; i < 1000; i++)
                t.Add(i * T / 1000);

            List<Tuple<double, List<double>>> experimentValues = new List<Tuple<double, List<double>>>();
            double? result = null;
            foreach (double i in t)
            {
                bool flag = true;
                List<double> laguerreValues = new List<double>();
                for (int j = 1; j <= N; j++)
                {
                    double laguerreVal = LaguerreFunction(i, j);
                    laguerreValues.Add(laguerreVal);
                    if (Math.Abs(laguerreVal) > epsilon)
                    {
                        flag = false;
                        break;
                    }
                }

                if (flag && !result.HasValue)
                {
                    result = i;
                }

                experimentValues.Add(Tuple.Create(i, laguerreValues));
            }

            experiment = result;

            return Tuple.Create(result, experimentValues);
        }


    }


    public class Function
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
        public static double F_n(double t)
        {
            if (t <= 2 * Math.PI && t >= 0) { return Math.Cos(t + (3 * Math.PI) / 2); }

            else { return 0; }
        }

        public static double F_v(double t)
        {
            return Math.Pow(t, 3) * Math.Cos(t);
        }

        public static double F_ma(double t)
        {
            return Math.Sin(2 * t) + 3 * t;
        }

        public static double F_r(double t)
        {
            return Math.Pow(t, 2) * Math.Sin(t);
        }
        public static Func<double, double> Gauss(double mu, double lambda)
        {
            if (lambda <= 0)
                throw new ArgumentException("lambda must be greater than 0");

            return t =>
            {
                double exponent = -Math.Pow((t - mu), 2) / (2 * Math.Pow(lambda, 2));
                double denominator = lambda * Math.Sqrt(2 * Math.PI);

                return Math.Exp(exponent) / denominator;
            };
        }
    }

    class Program
    {
        static void Main(string[] args)
        {


            var laguer = new Laguerre(2,4);

            // Laguerre tab 
            var tabulateLaguerre = laguer.TabulateLaguerre(10, 2);
            //foreach (var tab in tabulateLaguerre)
            //{
            //    Console.WriteLine(tab);
            //}
            writeInFiles(tabulateLaguerre, @"LaguerreTab.csv");

            // Experiment
            var experimentResult = laguer.Experiment(100);
            //if (experimentResult != null)
            //{
            //    Console.WriteLine($"Experiment result: {string.Join(", ", experimentResult.Item1)}, Time: {experimentResult.Item2}");
            //}
            //else
            //{
            //    Console.WriteLine("Experiment result: null");
            //}

            // Transformed Laguerre
            var transformLaguerre = laguer.TabulateTransform(Function.F,20);
            //foreach (var trans in transformLaguerre)
            //{
            //    Console.WriteLine(trans);
            //}
            writeInFiles(transformLaguerre, @"transformTab.csv");

            // Reversed Transform Laguerre
            var reversedTransformLaguerre = laguer.ReversedTransformLaguerre(transformLaguerre, Math.PI / 2);
           // Console.WriteLine($"Reversed Transform: {reversedTransformLaguerre}");
            var reversedTransformLaguerreTab = laguer.ReversedLaguerreTransformationTabulate(transformLaguerre, 0,Math.PI/2);
            writeInFiles(reversedTransformLaguerreTab, @"ReverseTransfTab.csv");


            // Transform Gauss
            var transformGauss = laguer.TabulateTransform(Function.Gauss(3, 5), 20);
            //foreach (var trans in transformGauss)
            //{
            //    Console.WriteLine(trans);
            //}
            writeInFiles(transformGauss, @"transformGauss.csv");


            //Gauss
            var reversedTransformGauss = laguer.ReversedTransformLaguerre(transformGauss, Math.PI);
            // Console.WriteLine($"Reversed Transform Gauss: {reversedTransformGauss}");
            var reversedTransformGaussTab = laguer.ReversedLaguerreTransformationTabulate(transformGauss, 0, Math.PI / 2);
            writeInFiles(reversedTransformGaussTab.Select(t => t.ToString()).ToList(), @"ReversedGaussTab.csv");


            //gauss for mu=5
            var transformGauss5 = laguer.TabulateTransform(Function.Gauss(5, 8), 20);
            //foreach (var trans in transformGauss)
            //{
            //    Console.WriteLine(trans);
            //}
            writeInFiles(transformGauss5, @"transformGauss5.csv");


            //Gauss
            var reversedTransformGauss5 = laguer.ReversedTransformLaguerre(transformGauss5, Math.PI);
            // Console.WriteLine($"Reversed Transform Gauss: {reversedTransformGauss}");
            var reversedTransformGaussTab5 = laguer.ReversedLaguerreTransformationTabulate(transformGauss5, 0, Math.PI / 2);
            writeInFiles(reversedTransformGaussTab5.Select(t => t.ToString()).ToList(), @"ReversedGaussTab5.csv");


            var valuesToWrite = new Dictionary<string, double>
        {
            //{ "experimentResult", experimentResult?.Item2 ?? 0 },
            { "reversedTransformLaguerre", reversedTransformLaguerre },
            { "reversedTransformGauss", reversedTransformGauss }
        };
            writeValuesToFile(valuesToWrite, @"otherValues.csv");


            static void writeInFiles<T>(List<T> data, string path)
            {
                using (StreamWriter writetext = new StreamWriter(path))
                {
                    foreach (T value in data)
                    {
                        writetext.WriteLine(value);
                    }
                }
            }

            static void writeValuesToFile(Dictionary<string, double> values, string path)
            {
                using (StreamWriter writetext = new StreamWriter(path))
                {
                    foreach (KeyValuePair<string, double> entry in values)
                    {
                        writetext.WriteLine($"{entry.Key}, {entry.Value}");
                    }
                }
            }


        }
    }

}
