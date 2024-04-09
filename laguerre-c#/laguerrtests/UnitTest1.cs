using System;
using System.Collections.Generic;
using System.Reflection.Metadata.Ecma335;
using Xunit;

class Integral
{
    //private double _a;
    //private double _b;
    //private double _e;

    public Integral(double a, double b, double e)
    {
        this.A = a;
        this.B = b;
        this.E = e;
    }

    public double A
    {
        get { return A; }
        set { A = value; }
    }

    public double B
    {
        get { return B; }
        set
        {
            if (value < this.A)
                throw new ArgumentException("b must be greater than a");
            B = value;
        }
    }

    public double E
    {
        get { return E; }
        set
        {
            if (value < 0)
                throw new ArgumentException("e must be greater than 0");
            E = value;
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


class Laguerre
{
    //private double _beta;
   // private double _sigma;
   public double Beta
    {
        get=>this.Beta;
        set {if  (value < 0)
            {
                throw new ArgumentException("beta must be greater than 0");
            }
        this.Beta = value;

        }

    }
   public double Sigma
    {
        get => this.Sigma;
        set { if (value < this.Beta)
            {
                throw new ArgumentException("sigma must be greater than beta");

            }
            this.Sigma = value;
        }
    }

    public Laguerre(double beta, double sigma)
    {
        this.Beta = beta;
        this.Sigma = sigma;
    }

    //public double Beta
    //{
    //    get { return _beta; }
    //    set
    //    {
    //        if (value < 0)
    //            throw new ArgumentException("beta must be greater than 0");
    //        _beta = value;
    //    }
    //}

    //public double Sigma
    //{
    //    get { return _sigma; }
    //    set
    //    {
    //        if (value < this.Beta)
    //            throw new ArgumentException("sigma must be greater than beta");
    //        _sigma = value;
    //    }
    //}

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

    public List<double> TransformLaguerre(Func<double, double> f, double T, int N)
    {
        List<double> ns = Enumerable.Range(0, N + 1).Select(x => (double)x).ToList();
        List<double> results = new List<double>();

        foreach (var i in ns)
        {
            Func<double, double> func = (double x) => f(x) * LaguerreFunction(x, (int)i) * Math.Exp(-(this.Sigma - this.Beta) * x);
            results.Add(new Integral(0, T, 0.0001).RectangleIntegral(func));
        }

        return results;
    }

    public double ReversedTransformLaguerre(List<double> seq, double t)
    {
        double sumRes = 0;

        for (int i = 0; i < seq.Count; i++)
        {
            sumRes += seq[i] * LaguerreFunction(t, i);
        }

        return sumRes;
    }
}

class Experiment
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


public class Program
{
    static void Main(string[] args)
    {
        Laguerre laguerre = new Laguerre(2, 4);

        Experiment exp = new Experiment(laguerre);

        var result = exp.RunExperiment(100);

        if (result != null)
        {
            List<double> ns = result.Item1;
            double ans = result.Item2;

            Func<double, double> f = x =>
            {
                if (x >= 2 * Math.PI)
                    return 0;
                else
                    return Math.Sin(x - Math.PI / 2) + 1;
            };

            var transformed = laguerre.TransformLaguerre(f, ans, 20);

            Console.WriteLine("Transformed Laguerre data:");
            foreach (var item in transformed)
            {
                Console.WriteLine($"N: {item}");
            }
        }
    }
}

public class LaguerreTests
{
    [Theory]
    [InlineData(2, 4)]
    [InlineData(3, 5)]
    public void LaguerreFunction_Returns_CorrectResult(double beta, double sigma)
    {
        Laguerre laguerre = new Laguerre(beta, sigma);

        Assert.Equal(1.0, laguerre.LaguerreFunction(0, 0));
        Assert.Equal(1.0, laguerre.LaguerreFunction(0, 1));
        Assert.Equal(0.0, laguerre.LaguerreFunction(0, 2));
    }

    [Fact]
    public void TabulateLaguerre_Returns_CorrectResult()
    {
        Laguerre laguerre = new Laguerre(2, 4);

        List<double> tabulatedValues = laguerre.TabulateLaguerre(10, 3, 0.1);

        Assert.NotNull(tabulatedValues);
        Assert.Equal(101, tabulatedValues.Count);
        Assert.Equal(1.0, tabulatedValues[0]);
        Assert.Equal(0.1, tabulatedValues[1]);
    }

    [Fact]
    public void TransformLaguerre_Returns_CorrectResult()
    {
        Laguerre laguerre = new Laguerre(2, 4);
        Func<double, double> f = x => Math.Sin(x);

        List<double> transformedValues = laguerre.TransformLaguerre(f, 10, 3);

        Assert.NotNull(transformedValues);
        Assert.Equal(4, transformedValues.Count);
    }
}

public class ExperimentTests
{
    [Fact]
    public void Experiment_Returns_Result()
    {
        Laguerre laguerre = new Laguerre(2, 4);
        Experiment exp = new Experiment(laguerre);

        var result = exp.RunExperiment(100);

        Assert.NotNull(result);
        Assert.IsType<Tuple<List<double>, double>>(result);
    }
}

