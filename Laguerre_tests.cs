using Laguerr;
using Xunit;
using static Laguerr.Laguerre;

public class TestData : IDisposable
{
    public Integral IntegralInstance { get; set; }
    public Experiment ExperimentInstance { get; set; }
    public Laguerre LaguerreInstance { get; set; }

    public TestData()
    {
        IntegralInstance = new Integral(1, 2, 0.1);
        LaguerreInstance = new Laguerre(2, 4);
        ExperimentInstance = new Experiment(LaguerreInstance);
    }

    public void Dispose()
    {
        IntegralInstance = null;
        LaguerreInstance = null;
        ExperimentInstance = null;
    }
}

public class LaguerreTests : IClassFixture<TestData>
{
    private readonly TestData _fixture;

    public LaguerreTests(TestData fixture)
    {
        _fixture = fixture;
    }

    [Theory]
    [InlineData(2, 4)]
    public void LaguerreFunctionCorrectResult(double beta, double sigma)
    {
        Laguerre laguerre = new Laguerre(beta, sigma);

        Assert.Equal(2, 2, Math.Round(laguerre.LaguerreFunction(0, 0), 1));
        Assert.Equal(1.06, Math.Round(laguerre.LaguerreFunction(1, 3), 2));
    }


}

public class ExperimentTests : IClassFixture<TestData>
{
    private readonly TestData fixture;

    public ExperimentTests(TestData fixture)
    {
        this.fixture = fixture;
    }

    [Fact]
    public void Experiment_Returns_Result()
    {
        Experiment exp = fixture.ExperimentInstance;

        var result = exp.RunExperiment(100);

        Assert.NotNull(result);
        Assert.IsType<Tuple<List<double>, double>>(result);
    }
}

public class IntegralTests : IClassFixture<TestData>
{
    private readonly TestData fixture;

    public IntegralTests(TestData fixture)
    {
        this.fixture = fixture;
    }

    [Fact]
    public void RectangleIntegralCorrectResult()
    {
        Integral integral = fixture.IntegralInstance;

        double result = integral.RectangleIntegral(x => x * x);

        Assert.Equal(2.3, result);
    }
}