package proteomics.Score;

public class Binomial {

    private static final int maxValue = 1000;
    private final double[] log10FactorialArray = new double[maxValue];

    public Binomial() {
        double temp = 0;
        log10FactorialArray[0] = 0;
        for (int i = 1; i < maxValue; ++i) {
            temp += Math.log10(i);
            log10FactorialArray[i] = temp;
        }
    }

    public double calPValue(int N, int K, double p) {
        double pValue = 0;
        for (int i = K; i <= N; ++i) {
            pValue += Math.pow(10, calLog10PMF(N, i, p));
        }
        return pValue;
    }

    private double calLog10PMF(int N, int i, double p) {
        return log10FactorialArray[N] - log10FactorialArray[i] - log10FactorialArray[N - i] + i * Math.log10(p) + (N - i) * Math.log10(1 - p);
    }
}
