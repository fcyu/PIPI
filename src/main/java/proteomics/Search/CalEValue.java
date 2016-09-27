package proteomics.Search;

import proteomics.PIPI;
import proteomics.Types.FinalResultEntry;

public class CalEValue {

    public CalEValue(FinalResultEntry resultEntry) {
        int[] scoreHistogram = resultEntry.getScoreHistogram();
        double histogramBinSize = resultEntry.getHistogramBinSize();
        double histogramBinOffset = resultEntry.getHistogramBinOffset();

        // find the max non-zero idx. It may be an outlier.
        int maxNonzeroIdx = 0;
        for (int i = 0; i < scoreHistogram.length; ++i) {
            if (scoreHistogram[i] > 0) {
                maxNonzeroIdx = i;
            }
        }

        // generate survival array
        int[] survivalCountArray = new int[maxNonzeroIdx + 1];
        survivalCountArray[maxNonzeroIdx] = scoreHistogram[maxNonzeroIdx];
        for (int i = maxNonzeroIdx - 1; i >= 0; --i) {
            survivalCountArray[i] = survivalCountArray[i + 1] + scoreHistogram[i];
        }

        // find the next max nonzero idx. from this down to the changing point are from null.
        int topCount = survivalCountArray[maxNonzeroIdx] + 9; // these are from non-null.
        int nullEndIdx = 0;
        for (int i = maxNonzeroIdx - 1; i >= 0; --i) {
            if (survivalCountArray[i] > topCount) {
                nullEndIdx = i;
                break;
            }
        }

        // log transform. eliminate the non-null points
        double[] lnSurvivalCountArray = new double[nullEndIdx + 1];
        for (int i = 0; i <= nullEndIdx; ++i) {
            lnSurvivalCountArray[i] = Math.log(survivalCountArray[i]);
        }

        // try and find the best linear regression result
        int startIdx = (int) (nullEndIdx * 0.5);
        double maxRSquare = 0;
        double optimalSlope = 1;
        double optimalIntercept = 0;
        int optimalStartIdx = 0;
        while (startIdx >= 0) {
            double xSum = 0;
            double ySum = 0;
            double xMean = 0;
            double yMean = 0;
            double xxMean = 0;
            double yyMean = 0;
            double xyMean = 0;
            int pointNum = 0;
            double rSquare;
            double slope;
            double intercept;

            for (int i = startIdx; i <= nullEndIdx; i++) {
                if (survivalCountArray[i] > 0) {
                    xSum += i;
                    ySum += lnSurvivalCountArray[i];
                    ++pointNum;
                }
            }

            if (pointNum > 0) {
                xMean = xSum / pointNum;
                yMean = ySum / pointNum;
            }

            for (int i = startIdx; i <= nullEndIdx; i++) {
                if (survivalCountArray[i] > 0) {
                    double dX = i - xMean;
                    double dY = lnSurvivalCountArray[i] - yMean;
                    xxMean += dX * dX;
                    yyMean += dY * dY;
                    xyMean += dX * dY;
                }
            }

            if (xxMean > 0) {
                slope = xyMean / xxMean;
            } else {
                slope = 0;
            }

            intercept = yMean - slope * xMean;

            double ssr = 0.0;
            for (int i = startIdx; i <= nullEndIdx; i++) {
                if (survivalCountArray[i] > 0) {
                    double fit = slope * i + intercept;
                    ssr += (fit - yMean) * (fit - yMean);
                }
            }

            rSquare = ssr / yyMean;
            if (rSquare > maxRSquare) {
                optimalSlope = slope;
                optimalIntercept = intercept;
                maxRSquare = rSquare;
                optimalStartIdx = startIdx;
            }
            --startIdx;
        }

        if (PIPI.DEV) {
            resultEntry.setLnSurvivalArray(lnSurvivalCountArray);
            resultEntry.setStartIdx(optimalStartIdx);
            resultEntry.setRSquare(maxRSquare);
        }

        if (optimalSlope >= 0) {
            resultEntry.setEValue(9999);
        } else {
            double eValue = Math.exp(optimalSlope * (resultEntry.getScore() / histogramBinSize + histogramBinOffset) + optimalIntercept);
            resultEntry.setEValue(eValue);
        }
    }
}
