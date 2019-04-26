
import java.util.HashMap;
import java.util.Map;
import java.util.Iterator;
import java.lang.Math;

import java.io.IOException;

public class StatsCalculator {

		private static boolean DEBUG = true;

		public static void anova(Double[][] samples) {
			Integer k = samples.length;
			Integer n = getN(samples);
			Double[] means = getMeans(samples, k);
			Double[] stdDevs = getStdDevs(k, n, samples, means);
			Double combinedMean = getCombinedMean(samples, k, n);
			Double FStat = getFStat(samples, means, k, n, stdDevs, combinedMean);
			
			if(DEBUG) {
				printChecksum(k, n, combinedMean, FStat);
			}
		}

		public static Integer getN(Double[][] samples) {
			Integer N = 0;
			for(int x = 0; x < samples.length; x++) {
				N += samples[x].length;
			}
			return N;
		}

		public static Double[] getMeans(Double[][] samples, Integer k) {
			Double[] means = new Double[k];
			for(int x = 0; x < k; x++) {
				Double thisMean = 0.0;
				int setSize = samples[x].length;
				for(int y = 0; y < setSize; y++) {
					thisMean += samples[x][y];
				}
				thisMean /= setSize;
				means[x] = thisMean;
			}
			return means;
		}

		public static Double getCombinedMean(Double[][] samples, Integer k, Integer n) {
			Double combinedMean = 0.0;
			for(int x = 0; x < k; x++) {
				for(int y = 0; y < samples[x].length; y++) {
					combinedMean += samples[x][y];
				}
			}
			return (combinedMean / n);
		}

		public static Double[] getStdDevs(Integer k, Integer n, Double[][] samples, Double[] means) {
			//Calculate the standard deviation for each sample group
			Double[] stdDevs = new Double[k];
			for(int i = 0; i < k; i++) {
				if(samples[i].length != 0) {
					Double stdDev = 0.0;
					for(int j = 0; j < samples[i].length; j++) {
						Double val = samples[i][j];
						if(val != null && means[i] != null) {
							stdDev += Math.pow((val - means[i]), 2);
						}
					}
					stdDev /= (n - 1);
					if(Double.compare(stdDev, 0.0) < 0) {
						stdDev *= 1;
					}
					stdDev = Math.sqrt(stdDev);
					stdDevs[i] = stdDev;
				}
			}
			return stdDevs;
		}


		public static Double getFStat(Double[][] samples, Double[] means, Integer k, Integer n, Double[] stdDevs, Double combinedMean) {
			//Performs ANOVA statistical analysis given a set of means, combined sample size, # of test groups, combined mean of all sets, and the raw set of values
			//The calculations performed here are taken straight from my Statistics textbook


			Double SSG = 0.0;
			//SSG = SUM(Nsub(j)(means[j] - combinedMean)^2)
			for(int j = 0; j < k; j++) {
				if(means[j] != null) {
					int NsubJ = samples[j].length;
					SSG += NsubJ * Math.pow((means[j] - combinedMean), 2);
				}
			}
			//SSE = SUM( (sigByVal[j] - 1) * stdDevs[j]^2
			Double SSE = 0.0;
			for(int j = 0; j < k; j++) {
				if(stdDevs[j] != null) {
					int NsubJ = samples[j].length;
					SSE += (NsubJ - 1) * (stdDevs[j] * stdDevs[j]);
				}
			}

			Double MSG = SSG / (k-1);
			Double MSE = SSE / (n - k);


			if(DEBUG) {
				printAnovaIntermediates(SSG, SSE, MSG, MSE);
			}

			Double F = MSG / MSE;
			return F;
		}

		public static void printAnovaIntermediates(Double SSG, Double SSE, Double MSG, Double MSE) {
			System.out.println("Sum of Squares (Group): " + SSG.toString());
			System.out.println("Sum of Squares (Error): " + SSE.toString());
			System.out.println("Mean Square (Group): " + MSG.toString());
			System.out.println("Mean Square (Error): " + MSE.toString());
		}

		public static void printChecksum(Integer k, Integer n, Double combinedMean, Double FStat) {
			System.out.println("k (number of sample groups): " + k.toString());
			System.out.println("n (total number of samples): " + n.toString());
			System.out.println("Combined Mean: " + combinedMean.toString());
			System.out.println("F Statistic: " + FStat.toString());
		}	
}
