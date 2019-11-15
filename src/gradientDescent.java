import java.util.ArrayList;

public class gradientDescent {
	
	private double alpha;
	private int maxIter;
	private double stoppingParameter;
	private int candidates;
	private int dimensionality;
	private int[] charWeights;
	private double penaltyTerm;
	private ArrayList<Double> charMatrix = new ArrayList<>();	//[rowLen][colLen]
	private ArrayList<Integer> ltrIds = new ArrayList<>();
	double[] C;
	int N;
	
	public gradientDescent(double _alpha, int _maxIter, double _stoppingParameter, int[] _charWeights, double _penalty) {
		this.alpha = _alpha;
		this.maxIter = _maxIter;
		this.stoppingParameter = _stoppingParameter;
		this.dimensionality = _charWeights.length;
		this.charWeights = _charWeights;
		this.penaltyTerm = _penalty;
	}

	public void addCandidate(double[] characteristics, Integer id) {
		if(characteristics.length != dimensionality) {
			throw new java.lang.RuntimeException("Dimensionality mismatch. Check number of characteristics.");
		}
		ltrIds.add(id);
		for(double c: characteristics) {
			charMatrix.add(c);
		}
		
	}
	
	private void calculateC() {
		this.C = new double[N];
		int offset = 0;
		for(int i = 0; i < C.length; i++) {
			double c = 0;
			for(int j = 0; j < charWeights.length; j++) {
				c += charWeights[j] * charMatrix.get(j + offset);
			}
			C[i] = c;
			offset += 3; 	//skip ahead to next row
		}
	}
	
	private double[] calculateB_i(double[] B) {
		double[] gradient = calculateGradient(B);
		double[] descent = new double[N];
		for(int i = 0; i < N; i++) {
			descent[i] = B[i] - (alpha * gradient[i]); 
		}	
		return descent;
	}
	
	private double calculateCost(double[] B) {
		return calculateCB(B) - calculatePenalty(B);
	}
	
	private double calculateCB(double[] B) {
		double CB = 0;
		if(B.length != N) {
			throw new java.lang.RuntimeException("Dimensionality mismatch.");
		}
		for(int i = 0; i < B.length; i++) {
			CB += B[i] * C[i];
		}
		return CB;
	}
	
	private double calculatePenalty(double[] B) {
		if((B.length != N)) {
			throw new java.lang.RuntimeException("Dimensionality mismatch.");
		}
		double penalty = 0;
		for(double b: B) {
			penalty += b;
		}
		return penalty / N;
	}
	
	private double[] calculateGradient(double[] B) {
		if(B.length != N) {
			throw new java.lang.RuntimeException("Dimensionality mismatch.");
		}
		double[] gradient = new double[N];
		for(int i = 0; i < C.length; i++) {
			gradient[i] = C[i] - (B[i] / N);
		}
		return gradient;
	}
	
	private double[] initB() {
		double[] initialB = new double[N];
		for(int i = 0; i < N; i++) {
			initialB[i] = Math.random(); 
		}
		return initialB;
	}
	
	private double[] normalizeB(double[] B) {
		double[] normalizedB = new double[N];
		double max = 0;
		for(double b:B) {
			if(b>=max) {
				max = b;
			}
		}
		double min = max;
		for(double b:B) {
			if(b<=min) {
				min = b;
			}
		}
		for(int i = 0; i < N; i++) {
			normalizedB[i] = (B[i] - min) / (max - min);
		}
		return normalizedB;
	}
	
	public ArrayList<Integer> runGradientDescent() {
		this.N = ltrIds.size();
		calculateC();
		ArrayList<Integer> chosenLTRs = new ArrayList<>();
		double[] B = initB();
		double prevCost = calculateCost(B) + stoppingParameter + 1;
		int iterations = 0;
		while((prevCost - calculateCost(B)) > stoppingParameter) {
			if(iterations >= maxIter) {
				break;
			}
			double[] B_i = calculateB_i(B);
			prevCost = calculateCost(B);
			if(prevCost > calculateCost(B_i)) {
				System.out.print("Warning, cost function is not decreasing.");
			}
			B = B_i;
			iterations++;
		}
		
		double[] normB = normalizeB(B);
		for(int i = 0; i < N; i++) {
			if(normB[i] >= 0.5) {
				chosenLTRs.add(ltrIds.get(i));
			}
		}
		
		return chosenLTRs;
	}
}
