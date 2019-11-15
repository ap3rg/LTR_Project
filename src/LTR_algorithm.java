import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Set;

import ngsep.gbs.VCFTranslator;
import ngsep.main.CommandsDescriptor;
import ngsep.sequences.DNASequence;
import ngsep.sequences.FMIndexSingleSequence;
import ngsep.sequences.QualifiedSequence;
import ngsep.sequences.QualifiedSequenceList;
import ngsep.sequences.io.FastaSequencesHandler;

public class LTR_algorithm {

	private int DMIN = 100;
	private int DMAX = 15000;
	private int LMIN = 250;
	private int LMAX = 600;
	private int LEX = 20;
	private int idCount = 1;
	String genomeFile;
	double mutationRate;
	private String outPrefix="./output";
	private FMIndexSingleSequence fmIndex;
	String genome = "";
	private static ArrayList<LTR> candidatePairs = new ArrayList<>();
	
	
	public static void main(String[] args) throws Exception {
		
		LTR_algorithm instance = new LTR_algorithm();
		int i = CommandsDescriptor.getInstance().loadOptions(instance, args);
		instance.genomeFile = args[i++];
		instance.outPrefix = args[i++];
		instance.mutationRate = Double.valueOf(args[i++]);
//		instance.DMIN = Integer.valueOf(args[i++]);
//		instance.DMAX = Integer.valueOf(args[i++]);
		instance.run();
	}

	private void run() throws IOException {
		this.fmIndex = preprocessGenome();
		System.out.println(System.currentTimeMillis());
		findCandidatePairs();
		System.out.println(System.currentTimeMillis());
		processCandidatePairs();
		printCandidatePositions();
		ArrayList<Integer> chosenLTRs = runModel();
		printFoundLTRs(chosenLTRs);
	}
	
	private void printFoundLTRs(ArrayList<Integer> chosenLTRs) throws FileNotFoundException {
		// TODO Auto-generated method stub
		try(PrintStream results = new PrintStream(outPrefix+"_results.txt");) {
			results.println("start_5\tend_5\tstart_3\tend_3\tLTR_length\tTotal_length\tpol\tgag\tTG\tCA");
			for(LTR ltr: candidatePairs) {
				if(Arrays.asList(chosenLTRs).contains(ltr.getID())) {
					results.println(
							ltr.getStart5() + "\t" +
							ltr.getEnd5() + "\t" +
							ltr.getStart3() + "\t" +
							ltr.getEnd3() + "\t" +
							ltr.getLTRLength() + "\t" +
							ltr.getTotalLength() + "\t" +
							ltr.getPol() + "\t" + 
							ltr.getGag() + "\t" + 
							ltr.getTG() + "\t" + 
							ltr.getCA() + "\t"
							);
				}
			}
		}
		
	}
	
	//Arrays.asList(yourArray).contains(yourValue)

	private ArrayList<Integer> runModel() {
		// TODO Auto-generated method stub
		double alpha = 0.005;
		int maxIter = 100000;
		double stoppingParameter = 0.1;
		int[] charWeights = {10, 5, 5, 5, 5};
		double penalty = 50;

		gradientDescent model = new gradientDescent(alpha, maxIter, stoppingParameter, charWeights, penalty);  
		for(LTR ltr: candidatePairs) {
			model.addCandidate(ltr.getCharVector(), ltr.getID());
		}
		 return model.runGradientDescent();
	}

	private void processCandidatePairs() {
		for(LTR ltr: candidatePairs) {
			setLTRLength(ltr);
			findPol(ltr);
		}
		
	}
	
	private void findPol(LTR ltr) {
		// TODO Auto-generated method stub
		
	}

	private void setLTRLength(LTR ltr) {
		int positionToCheck5 = ltr.getStart5() + LEX + 1;
		int positionToCheck3 = ltr.getStart3() + LEX + 1;
		int lengthLTR = LMAX;
		boolean mismatch = false;
		for(int i = 0; i < LMAX - LEX; i++) {
			//if two mismatches are found set length to current pos and stop.
			if((genome.charAt(positionToCheck5 + i) != genome.charAt(positionToCheck3 + i)) && mismatch) {
				lengthLTR = LEX + i;
				break;
			}
			// allow one mismatch
			if(genome.charAt(positionToCheck5 + i) != genome.charAt(positionToCheck3 + i)) {
				mismatch = true;
			}
		}
		ltr.setSequence5(genome.substring(ltr.getStart5(), (ltr.getStart5() + lengthLTR)));
		ltr.setSequence3(genome.substring(ltr.getStart3(), (ltr.getStart3() + lengthLTR)));
		ltr.setLength(lengthLTR);
	}

	/**
	 * TODO: would be interesting to explore how to simulate this mutation rate to align older LTRs (more degraded)
	 * lEx describes the length of a section between two LTRs that exactly matches. (the longest?)
	 * @param mutationRate this rate depends on the organism. 
	 * @return lEx
	 */
	private int lEx(double mutationRate) {
		return (int) Math.floor(LMIN / ((1 / mutationRate) * LMIN + 1));
	}
	
	private boolean leftMaximal(String s1, String s2) {
		return false;
	}
	
	private boolean rightMaximal(String s1, String s2) {
		return false;
	}
	
	private boolean maximalMatch(String s1, String s2) {
		return leftMaximal(s1, s2) && rightMaximal(s1, s2);
	}
	
	//

	
	/**
	 * This method checks the distance constraint between two genomic positions.
	 * @param i1 first genomic position
	 * @param i2 second genomic position
	 * @return 0 if the constraint is met, -1 if DMIN is met but not DMAX, 1 if DMAX is met but not DMIN
	 */
	private int distanceConstraint(int i1, int i2) {
		
		if(((i1 + DMIN) <= i2 ) && (i2 <= (i1 + DMAX))) {
			return 0;
		}
		
		if((i1 + DMIN) <= i2 ) {
			return -1;
		}
		if(i2 <= (i1 + DMAX)) {
			return 1;
		}
		
		return 42;
	}
	
	/**
	 * This step creates the FM-Index:
	 * @throws IOException 
	 */
	private FMIndexSingleSequence preprocessGenome() throws IOException {
		System.out.print("Using a lEx value of: " + LEX);
		
		FastaSequencesHandler openFile = new FastaSequencesHandler(); 
		QualifiedSequenceList genomeSequences = openFile.loadSequences(this.genomeFile);
		for (QualifiedSequence seq: genomeSequences) {
			genome += seq.getCharacters().toString();
		}
		FMIndexSingleSequence fmIndex = new FMIndexSingleSequence(genome);
		return fmIndex;
	}
	
	private void printProgress(int pos, int genomeLength) {
		double progress = (pos / genomeLength) * 100;
		if((progress % 10 == 0) && (progress != 0)) {
			System.out.println("Processed " + progress + " of the genome.");
		}
	}
	
	private void findCandidatePairs() {
		boolean nextPos = false;
		
		for(int i = 0; i < genome.length() - LEX; i++) {
			nextPos = false;
			
			String lookUpString = genome.substring(i, (i + LEX - 1));
			if(!DNASequence.isDNA(lookUpString)) {
				continue;
			}
			Object[] positions = fmIndex.search(lookUpString).toArray();
			if((positions.length < 2)) {
				// Must have at least two positions to build the pair.
				continue;
			}
			for(int j = 0; j < positions.length; j++) {
				int i1 = (int) positions[j];
				for(int k=j+1; k < positions.length - j; k++) {
					int i2 = (int) positions[k];
					int distance = distanceConstraint(i1, i2);
					if(distance == 0) {
						LTR candidate  = new LTR(i1, i2, lookUpString, lookUpString, idCount);
						candidatePairs.add(candidate);
						idCount++;
						nextPos = true;
					} else if(distance == 1) {
						/**
						 *  given that the positions are ordered, if the second position
						 *  is too far away, none of the others will be candidate pairs and thus
						 *  we can evaluate the next.
						 */
						break;
					}
				}
			}
			if(nextPos) {
				i = i + 20;
			}
		}
		
		return;
	}
	
	private void printCandidatePositions() {
		System.out.println(candidatePairs.size());
		for(LTR candidate: candidatePairs) {
			System.out.println(candidate.getStart5() + ", " + candidate.getStart3());
		}
	}
	
	
}
