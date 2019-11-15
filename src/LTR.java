
public class LTR {

	private String sequence5;
	private String sequence3;
	private int start5;		// global position to genome
	private int start3;	
	private int end5;
	private int end3;
	private int lengthLTR;
	private int id;
	private double homologyFraction;
	private boolean pol = false;
	private boolean gag = false;
	private boolean TG = false;
	private boolean CA = false;
	
	public LTR(int _start5, int _start3, String _sequence5, String _sequence3, int id) {
		this.start5 = _start5;
		this.start3 = _start3;
		this.sequence5 = _sequence5;
		this.sequence3 = _sequence3;
	}
	
	public void setPol(boolean _pol) {
		this.pol = _pol;
	}
	
	public void setGag(boolean _gag) {
		this.gag = _gag;
	}
	
	public int getStart5() {
		return start5;
	}
	
	public int getStart3() {
		return start3;
	}
	public int getEnd5() {
		return end5;
	}
	
	public int getEnd3() {
		return end3;
	}
	
	public int getLTRLength() {
		return lengthLTR;
	}
	
	public int getTotalLength() {
		return end3 - start5;
	}
	
	public boolean getPol() {
		return pol;
	}
	
	public boolean getGag() {
		return gag;
	}
	
	public boolean getTG() {
		return TG;
	}
	
	public boolean getCA() {
		return CA;
	}
	
	public double[] getCharVector() {
		double[] charVector = new double[5];
		charVector[0] = homologyFraction;
		if(pol) {
			charVector[1] = 1;
		} else {
			charVector[1] = 0;
		}
		if(gag) {
			charVector[2] = 1;
		} else {
			charVector[2] = 0;
		}
		if(TG) {
			charVector[3] = 1;
		} else {
			charVector[3] = 0;
		}
		if(CA) {
			charVector[4] = 1;
		} else {
			charVector[4] = 0;
		}
		
		return charVector;
	}
	
	public void setLength(int _length) {
		this.lengthLTR = _length;
		this.end3 = start3 + lengthLTR;
		this.end5 = start5 + lengthLTR;
		this.homologyFraction = 1 / lengthLTR;
		if((sequence5.substring(0,1) == "TG") || (sequence3.substring(0,1) == "TG")){
			this.TG = true;
		}
		if((sequence5.substring(lengthLTR-1, lengthLTR) == "CA") || (sequence3.substring(lengthLTR-1, lengthLTR) == "CA")) {
			this.CA = true;
		}
	}
	
	public void setSequence5(String seq) {
		this.sequence5 = seq;
	}
	public void setSequence3(String seq) {
		this.sequence3 = seq;
	}

	public int getID() {
		return id;
	}
	
	public String getSeq5() {
		return sequence5;
	}
	
}
