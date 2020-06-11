/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package needleman.wunsch;

/**
 * Needleman-Wunsch alignment of two sequences entered on command line The major
 * goal of computational sequence analysis is to predict the function and
 * structure of genes and proteins from their sequence. The algorithm consists
 * of two major steps: 1) generate of the alignment F-matrix 2) The trace-back
 * calls generating final alignment sequence.
 *
 * @author epetkova2
 */
public class Needleman_Wunsch {

    /**
     * @param GAP This is the gap penalty
     * @param MATCH This is the match award
     * @param MISMATCH This is the mismatch penalty
     * @param firstSequence This is the first sequence
     * @param secondSequence This is the second sequence
     * @param scoreArray This is the the score matrix
     * @param backtraceArray This is the the back trace matrix
     * @param firstSequenceAligned This is the first sequence after the
     * alignment
     * @param secondSequenceAligned This is the second sequence after the
     * alignment
     *
     */
    final static int GAP = -5;
    final static int MATCH = 10;
    final static int MISMATCH = -2;

    private String firstSequence;
    private String secondSequence;
    private int firstSequenceLength, secondSequenceLength; // their lengths
    private int[][] scoreArray;
    private String[][] backtraceArray;
    private String firstSequenceAligned = "";
    private String secondSequenceAligned = "";

    /**
     * This is the constructor which needs two parameters: the first sequence
     * the second sequence and run the method for filling the matrix
     */
    public Needleman_Wunsch(String firstSequence, String secondSequence) {
        this.firstSequence = firstSequence;
        this.secondSequence = secondSequence;
        this.firstSequenceLength = this.firstSequence.length();
        this.secondSequenceLength = this.secondSequence.length();
        this.scoreArray = new int[this.secondSequenceLength + 1][this.firstSequenceLength + 1];
        this.backtraceArray = new String[this.secondSequenceLength][this.firstSequenceLength];
        this.fillArrays();
        this.alignment();
    }

    /**
     * This method fills the score array as it first fills the first row second
     * fills the first column and after that fills the rest of the array
     *
     * @param diagScore this is the element on (row -1, col -1) possition from
     * the current possition
     * @param upScore this is the element upper from the current one
     * @param leftScore this is the element left from the current one
     * @param maxScore this is the maximum from diagScore, upScore, leftScore
     *
     */
    public void fillArrays() {
        int diagScore;
        int topScore;
        int leftScore;
        int maxScore;
        String step;

        // Fill the first row
        for (int col = 0; col <= this.firstSequenceLength; col++) {
            scoreArray[0][col] = GAP * col;
        }

        // Fill the first column  
        for (int row = 0; row <= this.secondSequenceLength; row++) {
            scoreArray[row][0] = GAP * row;
        }

        // Fill the rest of the array:
        for (int row = 1; row <= this.secondSequenceLength; row++) {
            for (int col = 1; col <= this.firstSequenceLength; col++) {
                if (this.firstSequence.charAt(col - 1) == this.secondSequence.charAt(row - 1)) {
                    diagScore = scoreArray[row - 1][col - 1] + MATCH;
                } else {
                    diagScore = scoreArray[row - 1][col - 1] + MISMATCH;
                }
                leftScore = scoreArray[row][col - 1] + GAP;
                topScore = scoreArray[row - 1][col] + GAP;
                if (diagScore >= leftScore && diagScore >= topScore) {
                    maxScore = diagScore;
                    step = "D";
                } else if (leftScore >= diagScore && leftScore >= topScore) {
                    maxScore = leftScore;
                    step = "L";
                } else {
                    maxScore = topScore;
                    step = "T";
                }
                backtraceArray[row - 1][col - 1] = step;
                scoreArray[row][col] = maxScore;
            }
        }
    }

    public void alignment() {
        int currentRow = this.secondSequenceLength - 1;
        int currentCol = this.firstSequenceLength - 1;

        while (currentRow >= 0 && currentCol >= 0) {
            String step = backtraceArray[currentRow][currentCol];
            switch (step) {
                case "D":
                    currentRow--;
                    currentCol--;
                    this.firstSequenceAligned += this.firstSequence.charAt(this.firstSequenceLength - 1);
                    this.firstSequence = this.firstSequence.substring(0, this.firstSequenceLength - 1);
                    this.firstSequenceLength = this.firstSequence.length();
                    this.secondSequenceAligned += this.secondSequence.charAt(this.secondSequenceLength - 1);
                    this.secondSequence = this.secondSequence.substring(0, this.secondSequenceLength - 1);
                    this.secondSequenceLength = this.secondSequence.length();
                    break;
                case "T":
                    currentRow--;
                    this.firstSequenceAligned += "_";
                    this.secondSequenceAligned += this.secondSequence.charAt(this.secondSequenceLength - 1);
                    this.secondSequence = this.secondSequence.substring(0, this.secondSequenceLength - 1);
                    this.secondSequenceLength = this.secondSequence.length();
                    break;
                case "L":
                    currentCol--;
                    this.firstSequenceAligned += this.firstSequence.charAt(this.firstSequenceLength - 1);
                    this.firstSequence = this.firstSequence.substring(0, this.firstSequenceLength - 1);
                    this.secondSequenceAligned += "_";
                    this.firstSequenceLength = this.firstSequence.length();
                    break;

            }
        }
        this.firstSequenceAligned = new StringBuilder(firstSequenceAligned).reverse().toString();
        this.secondSequenceAligned = new StringBuilder(secondSequenceAligned).reverse().toString();
    }

    /**
     * This method prints the score matrix
     */
    public void printScoreArray() {
        for (int row = 0; row < scoreArray.length; row++) {
            for (int col = 0; col < scoreArray[row].length; col++) {
                System.out.print(scoreArray[row][col] + " ");
            }
            System.out.println();
        }
    }

    /**
     * This method prints the back trace matrix
     */
    public void printBacktraceArray() {
        for (int row = 0; row < backtraceArray.length; row++) {
            for (int col = 0; col < backtraceArray[row].length; col++) {
                System.out.print(backtraceArray[row][col] + " ");
            }
            System.out.println();
        }
    }

    /**
     * This method prints the back trace matrix
     */
    public void printAlignedSequences() {
        System.out.println(this.firstSequenceAligned);
        System.out.println(this.secondSequenceAligned);
    }

}
