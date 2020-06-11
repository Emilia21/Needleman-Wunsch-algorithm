/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package needleman.wunsch;

import java.util.Scanner;

/**
 *
 * @author epetkova2
 */
public class NeedlemanWunsch {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        Scanner scanner = new Scanner(System.in);
        String firstSequence = scanner.nextLine();
        String secondSequence= scanner.nextLine();
        
        Needleman_Wunsch nw = new Needleman_Wunsch(firstSequence, secondSequence);
        nw.printBacktraceArray();
        nw.printScoreArray();
        nw.printAlignedSequences();
    }
    
}
