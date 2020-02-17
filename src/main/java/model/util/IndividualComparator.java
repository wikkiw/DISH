package model.util;

import java.util.Comparator;
import model.Individual;

/**
 *
 * @author wiki on 17/06/2016
 */
public class IndividualComparator implements Comparator<Individual> {
    
    @Override
        public int compare(Individual t, Individual t1) {
            
            if(t.fitness < t1.fitness) {
                return -1;
            }
            else if(t.fitness == t1.fitness){
                return 0;
            }
            else {
                return 1;
            }
            
        }
    
}
