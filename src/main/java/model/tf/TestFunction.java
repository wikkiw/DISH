package model.tf;

import model.Individual;

/**
 * Basic interface for test functions.
 * Created by wiki on 12/02/2018
 */
public interface TestFunction {

    double fitness(Individual individual);

    /**
     * Convenience override for fitness(Individiual individual)
     *
     * @param vector
     * @return fitness for given vector
     */
    double fitness(double[] vector);

    void constrain(Individual individual);

    double[] generateTrial(int dim);

    double fixedAccLevel();

    double optimum();

    double max(int dim);

    double min(int dim);
    
    String name();

}
