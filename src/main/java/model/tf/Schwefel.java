package model.tf;

import model.Individual;
import model.util.IndividualUtil;
import model.util.Random;
import model.util.UniformRandom;

/**
 * Created by jakub on 27/10/15.
 */
public class Schwefel implements TestFunction {

    @Override
    public double fitness(Individual individual) {
        return fitness(individual.vector);
    }

    @Override
    public double fitness(double[] vector) {
        double D = vector.length;
        double s1 = 0;

        for (int i = 0; i < vector.length; i++) {
            s1 += vector[i] * Math.sin(Math.sqrt(Math.abs(vector[i])));
        }
        return 418.9829*D - s1;
    }

    @Override
    public void constrain(Individual individual) {
        IndividualUtil.clipInBounds(individual, -500, 500);
    }

    @Override
    public double[] generateTrial(int dim) {
        double[] vector = new double[dim];
        Random rnd = new UniformRandom();
        for (int i = 0; i < dim; i++) vector[i] = rnd.nextDouble(-500, 500);
        return vector;
    }

    @Override
    public double fixedAccLevel() {
        return 10E-7;
    }

    @Override
    public double optimum() {
        return 0.0;
    }

    @Override
    public double max(int dim) {
        return 500;
    }

    @Override
    public double min(int dim) {
        return -500;
    }

    @Override
    public String name() {
        return "Schwefel";
    }
}
