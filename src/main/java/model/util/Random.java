package model.util;

/**
 *
 * @author wiki on 24/11/2015
 */
public interface Random {

    double nextDouble();

    int nextInt(int bound);

    default double nextDouble(double min, double max) {
        return (nextDouble() * (max - min) + min);
    }

    @Override
    public String toString();

}
