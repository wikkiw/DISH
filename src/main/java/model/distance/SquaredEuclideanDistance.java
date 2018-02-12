package model.distance;

/**
 *
 * This class serves for calculating the euclidean distance
 * 
 * @author wiki on 30/03/2017
 */
public class SquaredEuclideanDistance implements Distance {
    
    @Override
    public double getDistance(double[] point_a, double[] point_b) {
        
        if(point_a == null || point_b == null || point_a.length != point_b.length) {
            return Double.MAX_VALUE;
        }
        
        double distance = 0;
        
        for(int i = 0; i < point_a.length; i++) {
            distance += Math.pow(point_a[i]-point_b[i], 2);
        }
        
        return distance;
        
    }
    
}
