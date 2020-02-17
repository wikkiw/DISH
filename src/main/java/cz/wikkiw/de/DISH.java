package cz.wikkiw.de;

import java.util.ArrayList;
import java.util.Date;
import java.util.List;
import java.util.Locale;
import java.util.stream.DoubleStream;
import model.Individual;
import model.tf.TestFunction;
import org.apache.commons.math3.stat.descriptive.moment.Mean;
import org.apache.commons.math3.stat.descriptive.moment.StandardDeviation;
import org.apache.commons.math3.stat.descriptive.rank.Median;
import model.util.IndividualComparator;
import model.util.OtherDistributionsUtil;
import model.distance.EuclideanDistance;
import model.tf.Cec2015;
import model.util.UniformRandom;

/**
 * 
 * DISH algorithm - multithread ready
 * 
 * @author adam on 17/02/2020
 */
public class DISH implements Runnable {

    int D;
    int G;
    int NP;
    List<Individual> Aext;
    List<Individual> P;
    int FES;
    int MAXFES;
    TestFunction f;
    Individual best;
    List<Individual> bestHistory;
    double[] M_F;
    double[] M_CR;
    List<Double> S_F;
    List<Double> S_CR;
    int H;
    model.util.Random rndGenerator = new UniformRandom();
    int id;
    int Asize;
    List<double[]> M_Fhistory;
    List<double[]> M_CRhistory;
    
    protected final int minPopSize;
    protected final int maxPopSize;

    protected List<double[]> imp_hist;
    
    /**
     * 
     * @param D
     * @param MAXFES
     * @param f
     * @param H
     * @param NP
     * @param minPopSize 
     */
    public DISH(int D, int MAXFES, TestFunction f, int H, int NP, int minPopSize) {
        
        this.D = D;
        this.MAXFES = MAXFES;
        this.f = f;
        this.H = H;
        this.NP = NP;

        this.minPopSize = minPopSize;
        this.maxPopSize = NP;
        this.imp_hist = new ArrayList<>();
    }

    public String getName() {
        return "DISH";
    }
    
    /**
     *
     * @param ind
     * @return
     */
    protected boolean isBest(Individual ind) {

        if (this.best == null || ind.fitness < this.best.fitness) {
            this.best = ind;
            return true;
        }

        return false;

    }
    
    /**
     * Creation of initial population.
     */
    protected void initializePopulation(){
        
        /**
         * Initial population
         */
        id = 0;
        double[] features;
        this.P = new ArrayList<>();
        Individual ind;

        for (int i = 0; i < this.NP; i++) {
            id = i;
            features = this.f.generateTrial(this.D).clone();
            ind = new Individual(String.valueOf(id), features, this.f.fitness(features));
            this.P.add(ind);
            this.FES++;
            if(this.isBest(ind)) {
                this.writeHistory();
            }
            
        }
        
    }
    
    /**
     *
     */
    protected void writeHistory() {
        
        this.imp_hist.add(new double[]{this.FES, this.best.fitness});

    }
    
    /**
     * For multithread purposes
     */
    @Override
    public void run() {
        this.runAlgorithm();
    }
    
    /**
     * 
     * Main algorithm function
     * 
     * @return 
     */
    public Individual runAlgorithm() {

        /**
         * Initialization
         */
        this.G = 0;
        this.Aext = new ArrayList<>();
        this.best = null;
        this.bestHistory = new ArrayList<>();
        

        /**
         * Initial population
         */
        initializePopulation();

        this.M_F = new double[this.H];
        this.M_CR = new double[this.H];

        for (int h = 0; h < this.H; h++) {
            this.M_F[h] = 0.5;
            this.M_CR[h] = 0.8;
            
            if(h == this.H-1) {
                this.M_F[h] = 0.9;
                this.M_CR[h] = 0.9;
            }
        }

        int r, Psize, pbestIndex, memoryIndex, k = 0;
        double Fg, CRg, Fw, gg, pmin = 0.125, pmax = 0.25, p, wSsum, meanS_F1, meanS_F2, meanS_CR1, meanS_CR2;
        Individual trial, x;
        double[] v, pbest, pr1, pr2, u, wsList = new double[NP], SFlist = new double[NP], SCRlist = new double[NP];
        int[] rIndexes;
        double[][] parents;
        List<Individual> newPop, pBestArray;

        EuclideanDistance euclid = new EuclideanDistance();
        
        /**
         * Generation iteration;
         */
        while (true) {

            this.G++;
            newPop = new ArrayList<>();
                
            memoryIndex = 0;

            p = ((pmax - pmin)/(double) this.MAXFES) * this.FES + pmin;

            for (int i = 0; i < this.NP; i++) {

                x = this.P.get(i);
                r = rndGenerator.nextInt(this.H);
                Fg = OtherDistributionsUtil.cauchy(this.M_F[r], 0.1);
                while (Fg <= 0) {
                    Fg = OtherDistributionsUtil.cauchy(this.M_F[r], 0.1);
                }
                if(Fg > 1) {
                    Fg = 1;
                }


                CRg = OtherDistributionsUtil.normal(this.M_CR[r], 0.1);
                if(CRg > 1) {
                    CRg = 1;
                }
                else if(CRg < 0) {
                    CRg = 0;
                }

                /**
                 * new in jSO
                 */
                gg = (double) this.FES / (double) this.MAXFES;
                if(gg < 0.25) {

                    if(CRg < 0.7) {
                        CRg = 0.7;
                    }

                }
                else if(gg < 0.5) {

                    if(CRg < 0.6) {
                        CRg = 0.6;
                    }

                }

                if(gg < 0.6 && Fg > 0.7) {
                    Fg = 0.7;
                }

                Psize = (int) (rndGenerator.nextDouble(p, 0.2) * this.NP);
                if(Psize < 2) {
                    Psize = 2;
                }

                pBestArray = this.resize(this.P, Psize);

                /**
                 * Parent selection
                 */
                pbestIndex = this.getIndexOfRandBestFromList(pBestArray, x.id);
                pbest = this.P.get(pbestIndex).vector.clone();
                rIndexes = this.genRandIndexes(i, this.NP, this.NP + this.Aext.size(), pbestIndex);
                pr1 = this.P.get(rIndexes[0]).vector.clone();
                if(rIndexes[1] > this.NP - 1) {
                    pr2 = this.Aext.get(rIndexes[1] - this.NP).vector.clone();
                }else {
                    pr2 = this.P.get(rIndexes[1]).vector.clone();
                }
                parents = new double[4][D];
                parents[0] = x.vector;
                parents[1] = pbest;
                parents[2] = pr1;
                parents[3] = pr2;

                if(gg < 0.2) {
                    Fw = 0.7 * Fg;
                }
                else if(gg < 0.4) {
                    Fw = 0.8 * Fg;
                }
                else {
                    Fw = 1.2 * Fg;
                }

                /**
                 * Mutation
                 */               
                v = mutationDISH(parents, Fg, Fw);

                /**
                 * Crossover
                 */
                u = crossover(CRg, v, x.vector);

                /**
                 * Constrain check
                 */
                u = constrainCheck(u, x.vector);

                /**
                 * Trial ready
                 */
                id++;
                trial = new Individual(String.valueOf(id), u, f.fitness(u));

                /**
                 * Trial is better
                 */
                if(trial.fitness <= x.fitness) {
                    newPop.add(trial);
                    this.Aext.add(x);

                    SFlist[memoryIndex] = Fg;
                    SCRlist[memoryIndex] = CRg;
                    /**
                     * inverse distance
                     */
                    wsList[memoryIndex] = euclid.getDistance(x.vector, trial.vector);

                    memoryIndex++;

                }else {
                    newPop.add(x);
                }

                this.FES++;
                if(this.isBest(trial)) {
                    this.writeHistory();
                }
                if (this.FES >= this.MAXFES) {
                    break;
                }

                this.Aext = this.resizeAext(this.Aext, this.NP);

            }

            if(this.FES >= this.MAXFES) {
                break;
            }

            /**
             * Memories update - new
             */
            if(memoryIndex > 0) {
                wSsum = 0;
                for(int i = 0; i < memoryIndex; i++) {
                    wSsum += wsList[i];
                }

                meanS_F1 = 0;
                meanS_F2 = 0;
                meanS_CR1 = 0;
                meanS_CR2 = 0;

                for (int s = 0; s < memoryIndex; s++) {
                    meanS_F1 += (wsList[s] / wSsum) * SFlist[s] * SFlist[s];
                    meanS_F2 += (wsList[s] / wSsum) * SFlist[s];
                    meanS_CR1 += (wsList[s] / wSsum) * SCRlist[s] * SCRlist[s];
                    meanS_CR2 += (wsList[s] / wSsum) * SCRlist[s];
                }

                if(meanS_F2 != 0) {
                    this.M_F[k] = ((meanS_F1 / meanS_F2) + this.M_F[k])/2;
                }
                if(meanS_CR2 != 0) {
                    this.M_CR[k] = ((meanS_CR1 / meanS_CR2) + this.M_CR[k])/2;
                }

                k++;
                if (k >= (this.H - 1)) {
                    k = 0;
                }
            }
            
            /**
             * Resize of population and archives
             */
            this.P = new ArrayList<>();
            this.P.addAll(newPop);
            NP = (int) Math.round(this.maxPopSize - ((double) this.FES/(double) this.MAXFES)*(this.maxPopSize - this.minPopSize));
            P = this.resizePop(P, NP);

        }

        return this.best;

    }
    
    /**
     *
     * @param list
     * @param size
     * @return
     */
    protected List<Individual> resizeAext(List<Individual> list, int size) {

        if(size <= 0) {
            return new ArrayList<>();
        }
        
        if(size >= list.size()){
            return list;
        }

        List<Individual> toRet = new ArrayList<>();
        toRet.addAll(list);
        int index;

        while (toRet.size() > size) {

            index = rndGenerator.nextInt(toRet.size());
            toRet.remove(index);

        }

        return toRet;

    }
    
    /**
     * 
     * @param u
     * @param x
     * @return 
     */
    protected double[] constrainCheck(double[] u, double[] x){
        /**
         * Constrain check
         */
        for (int d = 0; d < this.D; d++) {
            if (u[d] < this.f.min(this.D)) {
                u[d] = (this.f.min(this.D) + x[d]) / 2.0;
            } else if (u[d] > this.f.max(this.D)) {
                u[d] = (this.f.max(this.D) + x[d]) / 2.0;
            }
        }
        
        return u;
    }
    
    /**
     * 
     * @param CR
     * @param v
     * @param x
     * @return 
     */
    protected double[] crossover(double CR, double[] v, double[] x){
        
        /**
         * Crossover
         */
        double[] u = new double[this.D];
        int jrand = rndGenerator.nextInt(this.D);

        for (int j = 0; j < this.D; j++) {
            if (getRandomCR() <= CR || j == jrand) {
                u[j] = v[j];
            } else {
                u[j] = x[j];
            }
        }
        
        return u;
        
    }
    
    /**
     * 
     * @return 
     */
    protected double getRandomCR(){
        return rndGenerator.nextDouble();
    }
    
    /**
     *
     * @param index
     * @param max1
     * @param max2
     * @param pbest
     * @return
     */
    protected int[] genRandIndexes(int index, int max1, int max2, int pbest) {

        int a, b;

        a = rndGenerator.nextInt(max1);
        
        while(a == pbest || a == index){
            a = rndGenerator.nextInt(max1);
        }
        
        b = rndGenerator.nextInt(max2);

        while (b == a || b == index || b == pbest) {
            b = rndGenerator.nextInt(max2);
        }

        return new int[]{a, b};
    }
    
    /**
     *
     * @param list
     * @param size
     * @return
     */
    protected List<Individual> resize(List<Individual> list, int size) {

        List<Individual> toRet = new ArrayList<>();
        toRet.addAll(list);
        
        toRet.sort(new IndividualComparator());
        toRet = toRet.subList(0, size);

        return toRet;

    }
    
    /**
     *
     * @param list
     * @param id
     * @return
     */
    protected int getIndexOfRandBestFromList(List<Individual> list, String id) {
        
        int index = rndGenerator.nextInt(list.size());
        
        while(list.get(index).id.equals(id)) {
            index = rndGenerator.nextInt(list.size());
        }

        return index;

    }
    
    /**
     *
     * @param parentArray
     * @param F
     * @return
     */
    protected double[] mutationRand1(double[][] parentArray, double F) {

        double[] u = new double[D];
        double[] a = parentArray[1];
        double[] b = parentArray[2];
        double[] c = parentArray[3];

        for (int i = 0; i < D; i++) {

            u[i] = a[i] + F * (b[i] - c[i]);

        }

        return u;

    }
    
    /**
     *
     * List of parents for mutation x, a, b, c
     *
     * @param xIndex
     * @return
     */
    protected double[][] getParents(int xIndex) {

        int r1, r2, r3;

        r1 = rndGenerator.nextInt(NP);
        
        while(r1 == xIndex){
            r1 = rndGenerator.nextInt(NP);
        }
        
        r2 = rndGenerator.nextInt(NP);

        while (r2 == r1 || r2 == xIndex) {
            r2 = rndGenerator.nextInt(NP);
        }
        
        r3 = rndGenerator.nextInt(NP);

        while (r3 == r2 || r3 == r1 || r3 == xIndex) {
            r3 = rndGenerator.nextInt(NP);
        }
        
        double[][] parrentArray = new double[4][D];

        parrentArray[0] = P.get(xIndex).vector;
        parrentArray[1] = P.get(r1).vector;
        parrentArray[2] = P.get(r2).vector;
        parrentArray[3] = P.get(r3).vector;

        return parrentArray;

    }
    
    /**
     * 
     * @param parents
     * @param F
     * @param Fw
     * @return 
     */
    protected double[] mutationDISH(double[][] parents, double F, double Fw){
        
        /**
         * Parents:
         * x
         * pbest
         * pr1
         * pr2
         */
        
        double[] v = new double[this.D];
        for (int j = 0; j < this.D; j++) {

            v[j] = parents[0][j] + Fw * (parents[1][j] - parents[0][j]) + F * (parents[2][j] - parents[3][j]);

        }
        
        return v;
        
    }
    
    /**
     *
     * @param list
     * @param size
     * @return
     */
    protected List<Individual> resizePop(List<Individual> list, int size) {

        if(size == list.size()){
            return list;
        }
        
        List<Individual> toRet = new ArrayList<>();
        toRet.addAll(list);
        toRet.sort(new IndividualComparator());
        toRet = toRet.subList(0, size);

        return toRet;

    }

    /**
     * 
     * @return 
     */
    public List<double[]> getImp_hist() {
        return imp_hist;
    }

    /**
     * 
     * @param vector
     * @return 
     */
    protected String getNiceVector(double[] vector) {
        
        String out = "";
        out += "{";
        
        for(int i = 0; i < vector.length; i++) {
            out += String.format(Locale.US, "%.4f", vector[i]);
            if(i != vector.length-1)
                out +=", ";
        }
        
        out += "}";
        
        return out;
    }
    
    /**
     * 
     * @return 
     */
    public Individual getBest() {
        return best;
    }
    
    
    
    /**
     * @param args the command line arguments
     * @throws java.lang.Exception
     */
    public static void main(String[] args) throws Exception {
        
        int dimension = 10;
        int NP = (int) (25*Math.log(dimension)*Math.sqrt(dimension));
        int minNP = (int) (25*Math.log(dimension)*Math.sqrt(dimension));
        int MAXFES = 10000 * dimension;
        int H = 5;
        
        int funcNumber = 6;
        TestFunction tf = new Cec2015(dimension, funcNumber);

        DISH shade;

        int runs = 10;
        double[] bestArray = new double[runs];

        System.out.println("START: " + new Date());
        
        for (int k = 0; k < runs; k++) {

            shade = new DISH(dimension, MAXFES, tf, H, NP, minNP);

            shade.runAlgorithm();

            bestArray[k] = shade.getBest().fitness - tf.optimum();
            System.out.println(shade.getBest().fitness - tf.optimum());
            System.out.println(shade.getNiceVector(shade.getBest().vector));
            System.out.println("Best solution found in: " + shade.getImp_hist().get(shade.getImp_hist().size()-1)[0]);

        }

        System.out.println("END: " + new Date());
        
        System.out.println("=================================");
        System.out.println("Min: " + DoubleStream.of(bestArray).min().getAsDouble());
        System.out.println("Max: " + DoubleStream.of(bestArray).max().getAsDouble());
        System.out.println("Mean: " + new Mean().evaluate(bestArray));
        System.out.println("Median: " + new Median().evaluate(bestArray));
        System.out.println("Std. Dev.: " + new StandardDeviation().evaluate(bestArray));
        System.out.println("=================================");

    }
    
}
