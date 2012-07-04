/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package net.epsilony.isomesh;

/**
 *
 * @author epsilon
 */
public interface IsoFunction {
    /**
     * determine the [f(x) grad(f(x)) ...] the content of return is set by {@code partDiffOrder}
     * @param 
     * @param partDiffOrder the highest partial differential order of results 
     * @param result for saving result, can be null
     * @return {@code result} or new array contains results if {@code result} is {@code null}
     */
    double [] value(double[] vec, int partDiffOrder, double[] result);
    
    int getDimension();
}
