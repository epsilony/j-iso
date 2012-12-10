/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package net.epsilony.levelset.ops;

import net.epsilony.math.CoordinatePartDiffFunction;
import net.epsilony.utils.geom.Coordinate;

/**
 *
 * @author epsilon
 */
public class LogisticNorm implements CoordinatePartDiffFunction {

    /**
     * $\frac{2s}{1+e^{-kt}}-s$
     *
     * @param ori
     * @param results
     * @param scale
     * @param k
     * @param partDiffOrder
     * @param dim
     */
    public static void valueOfLogistic(double[] ori, double[] results, double k, int partDiffOrder, int dim) {
        double t = ori[0];
        double e = Math.exp(-k * t);
        double b = 1 + e;
        double f = 2 / b - 1;
        results[0] = f;
        if (partDiffOrder > 0) {
            double f_t = e * k * 2 / (b * b);
            for (int i = 0; i < dim; i++) {
                results[i + 1] = ori[i + 1] * f_t;
            }
        }
    }
    CoordinatePartDiffFunction oriFun;
    double k;

    public LogisticNorm(CoordinatePartDiffFunction oriFun, double k) {
        this.oriFun = oriFun;
        this.k = k;
    }

    @Override
    public double[] values(Coordinate coord, double[] results) {
        results = oriFun.values(coord, results);
        valueOfLogistic(results, results, k, getDiffOrder(), getDim());
        return results;
    }

    @Override
    public void setDiffOrder(int order) {
        if (order > 1 || order < 0) {
            throw new UnsupportedOperationException("Only supports order 0 or 1.");
        }
        oriFun.setDiffOrder(order);
    }

    @Override
    public int getDiffOrder() {
        return oriFun.getDiffOrder();
    }

    @Override
    public int getDim() {
        return oriFun.getDim();
    }
    
}
