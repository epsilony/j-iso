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
public class Complement implements CoordinatePartDiffFunction {

    public static void valueOfComplement(double[] ori, double[] results) {
        for (int i = 0; i < ori.length; i++) {
            results[i] = -ori[i];
        }
    }
    private final CoordinatePartDiffFunction fun;

    public Complement(CoordinatePartDiffFunction fun) {
        this.fun = fun;
    }

    @Override
    public double[] values(Coordinate coord, double[] results) {
        fun.values(coord, results);
        for (int i = 0; i < results.length; i++) {
            results[i] *= -1;
        }
        return results;
    }

    @Override
    public void setDiffOrder(int order) {
        fun.setDiffOrder(order);
    }

    @Override
    public int getDiffOrder() {
        return fun.getDiffOrder();
    }

    @Override
    public int getDim() {
        return fun.getDim();
    }
    
}
