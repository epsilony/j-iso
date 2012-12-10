/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package net.epsilony.levelset.ops;

import net.epsilony.math.CoordinatePartDiffFunction;
import net.epsilony.spfun.CommonUtils;
import net.epsilony.utils.geom.Coordinate;

/**
 *
 * @author epsilon
 */
public class Scale implements CoordinatePartDiffFunction {

    /**
     *
     * @param ori
     * @param results
     * @param scale
     * @param k
     * @param partDiffOrder
     * @param dim
     */
    public static void valueOfScale(double[] ori, double[] results, double k, int partDiffOrder, int dim) {
        System.arraycopy(ori, 0, results, 0, CommonUtils.lenBase(dim, partDiffOrder));
    }
    CoordinatePartDiffFunction oriFun;
    double scale;

    public Scale(CoordinatePartDiffFunction oriFun, double scale) {
        this.oriFun = oriFun;
        this.scale = scale;
    }

    @Override
    public double[] values(Coordinate coord, double[] results) {
        results = oriFun.values(coord, results);
        valueOfScale(results, results, scale, getDiffOrder(), getDim());
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
