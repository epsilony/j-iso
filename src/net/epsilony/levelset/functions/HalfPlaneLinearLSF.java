/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package net.epsilony.levelset.functions;

import net.epsilony.math.CoordinatePartDiffFunction;
import net.epsilony.spfun.CommonUtils;
import net.epsilony.utils.geom.Coordinate;
import net.epsilony.utils.geom.GeometryMath;

/**
 *
 * @author epsilon
 */
public class HalfPlaneLinearLSF implements CoordinatePartDiffFunction{
    int diffOrder;
    Coordinate p1,p2;
    private final Coordinate normalize;
    private final double slope;

    public HalfPlaneLinearLSF(Coordinate p1, Coordinate p2,double slope) {
        this.p1 = new Coordinate(p1);
        this.p2 = new Coordinate(p2);
        this.normalize = GeometryMath.normalize(new Coordinate(p2.x-p1.x,p2.y-p1.y));
        this.slope = slope;
    }
    
    
    @Override
    public double[] values(Coordinate coord, double[] results) {
        if(null==results){
            results=new double[CommonUtils.len2DBase(diffOrder)];
        }
        double v=slope*GeometryMath.cross2D(normalize.x, normalize.y, coord.x-p1.x, coord.y-p1.y);
        results[0]=v;
        if(diffOrder>0){
            results[1]=-normalize.y*slope;
            results[2]=normalize.x*slope;
        }
        return results;
    }

    @Override
    public void setDiffOrder(int order) {
        if(diffOrder>1){
            throw new UnsupportedOperationException("Only support order 0 or 1");
        }
        this.diffOrder=order;
    }

    @Override
    public int getDiffOrder() {
        return diffOrder;
    }

    @Override
    public int getDim() {
        return 2;
    }
    
}
