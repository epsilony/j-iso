/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package net.epsilony.levelset.functions;

import net.epsilony.math.CoordinatePartDiffFunction;
import net.epsilony.spfun.CommonUtils;
import net.epsilony.utils.geom.Coordinate;

/**
 *
 * @author epsilon
 */
public class EllipseLSF implements CoordinatePartDiffFunction{
    private int diffOrder;
    double xc,yc,a,b;
    private final double scale;


    public EllipseLSF(double xc, double yc, double a, double b,double scale) {
        this.xc = xc;
        this.yc = yc;
        this.a = a;
        this.b = b;
        this.scale = scale;

    }

    @Override
    public double[] values(Coordinate coord, double[] results) {
        if(null==results){
            results=new double[CommonUtils.len2DBase(diffOrder)];
        }
        double x=coord.x-xc;
        double y=coord.y-yc;
        double f=1-x*x/(a*a)-y*y/(b*b);
        results[0]=f*scale;
        if(diffOrder>0){
            results[1]=-2*x/(a*a)*scale;
            results[2]=-2*y/(b*b)*scale;
        }
        return results;
    }

    @Override
    public void setDiffOrder(int order) {
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
