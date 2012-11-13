/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package net.epsilony.levelset.functions;

import net.epsilony.levelset.AtomOperations;
import net.epsilony.math.CoordinatePartDiffFunction;
import net.epsilony.utils.geom.Coordinate;
import net.epsilony.utils.geom.Quadrangle;

/**
 *
 * @author epsilon
 */
public class QuadrangleLSF extends Quadrangle implements CoordinatePartDiffFunction {

    private CoordinatePartDiffFunction fun;

    public QuadrangleLSF(int m,Coordinate[] vertes,double scale, double k) {
        super(vertes, true);
        init(m,scale,k);
    }
    
    public QuadrangleLSF(Coordinate[] vertes,double scale, double k) {
        super(vertes, true);
        init(0,scale,k);
    }
    
    private void init(int m,double scale,double k){
        CoordinatePartDiffFunction[] funs = new CoordinatePartDiffFunction[4];
        for (int i = 0; i < vertes.length; i++) {
            Coordinate p1 = vertes[i % 4];
            Coordinate p2 = vertes[(i + 1) % 4];
            funs[i]=new HalfPlaneLogisticLST(p1, p2, scale, k);
        }
        CoordinatePartDiffFunction itsec1 = AtomOperations.intersection(m,funs[0],1 ,funs[2],1, 2);
        CoordinatePartDiffFunction itsec2 = AtomOperations.intersection(m,funs[1],1 ,funs[3],1, 2);
        fun = AtomOperations.intersection(m,itsec1, itsec2, 2);
    }

    @Override
    public double[] values(Coordinate coord, double[] results) {
        return fun.values(coord, results);
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
