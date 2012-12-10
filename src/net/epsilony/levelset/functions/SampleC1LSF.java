/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package net.epsilony.levelset.functions;

import net.epsilony.levelset.ops.AtomOperations;
import net.epsilony.math.CoordinatePartDiffFunction;
import net.epsilony.utils.geom.Coordinate;

/**
 *
 * @author epsilon
 */
public class SampleC1LSF extends SampleLSF{
    public SampleC1LSF() {
        super(false);
        generateC1();
    }
    
    public final void generateC1(){
                Coordinate[] vertes = new Coordinate[4];
        vertes[0] = new Coordinate(x0, y0);
        vertes[1] = new Coordinate(x0 + w, y0);
        vertes[2] = new Coordinate(x0 + w, y0 + h);
        vertes[3] = new Coordinate(x0, y0 + h);


        quadF = new QuadrangleLSF(1,vertes, 0.3);


        ellF1 = new EllipseLSF(xe1 + x0, ye1 + y0, a1, b1);


        ellF2 = new EllipseLSF(xe2 + x0, ye2 + y0, a2, b2);


        ellF3 = new EllipseLSF(xe3 + x0, ye3 + y0, a3, b3);


        cirF = new EllipseLSF(xc + x0, yc + y0, r, r);
//        CoordinatePartDiffFunction t1 = AtomOperations.union(1,ellF1,ellF2,1, 2);
//        t1 = AtomOperations.union(1,t1, 0.5,ellF3,1, 2);
//        t1 = AtomOperations.union(1,t1, 0.5/1250,cirF,1, 2);
//
//        this.fun = AtomOperations.intersection(1,quadF, 1/50.0, t1, -1/1250.0, 2);
    }
}
