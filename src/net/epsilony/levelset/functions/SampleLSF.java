/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package net.epsilony.levelset.functions;

import net.epsilony.levelset.AtomOperations;
import net.epsilony.math.CoordinatePartDiffFunction;
import net.epsilony.utils.geom.Coordinate;

/**
 *
 * @author epsilon
 */
public class SampleLSF {

    public CoordinatePartDiffFunction fun;
    public double x0, y0, w, h;
    public int b1;
    public int a1;
    public double ye1;
    public double xe1;
    public int b2;
    public int a2;
    public double ye2;
    public double xe2;
    public int b3;
    public int a3;
    public double ye3;
    public double xe3;
    public int r;
    public double yc;
    public double xc;
    public int ellScale;
    public EllipseLSF cirF;
    public EllipseLSF ellF3;
    public EllipseLSF ellF2;
    public EllipseLSF ellF1;
    public QuadrangleLSF quadF;

    private void init() {
        x0 = -1;
        y0 = 2.5;
        w = 100;
        h = 60;

        ellScale = 1;

        xe1 = 40;
        ye1 = 40;
        a1 = 20;
        b1 = 8;

        xe2 = 50;
        ye2 = 20;
        a2 = 30;
        b2 = 8;

        xe3 = 80;
        ye3 = 25;
        a3 = 15;
        b3 = 20;

        xc = 5;
        yc = 30;
        r = 10;
    }

    public final void generate() {
        Coordinate[] vertes = new Coordinate[4];
        vertes[0] = new Coordinate(x0, y0);
        vertes[1] = new Coordinate(x0 + w, y0);
        vertes[2] = new Coordinate(x0 + w, y0 + h);
        vertes[3] = new Coordinate(x0, y0 + h);


        quadF = new QuadrangleLSF(vertes, 5, 0.3);


        ellF1 = new EllipseLSF(xe1 + x0, ye1 + y0, a1, b1);


        ellF2 = new EllipseLSF(xe2 + x0, ye2 + y0, a2, b2);


        ellF3 = new EllipseLSF(xe3 + x0, ye3 + y0, a3, b3);


        cirF = new EllipseLSF(xc + x0, yc + y0, r, r);
        CoordinatePartDiffFunction t1 = AtomOperations.union(0,ellF1, ellF2, 2);
        t1 = AtomOperations.union(0,t1, ellF3, 2);
        t1 = AtomOperations.union(0,t1, cirF, 2);

        this.fun = AtomOperations.intersection(0,quadF, 1, t1, -1, 2);

    }

    public double[] values(Coordinate coord, double[] results) {
        return fun.values(coord, results);
    }

    public void setDiffOrder(int order) {
        fun.setDiffOrder(order);
    }

    public int getDiffOrder() {
        return fun.getDiffOrder();
    }

    public int getDim() {
        return fun.getDim();
    }

    public SampleLSF() {
        init();
        generate();
    }

    protected SampleLSF(boolean generate) {
        init();
        if (generate) {
            generate();
        }
    }
}
