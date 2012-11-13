/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package net.epsilony.levelset.functions;

import net.epsilony.utils.geom.Coordinate;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 *
 * @author epsilon
 */
public class EllipseLSFTest {

    public EllipseLSFTest() {
    }

    @Test
    public void testSomeMethod() {
        double x0 = 11.1;
        double y0 = -2.3;
        double a = 44;
        double b = 12;
        double scale = 3.3;

        EllipseLSF ell = new EllipseLSF(x0, y0, a, b, scale);
        ell.setDiffOrder(1);
        double[] samples = new double[]{x0, y0, x0 + a, y0, x0, y0 - b, 10, 4};
        double[][] exps = new double[][]{
            {scale, 0, 0},
            {0, -2 / a * scale, 0},
            {0, 0, 2 / b * scale},
            {0.72375*scale, 0.003749999999999999, -0.28874999999999995
            }};
        for (int i = 0; i < exps.length; i++) {
            double[] acts = ell.values(new Coordinate(samples[i * 2], samples[i * 2 + 1]), null);
            assertArrayEquals(exps[i], acts, 1e-10);
        }
    }
}
