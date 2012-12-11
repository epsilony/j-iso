/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package net.epsilony.levelset.ops;

import net.epsilony.levelset.functions.EllipseLSF;
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

        EllipseLSF ell = new EllipseLSF(x0, y0, a, b);
        ell.setDiffOrder(1);
        double[] samples = new double[]{x0, y0, x0 + a, y0, x0, y0 - b, 10, 4};
        double[][] exps = new double[][]{
            {1, 0, 0},
            {0, -0.0315066900254521, 0},
            {0, 0, 0.115524530093324 },
            { 0.651469127282252, 0.00130080814699714, -0.100162227318780}};
        for (int i = 0; i < exps.length; i++) {
            double[] acts = ell.values(new Coordinate(samples[i * 2], samples[i * 2 + 1]), null);
            assertArrayEquals(exps[i], acts, 1e-10);
        }
    }
}
