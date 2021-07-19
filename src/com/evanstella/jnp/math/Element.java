package com.evanstella.jnp.math;

import com.evanstella.jnp.core.IllegalDimensionException;
import com.evanstella.jnp.core.NDArray;
import com.evanstella.jnp.core.Numeric;

public final class Element {

    // no instances for you
    private Element () {}


    /**************************************************************************
     * <p>Calculates W^Z for each element in W and Z. For each element, if the
     * base is positive and real and the exponent is real, the function just
     * computes W^N using builtin Math.pow to avoid extra computation.
     * Otherwise the following is used:
     *
     * <p>For w = a+bi, z = c+di: w^z = e^(z*ln(r)+i*theta) for r = sqrt(a^2+b^2)
     * and theta = atan2(b,a) ...
     * <p>w^z = r^c * e^(-d*theta) *
     *          [ cos(d*ln(r) + c*theta) + i*sin(d*ln(r) + c*theta) ]
     *
     * <p>Dimensions must agree according to element wise operation rules.
     *
     * @param W         power base
     * @param Z         power exponent
     *
     * @return a Numeric with elements W^Z
     *************************************************************************/
    public static Numeric pow ( Numeric W, Numeric Z ) {
        validateDimensionsFatal( W, Z );
        double[] real1 = W.getDataReal();
        double[] imag1 = W.getDataImag();
        double[] real2 = Z.getDataReal();
        double[] imag2 = Z.getDataImag();
        double a,b,c,d,r,rc,theta,omega;
        boolean scalarW = W.isScalar(), scalarZ = Z.isScalar();

        a = real1[0]; c = real2[0];
        if ( imag1 == null ) b = 0.0;
        else b = imag1[0];
        if ( imag2 == null ) d = 0.0;
        else d = imag2[0];
        r = Math.sqrt( a*a + b*b );
        theta = Math.atan2(b, a);

        int[] newShape = ( real1.length > real2.length ) ? W.shape() : Z.shape();
        Numeric result = new Numeric( newShape );
        double[] resultReal = result.getDataReal();
        double[] resultImag = result.getDataImag();

        for ( int i = 0; i < resultReal.length; i++ ) {
            if ( !scalarW ) {
                a = real1[i];
                if ( imag1 == null ) b = 0.0;
                else b = imag1[i];
            }
            if ( !scalarZ ) {
                c = real2[i];
                if ( imag2 == null ) d = 0.0;
                else d = imag2[i];
            }
            // only do all the extra computation if we have to
            if ( a >= 0 && b == 0 && d == 0 ) {
                resultReal[i] = Math.pow( a, c );
                continue;
            }
            // we'll only get here if the result must be complex
            if ( resultImag == null ) resultImag = result.initializeDataImag();
            if ( !scalarW ) {
                r = Math.sqrt( a*a + b*b );
                theta = Math.atan2(b, a);
            }
            rc = Math.pow( r, c ) * Math.exp( -d * theta );
            omega = ( d * Math.log(r) ) + ( c * theta );
            resultReal[i] = rc * Math.cos(omega);
            resultImag[i] = rc * Math.sin(omega);
        }

        return result;
    }

    /*TODO*/
    private static void validateDimensionsFatal ( Numeric N1, Numeric N2 ) {
        if ( NDArray.dimensionsMatch(N1,N2) )
            return;
        if ( N1.isScalar() || N2.isScalar() )
            return;
        throw new IllegalDimensionException(
                "Element wise operation: dimensions must be equal or one Numeric must be scalar"
        );
    }

}
