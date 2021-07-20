package com.evanstella.jnp.math;

import com.evanstella.jnp.core.IllegalDimensionException;
import com.evanstella.jnp.core.NDArray;
import com.evanstella.jnp.core.Numeric;

public final class Element {

    // no instances for you
    private Element () {}


    /**************************************************************************
     * <p>Add A and B element-wise. If A or B is scalar, add the scalar to the
     * elements of the other.
     *
     * @param A the first Numeric
     * @param B the second Numeric
     *
     * @return a Numeric with elements A+B
     *************************************************************************/
    public static Numeric add ( Numeric A, Numeric B ) {
        validateDimensionsFatal( A, B );
        double[] realA = A.getDataReal();
        double[] imagA = A.getDataImag();
        double[] realB = B.getDataReal();
        double[] imagB = B.getDataImag();
        double a,b,c,d;
        boolean scalarA = A.isScalar(), scalarB = B.isScalar();
        a = realA[0]; c = realB[0];
        if ( imagA == null ) b = 0.0; else b = imagA[0];
        if ( imagB == null ) d = 0.0; else d = imagB[0];

        int[] newShape = ( realA.length > realB.length ) ? A.shape() : B.shape();
        Numeric result = new Numeric( newShape );
        double[] resultReal = result.getDataReal();
        double[] resultImag = result.getDataImag();

        for ( int i = 0; i < resultReal.length; i++ ) {
            if ( !scalarA ) {
                a = realA[i];
                if ( imagA == null ) b = 0.0; else b = imagA[i];
            }
            if ( !scalarB ) {
                c = realB[i];
                if ( imagB == null ) d = 0.0; else d = imagB[i];
            }
            resultReal[i] = a + c;
            // only do all the extra computation if we have to
            if ( b == 0 && d == 0 )
                continue;
            // we'll only get here if the result must be complex
            if ( resultImag == null ) resultImag = result.initializeDataImag();
            resultImag[i] = b + d;
        }
        return result;
    }

    /**************************************************************************
     * <p>Subtract A and B element-wise. If A or B is scalar, subtract the
     * scalar from the elements of the other.
     *
     * @param A the first Numeric
     * @param B the second Numeric
     *
     * @return a Numeric with elements A-B
     *************************************************************************/
    public static Numeric sub ( Numeric A, Numeric B ) {
        validateDimensionsFatal( A, B );
        double[] realA = A.getDataReal();
        double[] imagA = A.getDataImag();
        double[] realB = B.getDataReal();
        double[] imagB = B.getDataImag();
        double a,b,c,d;
        boolean scalarA = A.isScalar(), scalarB = B.isScalar();
        a = realA[0]; c = realB[0];
        if ( imagA == null ) b = 0.0; else b = imagA[0];
        if ( imagB == null ) d = 0.0; else d = imagB[0];

        int[] newShape = ( realA.length > realB.length ) ? A.shape() : B.shape();
        Numeric result = new Numeric( newShape );
        double[] resultReal = result.getDataReal();
        double[] resultImag = result.getDataImag();

        for ( int i = 0; i < resultReal.length; i++ ) {
            if ( !scalarA ) {
                a = realA[i];
                if ( imagA == null ) b = 0.0; else b = imagA[i];
            }
            if ( !scalarB ) {
                c = realB[i];
                if ( imagB == null ) d = 0.0; else d = imagB[i];
            }
            resultReal[i] = a - c;
            // only do all the extra computation if we have to
            if ( b == 0 && d == 0 )
                continue;
            // we'll only get here if the result must be complex
            if ( resultImag == null ) resultImag = result.initializeDataImag();
            resultImag[i] = b - d;
        }
        return result;
    }

    /**************************************************************************
     * <p>Multiply A and B element-wise. If A or B is scalar, multiply the
     * scalar with the elements of the other.
     *
     * @param A the first Numeric
     * @param B the second Numeric
     *
     * @return a Numeric with elements A*B
     *************************************************************************/
    public static Numeric mul ( Numeric A, Numeric B ) {
        validateDimensionsFatal( A, B );
        double[] realA = A.getDataReal();
        double[] imagA = A.getDataImag();
        double[] realB = B.getDataReal();
        double[] imagB = B.getDataImag();
        double a,b,c,d;
        boolean scalarA = A.isScalar(), scalarB = B.isScalar();
        a = realA[0]; c = realB[0];
        if ( imagA == null ) b = 0.0; else b = imagA[0];
        if ( imagB == null ) d = 0.0; else d = imagB[0];

        int[] newShape = ( realA.length > realB.length ) ? A.shape() : B.shape();
        Numeric result = new Numeric( newShape );
        double[] resultReal = result.getDataReal();
        double[] resultImag = result.getDataImag();

        for ( int i = 0; i < resultReal.length; i++ ) {
            if ( !scalarA ) {
                a = realA[i];
                if ( imagA == null ) b = 0.0; else b = imagA[i];
            }
            if ( !scalarB ) {
                c = realB[i];
                if ( imagB == null ) d = 0.0; else d = imagB[i];
            }
            resultReal[i] = a*c - b*d;
            // only do all the extra computation if we have to
            if ( b == 0 && d == 0 )
                continue;
            // we'll only get here if the result must be complex
            if ( resultImag == null ) resultImag = result.initializeDataImag();
            resultImag[i] = a*d + b*c;
        }
        return result;
    }

    /**************************************************************************
     * <p>Multiply A and B element-wise. If A or B is scalar, multiply the
     * scalar with the elements of the other.
     *
     * @param A the first Numeric
     * @param B the second Numeric
     *
     * @return a Numeric with elements A*B
     *************************************************************************/
    public static Numeric mul ( Numeric A, Numeric B ) {
        validateDimensionsFatal( A, B );
        double[] realA = A.getDataReal();
        double[] imagA = A.getDataImag();
        double[] realB = B.getDataReal();
        double[] imagB = B.getDataImag();
        double a,b,c,d;
        boolean scalarA = A.isScalar(), scalarB = B.isScalar();
        a = realA[0]; c = realB[0];
        if ( imagA == null ) b = 0.0; else b = imagA[0];
        if ( imagB == null ) d = 0.0; else d = imagB[0];

        int[] newShape = ( realA.length > realB.length ) ? A.shape() : B.shape();
        Numeric result = new Numeric( newShape );
        double[] resultReal = result.getDataReal();
        double[] resultImag = result.getDataImag();

        for ( int i = 0; i < resultReal.length; i++ ) {
            if ( !scalarA ) {
                a = realA[i];
                if ( imagA == null ) b = 0.0; else b = imagA[i];
            }
            if ( !scalarB ) {
                c = realB[i];
                if ( imagB == null ) d = 0.0; else d = imagB[i];
            }
            resultReal[i] = a*c - b*d;
            // only do all the extra computation if we have to
            if ( b == 0 && d == 0 )
                continue;
            // we'll only get here if the result must be complex
            if ( resultImag == null ) resultImag = result.initializeDataImag();
            resultImag[i] = a*d + b*c;
        }
        return result;
    }

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
        if ( imag1 == null ) b = 0.0; else b = imag1[0];
        if ( imag2 == null ) d = 0.0; else d = imag2[0];
        r = Math.sqrt( a*a + b*b );
        theta = Math.atan2(b, a);

        int[] newShape = ( real1.length > real2.length ) ? W.shape() : Z.shape();
        Numeric result = new Numeric( newShape );
        double[] resultReal = result.getDataReal();
        double[] resultImag = result.getDataImag();

        for ( int i = 0; i < resultReal.length; i++ ) {
            if ( !scalarW ) {
                a = real1[i];
                if ( imag1 == null ) b = 0.0; else b = imag1[i];
            }
            if ( !scalarZ ) {
                c = real2[i];
                if ( imag2 == null ) d = 0.0; else d = imag2[i];
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
                theta = Math.atan2( b, a );
            }
            rc = Math.pow( r, c ) * Math.exp( -d * theta );
            omega = ( d * Math.log(r) ) + ( c * theta );
            resultReal[i] = rc * Math.cos( omega );
            resultImag[i] = rc * Math.sin( omega );
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
