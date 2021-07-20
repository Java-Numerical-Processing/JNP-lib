package com.evanstella.jnp.math;

import com.evanstella.jnp.core.IllegalDimensionException;
import com.evanstella.jnp.core.Logical;
import com.evanstella.jnp.core.NDArray;
import com.evanstella.jnp.core.Numeric;

public final class Element {

    // no instances for you
    private Element () {}


    /**************************************************************************
     * <p>Take the negative of A element wise
     *
     * @param A the Numeric
     *
     * @return a Numeric with elements -A
     *************************************************************************/
    public static Numeric neg ( Numeric A ) {
        Numeric result = A.copy();
        double[] resultReal = result.getDataReal();
        double[] resultImag = result.getDataImag();

        for ( int i = 0; i < resultReal.length; i++ ) {
            resultReal[i] = -resultReal[i];
            if ( resultImag != null )
                resultImag[i] = -resultImag[i];
        }
        return result;
    }

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
     * <p>Divide A by B element-wise. If A or B is scalar, divide the
     * scalar by the elements of the other or vice-versa.
     *
     * @param A the first Numeric
     * @param B the second Numeric
     *
     * @return a Numeric with elements A/B
     *************************************************************************/
    public static Numeric div ( Numeric A, Numeric B ) {
        validateDimensionsFatal( A, B );
        double[] realA = A.getDataReal();
        double[] imagA = A.getDataImag();
        double[] realB = B.getDataReal();
        double[] imagB = B.getDataImag();
        double a,b,c,d,c2d2;
        boolean scalarA = A.isScalar(), scalarB = B.isScalar();
        a = realA[0]; c = realB[0];
        if ( imagA == null ) b = 0.0; else b = imagA[0];
        if ( imagB == null ) d = 0.0; else d = imagB[0];
        c2d2 = c*c + d*d;

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
                c2d2 = c*c + d*d;
            }
            // only do all the extra computation if we have to
            if ( b == 0 && d == 0 ) {
                resultReal[i] = a/c;
                continue;
            }
            // we'll only get here if the result must be complex
            if ( resultImag == null ) resultImag = result.initializeDataImag();
            resultReal[i] = (a*c + b*d) / c2d2;
            resultImag[i] = (b*c + a*d) / c2d2;
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

    /**************************************************************************
     * <p>Compute the natural log of A element-wise.
     *
     * Note for complex z = a+bi, ln(z) = ln(mag(z)) + phase(z)*i
     *
     * @param A Numeric of radians
     *
     * @return a Numeric with elements ln(A)
     *************************************************************************/
    public static Numeric log ( Numeric A ) {
        double[] realA = A.getDataReal();
        double[] imagA = A.getDataImag();
        double a,b,r;
        boolean scalarA = A.isScalar();

        Numeric result = new Numeric( A.shape() );
        double[] resultReal = result.getDataReal();
        double[] resultImag = result.getDataImag();

        for ( int i = 0; i < resultReal.length; i++ ) {
            a = realA[i];
            if ( imagA == null ) b = 0.0; else b = imagA[0];
            // only do all the extra computation if we have to
            if ( b == 0 && a >=0 ) {
                resultReal[i] = Math.log(a);
                continue;
            }
            // we'll only get here if the result must be complex
            if ( resultImag == null ) resultImag = result.initializeDataImag();
            r = Math.sqrt( a*a + b*b );
            resultReal[i] = Math.log(r);
            resultImag[i] = Math.atan2(b,a);
        }
        return result;
    }

    /**************************************************************************
     * <p>Compute the sine of A element-wise. A is in radians.
     *
     * @param A Numeric of radians
     *
     * @return a Numeric with elements sin(A)
     *************************************************************************/
    public static Numeric sin ( Numeric A ) {
        double[] realA = A.getDataReal();
        double[] imagA = A.getDataImag();
        double a,b;
        boolean scalarA = A.isScalar();

        Numeric result = new Numeric( A.shape() );
        double[] resultReal = result.getDataReal();
        double[] resultImag = result.getDataImag();

        for ( int i = 0; i < resultReal.length; i++ ) {
            a = realA[i];
            if ( imagA == null ) b = 0.0; else b = imagA[i];
            // only do all the extra computation if we have to
            if ( b == 0 ) {
                resultReal[i] = Math.sin(a);
                continue;
            }
            // we'll only get here if the result must be complex
            if ( resultImag == null ) resultImag = result.initializeDataImag();
            resultReal[i] = Math.sin(a) * Math.cosh(b);
            resultImag[i] = Math.cos(a) * Math.sinh(b);
        }
        return result;
    }

    /**************************************************************************
     * <p>Compute the sine of A element-wise. A is in degrees. Note that the
     * conversion from radians to degrees is not exact.
     *
     * @param A Numeric of radians
     *
     * @return a Numeric with elements sin(A)
     *************************************************************************/
    public static Numeric sind ( Numeric A ) {
        double[] realA = A.getDataReal();
        double[] imagA = A.getDataImag();
        double a,b, toRad;
        boolean scalarA = A.isScalar();

        toRad = 0.017453292519943295;

        Numeric result = new Numeric( A.shape() );
        double[] resultReal = result.getDataReal();
        double[] resultImag = result.getDataImag();

        for ( int i = 0; i < resultReal.length; i++ ) {
            a = realA[i] * toRad;
            if ( imagA == null ) b = 0.0; else b = imagA[i] * toRad;
            // only do all the extra computation if we have to
            if ( b == 0 ) {
                resultReal[i] = Math.sin(a);
                continue;
            }
            // we'll only get here if the result must be complex
            if ( resultImag == null ) resultImag = result.initializeDataImag();
            resultReal[i] = Math.sin(a) * Math.cosh(b);
            resultImag[i] = Math.cos(a) * Math.sinh(b);
        }
        return result;
    }

    /**************************************************************************
     * <p>Compute the cosine of A element-wise. A is in radians.
     *
     * @param A Numeric of radians
     *
     * @return a Numeric with elements cos(A)
     *************************************************************************/
    public static Numeric cos ( Numeric A ) {
        double[] realA = A.getDataReal();
        double[] imagA = A.getDataImag();
        double a,b;
        boolean scalarA = A.isScalar();

        Numeric result = new Numeric( A.shape() );
        double[] resultReal = result.getDataReal();
        double[] resultImag = result.getDataImag();

        for ( int i = 0; i < resultReal.length; i++ ) {
            a = realA[i];
            if ( imagA == null ) b = 0.0; else b = imagA[i];
            // only do all the extra computation if we have to
            if ( b == 0 ) {
                resultReal[i] = Math.cos(a);
                continue;
            }
            // we'll only get here if the result must be complex
            if ( resultImag == null ) resultImag = result.initializeDataImag();
            resultReal[i] = Math.cos(a) * Math.cosh(b);
            resultImag[i] = -1 * Math.sin(a) * Math.sinh(b);
        }
        return result;
    }

    /**************************************************************************
     * <p>Compute the cosine of A element-wise. A is in degrees. Note that the
     * conversion from radians to degrees is not exact
     *
     * @param A Numeric of radians
     *
     * @return a Numeric with elements cos(A)
     *************************************************************************/
    public static Numeric cosd ( Numeric A ) {
        double[] realA = A.getDataReal();
        double[] imagA = A.getDataImag();
        double a,b, toRad;
        boolean scalarA = A.isScalar();

        toRad = 0.017453292519943295;

        Numeric result = new Numeric( A.shape() );
        double[] resultReal = result.getDataReal();
        double[] resultImag = result.getDataImag();

        for ( int i = 0; i < resultReal.length; i++ ) {
            a = realA[i] * toRad;
            if ( imagA == null ) b = 0.0; else b = imagA[i] * toRad;
            // only do all the extra computation if we have to
            if ( b == 0 ) {
                resultReal[i] = Math.cos(a);
                continue;
            }
            // we'll only get here if the result must be complex
            if ( resultImag == null ) resultImag = result.initializeDataImag();
            resultReal[i] = Math.cos(a) * Math.cosh(b);
            resultImag[i] = -1 * Math.sin(a) * Math.sinh(b);
        }
        return result;
    }

    /**************************************************************************
     * <p>Compute the tangent of A element-wise. A is in radians.
     *
     * @param A Numeric of radians
     *
     * @return a Numeric with elements tan(A)
     *************************************************************************/
    public static Numeric tan ( Numeric A ) {
        double[] realA = A.getDataReal();
        double[] imagA = A.getDataImag();
        double a,b,cos2acosh2b;
        boolean scalarA = A.isScalar();

        Numeric result = new Numeric( A.shape() );
        double[] resultReal = result.getDataReal();
        double[] resultImag = result.getDataImag();

        for ( int i = 0; i < resultReal.length; i++ ) {
            a = realA[i];
            if ( imagA == null ) b = 0.0; else b = imagA[i] * 2;
            // only do all the extra computation if we have to
            if ( b == 0 ) {
                resultReal[i] = Math.tan(a);
                continue;
            }
            // we'll only get here if the result must be complex
            if ( resultImag == null ) resultImag = result.initializeDataImag();
            a *= 2;  // note that b has already been x2
            cos2acosh2b = Math.cos(a) + Math.cosh(b);
            resultReal[i] = Math.sin(a)  / cos2acosh2b;
            resultImag[i] = Math.sinh(b) / cos2acosh2b;
        }
        return result;
    }

    /**************************************************************************
     * <p>Compute the tangent of A element-wise. A is in degrees, Not that the
     * conversion to radians is not exact.
     *
     * @param A Numeric of radians
     *
     * @return a Numeric with elements tan(A)
     *************************************************************************/
    public static Numeric tand ( Numeric A ) {
        double[] realA = A.getDataReal();
        double[] imagA = A.getDataImag();
        double a,b,cos2acosh2b,toRad;
        boolean scalarA = A.isScalar();

        toRad = 0.017453292519943295;

        Numeric result = new Numeric( A.shape() );
        double[] resultReal = result.getDataReal();
        double[] resultImag = result.getDataImag();

        for ( int i = 0; i < resultReal.length; i++ ) {
            a = realA[i] * toRad;
            if ( imagA == null ) b = 0.0; else b = imagA[i] * 2 * toRad;
            // only do all the extra computation if we have to
            if ( b == 0 ) {
                resultReal[i] = Math.tan(a);
                continue;
            }
            // we'll only get here if the result must be complex
            if ( resultImag == null ) resultImag = result.initializeDataImag();
            a *= 2;  // note that b has already been x2
            cos2acosh2b = Math.cos(a) + Math.cosh(b);
            resultReal[i] = Math.sin(a)  / cos2acosh2b;
            resultImag[i] = Math.sinh(b) / cos2acosh2b;
        }
        return result;
    }

    /**************************************************************************
     * <p>Compute the hyperbolic cosine of A element-wise. A is in radians.
     *
     * @param A Numeric of radians
     *
     * @return a Numeric with elements cosh(A)
     *************************************************************************/
    public static Numeric cosh ( Numeric A ) {
        double[] realA = A.getDataReal();
        double[] imagA = A.getDataImag();
        double a,b;
        boolean scalarA = A.isScalar();

        Numeric result = new Numeric( A.shape() );
        double[] resultReal = result.getDataReal();
        double[] resultImag = result.getDataImag();

        for ( int i = 0; i < resultReal.length; i++ ) {
            a = realA[i];
            if ( imagA == null ) b = 0.0; else b = imagA[i];
            // only do all the extra computation if we have to
            if ( b == 0 ) {
                resultReal[i] = Math.cosh(a);
                continue;
            }
            // we'll only get here if the result must be complex
            if ( resultImag == null ) resultImag = result.initializeDataImag();
            resultReal[i] = Math.cosh(a) * Math.cos(b);
            resultImag[i] = Math.sinh(a) * Math.sin(b);
        }
        return result;
    }

    /**************************************************************************
     * <p>Compute the hyperbolic sine of A element-wise. A is in radians.
     *
     * @param A Numeric of radians
     *
     * @return a Numeric with elements sinh(A)
     *************************************************************************/
    public static Numeric sinh ( Numeric A ) {
        double[] realA = A.getDataReal();
        double[] imagA = A.getDataImag();
        double a,b;
        boolean scalarA = A.isScalar();

        Numeric result = new Numeric( A.shape() );
        double[] resultReal = result.getDataReal();
        double[] resultImag = result.getDataImag();

        for ( int i = 0; i < resultReal.length; i++ ) {
            a = realA[i];
            if ( imagA == null ) b = 0.0; else b = imagA[i];
            // only do all the extra computation if we have to
            if ( b == 0 ) {
                resultReal[i] = Math.sinh(a);
                continue;
            }
            // we'll only get here if the result must be complex
            if ( resultImag == null ) resultImag = result.initializeDataImag();
            resultReal[i] = Math.sinh(a) * Math.cos(b);
            resultImag[i] = Math.cosh(a) * Math.sin(b);
        }
        return result;
    }

    /**************************************************************************
     * <p>Compute the hyperbolic tangent of A element-wise. A is in radians.
     *
     * @param A Numeric of radians
     *
     * @return a Numeric with elements tanh(A)
     *************************************************************************/
    public static Numeric tanh ( Numeric A ) {
        double[] realA = A.getDataReal();
        double[] imagA = A.getDataImag();
        double a,b,cosh2acos2b;
        boolean scalarA = A.isScalar();

        Numeric result = new Numeric( A.shape() );
        double[] resultReal = result.getDataReal();
        double[] resultImag = result.getDataImag();

        for ( int i = 0; i < resultReal.length; i++ ) {
            a = realA[i];
            if ( imagA == null ) b = 0.0; else b = imagA[i] * 2;
            // only do all the extra computation if we have to
            if ( b == 0 ) {
                resultReal[i] = Math.tanh(a);
                continue;
            }
            // we'll only get here if the result must be complex
            if ( resultImag == null ) resultImag = result.initializeDataImag();
            a *= 2;  // note that b has already been x2
            cosh2acos2b = Math.cosh(a) + Math.cos(b);
            resultReal[i] = Math.sinh(a) / cosh2acos2b;
            resultImag[i] = Math.sin(b)  / cosh2acos2b;
        }
        return result;
    }

    /**************************************************************************
     * <p>Compare A and B element wise. If one of the values is scalar, compare
     * that value to the other numeric. Only compares real components as
     * comparison is not well defined for complex values.
     *
     * @param A the first Numeric
     * @param B the second Numeric
     *
     * @return a Logical that indexes A > B
     *************************************************************************/
    public static Logical gre ( Numeric A, Numeric B ) {
        validateDimensionsFatal( A, B );
        double[] realA = A.getDataReal();
        double[] realB = B.getDataReal();
        double a,b;
        boolean scalarA = A.isScalar(), scalarB = B.isScalar();
        a = realA[0]; b = realB[0];

        int[] newShape = ( realA.length > realB.length ) ? A.shape() : B.shape();
        Logical result = new Logical( newShape );
        boolean[] resultData = result.getData();

        for ( int i = 0; i < resultData.length; i++ ) {
            if ( !scalarA ) a = realA[i];
            if ( !scalarB ) b = realB[i];
            resultData[i] = a > b;
        }
        return result;
    }

    /**************************************************************************
     * <p>Compare A and B element wise. If one of the values is scalar, compare
     * that value to the other numeric. Only compares real components as
     * comparison is not well defined for complex values.
     *
     * @param A the first Numeric
     * @param B the second Numeric
     *
     * @return a Logical that indexes A >= B
     *************************************************************************/
    public static Logical greq ( Numeric A, Numeric B ) {
        validateDimensionsFatal( A, B );
        double[] realA = A.getDataReal();
        double[] realB = B.getDataReal();
        double a,b;
        boolean scalarA = A.isScalar(), scalarB = B.isScalar();
        a = realA[0]; b = realB[0];

        int[] newShape = ( realA.length > realB.length ) ? A.shape() : B.shape();
        Logical result = new Logical( newShape );
        boolean[] resultData = result.getData();

        for ( int i = 0; i < resultData.length; i++ ) {
            if ( !scalarA ) a = realA[i];
            if ( !scalarB ) b = realB[i];
            resultData[i] = a >= b;
        }
        return result;
    }

    /**************************************************************************
     * <p>Compare A and B element wise. If one of the values is scalar, compare
     * that value to the other numeric. Only compares real components as
     * comparison is not well defined for complex values.
     *
     * @param A the first Numeric
     * @param B the second Numeric
     *
     * @return a Logical that indexes A < B
     *************************************************************************/
    public static Logical less ( Numeric A, Numeric B ) {
        validateDimensionsFatal( A, B );
        double[] realA = A.getDataReal();
        double[] realB = B.getDataReal();
        double a,b;
        boolean scalarA = A.isScalar(), scalarB = B.isScalar();
        a = realA[0]; b = realB[0];

        int[] newShape = ( realA.length > realB.length ) ? A.shape() : B.shape();
        Logical result = new Logical( newShape );
        boolean[] resultData = result.getData();

        for ( int i = 0; i < resultData.length; i++ ) {
            if ( !scalarA ) a = realA[i];
            if ( !scalarB ) b = realB[i];
            resultData[i] = a < b;
        }
        return result;
    }

    /**************************************************************************
     * <p>Compare A and B element wise. If one of the values is scalar, compare
     * that value to the other numeric. Only compares real components as
     * comparison is not well defined for complex values.
     *
     * @param A the first Numeric
     * @param B the second Numeric
     *
     * @return a Logical that indexes A <= B
     *************************************************************************/
    public static Logical leq ( Numeric A, Numeric B ) {
        validateDimensionsFatal(A, B);
        double[] realA = A.getDataReal();
        double[] realB = B.getDataReal();
        double a, b;
        boolean scalarA = A.isScalar(), scalarB = B.isScalar();
        a = realA[0];
        b = realB[0];

        int[] newShape = (realA.length > realB.length) ? A.shape() : B.shape();
        Logical result = new Logical(newShape);
        boolean[] resultData = result.getData();

        for (int i = 0; i < resultData.length; i++) {
            if (!scalarA) a = realA[i];
            if (!scalarB) b = realB[i];
            resultData[i] = a <= b;
        }
        return result;
    }

    /**************************************************************************
     * <p>Compare A and B element wise for equality with a tolerance. If one
     * of the values is scalar, compare that value with the elements of the
     * other
     *
     * @param A the first Numeric
     * @param B the second Numeric
     *
     * @return a Logical that indexes A == B
     *************************************************************************/
    public static Logical equal ( Numeric A, Numeric B, double tolerance ) {
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
        Logical result = new Logical( newShape );
        boolean[] resultData = result.getData();

        for ( int i = 0; i < resultData.length; i++ ) {
            if ( !scalarA ) {
                a = realA[i];
                if ( imagA == null ) b = 0.0; else b = imagA[i];
            }
            if ( !scalarB ) {
                c = realB[i];
                if ( imagB == null ) d = 0.0; else d = imagB[i];
            }
            resultData[i] =
                (Math.abs(a - c) <= tolerance) && (Math.abs(b - d) <= tolerance);
        }
        return result;
    }

    /**************************************************************************
     * <p>Take !A for each element in the Logical A
     *
     * @param A the Logical
     *
     * @return a Logical with elements !A
     *************************************************************************/
    public static Logical not ( Logical A ) {
        Logical result = A.copy();
        boolean[] resultData = result.getData();

        for ( int i = 0; i < resultData.length; i++ )
            resultData[i] = !resultData[i];

        return result;
    }

    /**************************************************************************
     * <p>Take A && B element wise
     *
     * @param A the first Numeric
     * @param B the second Numeric
     *
     * @return a Logical with elements A && b
     *************************************************************************/
    public static Logical and ( Logical A, Logical B ) {
        validateDimensionsFatal( A, B );
        boolean[] dataA = A.getData();
        boolean[] dataB = B.getData();
        boolean a,b;
        boolean oneDA = dataA.length == 1, oneDB = dataB.length == 1;
        a = dataA[0]; b = dataA[0];

        int[] newShape = ( dataA.length > dataB.length ) ? A.shape() : B.shape();
        Logical result = new Logical( newShape );
        boolean[] resultData = result.getData();

        for ( int i = 0; i < resultData.length; i++ ) {
            if ( !oneDA )
                a = dataA[i];
            if ( !oneDB )
                b = dataB[i];
            resultData[i] = a && b;
        }
        return result;
    }

    /**************************************************************************
     * <p>Take A || B element wise
     *
     * @param A the first Numeric
     * @param B the second Numeric
     *
     * @return a Logical with elements A || b
     *************************************************************************/
    public static Logical or ( Logical A, Logical B ) {
        validateDimensionsFatal( A, B );
        boolean[] dataA = A.getData();
        boolean[] dataB = B.getData();
        boolean a,b;
        boolean oneDA = dataA.length == 1, oneDB = dataB.length == 1;
        a = dataA[0]; b = dataA[0];

        int[] newShape = ( dataA.length > dataB.length ) ? A.shape() : B.shape();
        Logical result = new Logical( newShape );
        boolean[] resultData = result.getData();

        for ( int i = 0; i < resultData.length; i++ ) {
            if ( !oneDA )
                a = dataA[i];
            if ( !oneDB )
                b = dataB[i];
            resultData[i] = a || b;
        }
        return result;
    }


    /*Ensure that either one Numeric is scalar or both have the same shape*/
    private static void validateDimensionsFatal ( Numeric N1, Numeric N2 ) {
        if ( NDArray.dimensionsMatch(N1,N2) )
            return;
        if ( N1.isScalar() || N2.isScalar() )
            return;
        throw new IllegalDimensionException(
            "Element wise operation: dimensions must be equal or one Numeric must be scalar"
        );
    }

    /*Ensure that either one Logical is 1x1 or both have the same shape*/
    private static void validateDimensionsFatal ( Logical N1, Logical N2 ) {
        if ( NDArray.dimensionsMatch(N1,N2) )
            return;
        if ( N1.getData().length == 1 || N1.getData().length == 1 )
            return;
        throw new IllegalDimensionException(
            "Element wise operation: dimensions must be equal or one Logical must be 1x1"
        );
    }



}
