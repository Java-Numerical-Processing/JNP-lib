/*
 * MIT License
 *
 * Copyright (c) 2021 Evan Stella
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject to
 * the following conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
 * LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
 * OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
 * WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 *
 */

package com.evanstella.jnp.math;

import com.evanstella.jnp.core.*;

/******************************************************************************
 * Element encapsulates all of the element-wise operations that can be done on
 * NDArrays.
 *
 * @author Evan Stella
 *****************************************************************************/
public final class Element {

    // no instances for you
    private Element () {}


    /**************************************************************************
     * <p>Take the negative of A element wise
     *
     * @param A the Complex
     *
     * @return a Complex with elements -A
     *************************************************************************/
    public static Complex neg ( Complex A ) {
        Complex result = A.copy();
        double[] resultReal = result.getDataReal();
        double[] resultImag = result.getDataImag();

        for ( int i = 0; i < resultReal.length; i++ ) {
            resultReal[i] = -resultReal[i];
            resultImag[i] = -resultImag[i];
        }
        return result;
    }

    /**************************************************************************
     * <p>Take the negative of A element wise
     *
     * @param A the Complex
     *
     * @return a Complex with elements -A
     *************************************************************************/
    public static Numeric neg ( Numeric A ) {
        Numeric result = A.copy();
        double[] resultData = result.getData();

        for ( int i = 0; i < resultData.length; i++ ) {
            resultData[i] = -resultData[i];
        }
        return result;
    }

    /**************************************************************************
     * <p>Calculate the phase of A element wise
     *
     * @param A the Complex
     *
     * @return a Complex with elements phase(A)
     *************************************************************************/
    public static Complex phase ( Complex A ) {
        Complex result = new Complex( A.shape() );
        double[] resultReal = result.getDataReal();
        double[] dataReal = A.getDataReal();
        double[] dataImag = A.getDataImag();

        double a,b;
        for ( int i = 0; i < resultReal.length; i++ ) {
            a = dataReal[i];
            b = dataImag[i];
            resultReal[i] = Math.atan2( b, a );
        }
        return result;
    }

    /**************************************************************************
     * <p>Calculate the magnitude of A element wise
     *
     * @param A the Complex
     *
     * @return a Complex with elements magnitude(A)
     *************************************************************************/
    public static Complex mag ( Complex A ) {
        Complex result = new Complex( A.shape() );
        double[] resultReal = result.getDataReal();
        double[] dataReal = A.getDataReal();
        double[] dataImag = A.getDataImag();

        double a,b;
        for ( int i = 0; i < resultReal.length; i++ ) {
            a = dataReal[i];
            b = dataImag[i];
            resultReal[i] = Math.sqrt( a*a + b*b );
        }
        return result;
    }

    /**************************************************************************
     * <p>Add A and B element-wise. If A or B is scalar, add the scalar to the
     * elements of the other.
     *
     * @param A the first Complex
     * @param B the second Complex
     *
     * @return a Complex with elements A+B
     *************************************************************************/
    public static Complex add ( Complex A, Complex B ) {
        validateDimensionsFatal( A, B );
        double[] realA = A.getDataReal();
        double[] imagA = A.getDataImag();
        double[] realB = B.getDataReal();
        double[] imagB = B.getDataImag();
        double a,b,c,d;
        boolean scalarA = A.isScalar(), scalarB = B.isScalar();
        a = realA[0]; c = realB[0];
        b = imagA[0]; d = imagB[0];

        int[] newShape = ( realA.length > realB.length ) ? A.shape() : B.shape();
        Complex result = new Complex( newShape );
        double[] resultReal = result.getDataReal();
        double[] resultImag = result.getDataImag();

        for ( int i = 0; i < resultReal.length; i++ ) {
            if ( !scalarA ) {
                a = realA[i];
                b = imagA[i];
            }
            if ( !scalarB ) {
                c = realB[i];
                d = imagB[i];
            }
            resultReal[i] = a + c;
            resultImag[i] = b + d;
        }
        return result;
    }

    /**************************************************************************
     * <p>Add A and B element-wise. If A or B is scalar, add the scalar to the
     * elements of the other.
     *
     * @param A the first Complex
     * @param B the second Complex
     *
     * @return a Complex with elements A+B
     *************************************************************************/
    public static Numeric add ( Numeric A, Numeric B ) {
        validateDimensionsFatal( A, B );
        double[] realA = A.getData();
        double[] realB = B.getData();
        double a,b;
        boolean scalarA = A.isScalar(), scalarB = B.isScalar();
        a = realA[0]; b = realB[0];

        int[] newShape = ( realA.length > realB.length ) ? A.shape() : B.shape();
        Numeric result = new Numeric( newShape );
        double[] resultReal = result.getData();

        for ( int i = 0; i < resultReal.length; i++ ) {
            if ( !scalarA ) a = realA[i];
            if ( !scalarB ) b = realB[i];
            resultReal[i] = a + b;
        }
        return result;
    }

    /**************************************************************************
     * <p>Subtract B from A element-wise. If A or B is scalar, subtract the
     * scalar from the elements of the other.
     *
     * @param A the first Complex
     * @param B the second Complex
     *
     * @return a Complex with elements A-B
     *************************************************************************/
    public static Complex sub ( Complex A, Complex B ) {
        validateDimensionsFatal( A, B );
        double[] realA = A.getDataReal();
        double[] imagA = A.getDataImag();
        double[] realB = B.getDataReal();
        double[] imagB = B.getDataImag();
        double a,b,c,d;
        boolean scalarA = A.isScalar(), scalarB = B.isScalar();
        a = realA[0]; c = realB[0];
        b = imagA[0]; d = imagB[0];

        int[] newShape = ( realA.length > realB.length ) ? A.shape() : B.shape();
        Complex result = new Complex( newShape );
        double[] resultReal = result.getDataReal();
        double[] resultImag = result.getDataImag();

        for ( int i = 0; i < resultReal.length; i++ ) {
            if ( !scalarA ) {
                a = realA[i];
                b = imagA[i];
            }
            if ( !scalarB ) {
                c = realB[i];
                d = imagB[i];
            }
            resultReal[i] = a - c;
            resultImag[i] = b - d;
        }
        return result;
    }

    /**************************************************************************
     * <p>Subtract B from A element-wise. If A or B is scalar, subtract the
     * scalar from the elements of the other.
     *
     * @param A the first Complex
     * @param B the second Complex
     *
     * @return a Complex with elements A-B
     *************************************************************************/
    public static Numeric sub ( Numeric A, Numeric B ) {
        validateDimensionsFatal( A, B );
        double[] realA = A.getData();
        double[] realB = B.getData();
        double a,b;
        boolean scalarA = A.isScalar(), scalarB = B.isScalar();
        a = realA[0]; b = realB[0];

        int[] newShape = ( realA.length > realB.length ) ? A.shape() : B.shape();
        Numeric result = new Numeric( newShape );
        double[] resultReal = result.getData();

        for ( int i = 0; i < resultReal.length; i++ ) {
            if ( !scalarA ) a = realA[i];
            if ( !scalarB ) b = realB[i];
            resultReal[i] = a - b;
        }
        return result;
    }

    /**************************************************************************
     * <p>Multiply A and B element-wise. If A or B is scalar, multiply the
     * scalar with the elements of the other.
     *
     * @param A the first Complex
     * @param B the second Complex
     *
     * @return a Complex with elements A*B
     *************************************************************************/
    public static Complex mul ( Complex A, Complex B ) {
        validateDimensionsFatal( A, B );
        double[] realA = A.getDataReal();
        double[] imagA = A.getDataImag();
        double[] realB = B.getDataReal();
        double[] imagB = B.getDataImag();
        double a,b,c,d;
        boolean scalarA = A.isScalar(), scalarB = B.isScalar();
        a = realA[0]; c = realB[0];
        b = imagA[0]; d = imagB[0];

        int[] newShape = ( realA.length > realB.length ) ? A.shape() : B.shape();
        Complex result = new Complex( newShape );
        double[] resultReal = result.getDataReal();
        double[] resultImag = result.getDataImag();

        for ( int i = 0; i < resultReal.length; i++ ) {
            if ( !scalarA ) {
                a = realA[i];
                b = imagA[i];
            }
            if ( !scalarB ) {
                c = realB[i];
                d = imagB[i];
            }
            resultReal[i] = a*c - b*d;
            resultImag[i] = a*d + b*c;
        }
        return result;
    }

    /**************************************************************************
     * <p>Multiply A and B element-wise. If A or B is scalar, multiply the
     * scalar with the elements of the other.
     *
     * @param A the first Complex
     * @param B the second Complex
     *
     * @return a Complex with elements A*B
     *************************************************************************/
    public static Numeric mul ( Numeric A, Numeric B ) {
        validateDimensionsFatal( A, B );
        double[] realA = A.getData();
        double[] realB = B.getData();
        double a,b;
        boolean scalarA = A.isScalar(), scalarB = B.isScalar();
        a = realA[0]; b = realB[0];

        int[] newShape = ( realA.length > realB.length ) ? A.shape() : B.shape();
        Numeric result = new Numeric( newShape );
        double[] resultReal = result.getData();

        for ( int i = 0; i < resultReal.length; i++ ) {
            if ( !scalarA ) a = realA[i];
            if ( !scalarB ) b = realB[i];
            resultReal[i] = a*b;
        }
        return result;
    }

    /**************************************************************************
     * <p>Divide A by B element-wise. If A or B is scalar, divide the
     * scalar by the elements of the other or vice-versa.
     *
     * @param A the first Complex
     * @param B the second Complex
     *
     * @return a Complex with elements A/B
     *************************************************************************/
    public static Complex div ( Complex A, Complex B ) {
        validateDimensionsFatal( A, B );
        double[] realA = A.getDataReal();
        double[] imagA = A.getDataImag();
        double[] realB = B.getDataReal();
        double[] imagB = B.getDataImag();
        double a,b,c,d,c2d2;
        boolean scalarA = A.isScalar(), scalarB = B.isScalar();
        a = realA[0]; c = realB[0];
        b = imagA[0]; d = imagB[0];
        c2d2 = c*c + d*d;

        int[] newShape = ( realA.length > realB.length ) ? A.shape() : B.shape();
        Complex result = new Complex( newShape );
        double[] resultReal = result.getDataReal();
        double[] resultImag = result.getDataImag();

        for ( int i = 0; i < resultReal.length; i++ ) {
            if ( !scalarA ) {
                a = realA[i];
                b = imagA[i];
            }
            if ( !scalarB ) {
                c = realB[i];
                d = imagB[i];
                c2d2 = c*c + d*d;
            }
            resultReal[i] = (a*c + b*d) / c2d2;
            resultImag[i] = (b*c + a*d) / c2d2;
        }
        return result;
    }

    /**************************************************************************
     * <p>Divide A by B element-wise. If A or B is scalar, divide the
     * scalar by the elements of the other or vice-versa.
     *
     * @param A the first Complex
     * @param B the second Complex
     *
     * @return a Complex with elements A/B
     *************************************************************************/
    public static Numeric div ( Numeric A, Numeric B ) {
        validateDimensionsFatal( A, B );
        double[] realA = A.getData();
        double[] realB = B.getData();
        double a,b;
        boolean scalarA = A.isScalar(), scalarB = B.isScalar();
        a = realA[0]; b = realB[0];

        int[] newShape = ( realA.length > realB.length ) ? A.shape() : B.shape();
        Numeric result = new Numeric( newShape );
        double[] resultReal = result.getData();

        for ( int i = 0; i < resultReal.length; i++ ) {
            if (!scalarA) a = realA[i];
            if (!scalarB) b = realB[i];
            resultReal[i] = a / b;
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
     * @return a Complex with elements W^Z
     *************************************************************************/
    public static Complex pow ( Complex W, Complex Z ) {
        validateDimensionsFatal( W, Z );
        double[] real1 = W.getDataReal();
        double[] imag1 = W.getDataImag();
        double[] real2 = Z.getDataReal();
        double[] imag2 = Z.getDataImag();
        double a,b,c,d,r,rc,theta,omega;
        boolean scalarW = W.isScalar(), scalarZ = Z.isScalar();

        a = real1[0]; c = real2[0];
        b = imag1[0]; d = imag2[0];
        r = Math.sqrt( a*a + b*b );
        theta = Math.atan2(b, a);

        int[] newShape = ( real1.length > real2.length ) ? W.shape() : Z.shape();
        Complex result = new Complex( newShape );
        double[] resultReal = result.getDataReal();
        double[] resultImag = result.getDataImag();

        for ( int i = 0; i < resultReal.length; i++ ) {
            if ( !scalarW ) {
                a = real1[i];
                b = imag1[i];
                r = Math.sqrt( a*a + b*b );
                theta = Math.atan2( b, a );
            }
            if ( !scalarZ ) {
                c = real2[i];
                d = imag2[i];
            }
            rc = Math.pow( r, c ) * Math.exp( -d * theta );
            omega = ( d * Math.log(r) ) + ( c * theta );
            resultReal[i] = rc * Math.cos( omega );
            resultImag[i] = rc * Math.sin( omega );
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
     * @return a Complex with elements W^Z
     *************************************************************************/
    public static Numeric pow ( Numeric W, Numeric Z ) {
        validateDimensionsFatal( W, Z );
        double[] real1 = W.getData();
        double[] real2 = Z.getData();
        double a,b;
        boolean scalarW = W.isScalar(), scalarZ = Z.isScalar();

        a = real1[0]; b = real2[0];

        int[] newShape = ( real1.length > real2.length ) ? W.shape() : Z.shape();
        Numeric result = new Numeric( newShape );
        double[] resultReal = result.getData();

        for ( int i = 0; i < resultReal.length; i++ ) {
            if ( !scalarW ) a = real1[i];
            if ( !scalarZ ) b = real2[i];
            resultReal[i] = Math.pow( a, b );
        }
        return result;
    }

    /**************************************************************************
     * <p>Compute the natural log of A element-wise.
     *
     * Note for complex z = a+bi, ln(z) = ln(mag(z)) + phase(z)*i
     *
     * @param A Complex of radians
     *
     * @return a Complex with elements ln(A)
     *************************************************************************/
    public static Complex log ( Complex A ) {
        double[] realA = A.getDataReal();
        double[] imagA = A.getDataImag();
        double a,b,r;

        Complex result = new Complex( A.shape() );
        double[] resultReal = result.getDataReal();
        double[] resultImag = result.getDataImag();

        for ( int i = 0; i < resultReal.length; i++ ) {
            a = realA[i];
            b = imagA[i];
            r = Math.sqrt( a*a + b*b );
            resultReal[i] = Math.log(r);
            resultImag[i] = Math.atan2(b,a);
        }
        return result;
    }

    /**************************************************************************
     * <p>Compute the natural log of A element-wise.
     *
     * Note for complex z = a+bi, ln(z) = ln(mag(z)) + phase(z)*i
     *
     * @param A Complex of radians
     *
     * @return a Complex with elements ln(A)
     *************************************************************************/
    public static Numeric log ( Numeric A ) {
        Numeric result = A.copy();
        double[] resultReal = result.getData();

        for ( int i = 0; i < resultReal.length; i++ )
            resultReal[i] = Math.log( resultReal[i] );

        return result;
    }

    /**************************************************************************
     * <p>Compute the natural log of A element-wise.
     *
     * @param A Complex of radians
     *
     * @return a Complex with elements ln(A)
     *************************************************************************/
    public static Complex sum ( Complex A ) {
        double[] realA = A.getDataReal();
        double[] imagA = A.getDataImag();
        double sumReal = 0, sumImag = 0;

        for ( int i = 0; i < realA.length; i++ ) {
            sumReal += realA[i];
            sumImag += imagA[i];
        }

        return Complex.Scalar(sumReal, sumImag);
    }

    /**************************************************************************
     * <p>Compute the natural log of A element-wise.
     *
     * @param A Complex of radians
     *
     * @return a Complex with elements ln(A)
     *************************************************************************/
    public static Numeric sum ( Numeric A ) {
        double[] data = A.getData();
        double sumReal = 0, sumImag = 0;

        for (double d : data) sumReal += d;

        return Numeric.Scalar( sumReal );
    }

    /**************************************************************************
     * <p>Compute the complex conjugate of A element-wise.
     *
     * @param A a Complex
     *
     * @return a Complex with elements complex conjugate(A)
     *************************************************************************/
    public static Complex conjugate ( Complex A ) {

        Complex result = new Complex( A.shape() );
        double[] resultData = result.getDataImag();
        double[] imagA = result.getDataImag();

        for ( int i = 0; i < imagA.length; i++ ) {
            resultData[i] = -imagA[i];
        }

        return result;
    }

    /**************************************************************************
     * <p>Compute the sine of A element-wise. A is in radians.
     *
     * @param A Complex of radians
     *
     * @return a Complex with elements sin(A)
     *************************************************************************/
    public static Complex sin ( Complex A ) {
        double[] realA = A.getDataReal();
        double[] imagA = A.getDataImag();
        double a,b;

        Complex result = new Complex( A.shape() );
        double[] resultReal = result.getDataReal();
        double[] resultImag = result.getDataImag();

        for ( int i = 0; i < resultReal.length; i++ ) {
            a = realA[i];
            b = imagA[i];
            resultReal[i] = Math.sin(a) * Math.cosh(b);
            resultImag[i] = Math.cos(a) * Math.sinh(b);
        }
        return result;
    }

    /**************************************************************************
     * <p>Compute the sine of A element-wise. A is in radians.
     *
     * @param A Complex of radians
     *
     * @return a Complex with elements sin(A)
     *************************************************************************/
    public static Numeric sin ( Numeric A ) {
        double[] data = A.getData();

        Numeric result = new Numeric( A.shape() );
        double[] resultData = result.getData();

        for ( int i = 0; i < resultData.length; i++ ) {
            resultData[i] = Math.sin( data[i] );
        }
        return result;
    }

    /**************************************************************************
     * <p>Compute the sine of A element-wise. A is in degrees. Note that the
     * conversion from radians to degrees is not exact.
     *
     * @param A Complex of radians
     *
     * @return a Complex with elements sin(A)
     *************************************************************************/
    public static Complex sind ( Complex A ) {
        double[] realA = A.getDataReal();
        double[] imagA = A.getDataImag();
        double a,b, toRad = 0.017453292519943295;

        Complex result = new Complex( A.shape() );
        double[] resultReal = result.getDataReal();
        double[] resultImag = result.getDataImag();

        for ( int i = 0; i < resultReal.length; i++ ) {
            a = realA[i] * toRad;
            b = imagA[i] * toRad;
            resultReal[i] = Math.sin(a) * Math.cosh(b);
            resultImag[i] = Math.cos(a) * Math.sinh(b);
        }
        return result;
    }

    /**************************************************************************
     * <p>Compute the sine of A element-wise. A is in degrees. Note that the
     * conversion from radians to degrees is not exact.
     *
     * @param A Complex of radians
     *
     * @return a Complex with elements sin(A)
     *************************************************************************/
    public static Numeric sind ( Numeric A ) {
        double[] realA = A.getData();
        double toRad = 0.017453292519943295;

        Numeric result = new Numeric( A.shape() );
        double[] resultData = result.getData();

        for ( int i = 0; i < resultData.length; i++ )
            resultData[i] = Math.sin( realA[i] * toRad );

        return result;
    }

    /**************************************************************************
     * <p>Compute the cosine of A element-wise. A is in radians.
     *
     * @param A Complex of radians
     *
     * @return a Complex with elements cos(A)
     *************************************************************************/
    public static Complex cos ( Complex A ) {
        double[] realA = A.getDataReal();
        double[] imagA = A.getDataImag();
        double a,b;

        Complex result = new Complex( A.shape() );
        double[] resultReal = result.getDataReal();
        double[] resultImag = result.getDataImag();

        for ( int i = 0; i < resultReal.length; i++ ) {
            a = realA[i];
            b = imagA[i];
            resultReal[i] = Math.cos(a) * Math.cosh(b);
            resultImag[i] = -1 * Math.sin(a) * Math.sinh(b);
        }
        return result;
    }

    /**************************************************************************
     * <p>Compute the cosine of A element-wise. A is in radians.
     *
     * @param A Complex of radians
     *
     * @return a Complex with elements cos(A)
     *************************************************************************/
    public static Numeric cos ( Numeric A ) {
        double[] realA = A.getData();

        Numeric result = new Numeric( A.shape() );
        double[] resultReal = result.getData();

        for ( int i = 0; i < resultReal.length; i++ )
            resultReal[i] = Math.cos( realA[i] );
        return result;
    }

    /**************************************************************************
     * <p>Compute the cosine of A element-wise. A is in degrees. Note that the
     * conversion from radians to degrees is not exact
     *
     * @param A Complex of radians
     *
     * @return a Complex with elements cos(A)
     *************************************************************************/
    public static Complex cosd ( Complex A ) {
        double[] realA = A.getDataReal();
        double[] imagA = A.getDataImag();
        double a,b, toRad = 0.017453292519943295;

        Complex result = new Complex( A.shape() );
        double[] resultReal = result.getDataReal();
        double[] resultImag = result.getDataImag();

        for ( int i = 0; i < resultReal.length; i++ ) {
            a = realA[i] * toRad;
            b = imagA[i] * toRad;
            resultReal[i] = Math.cos(a) * Math.cosh(b);
            resultImag[i] = -1 * Math.sin(a) * Math.sinh(b);
        }
        return result;
    }

    /**************************************************************************
     * <p>Compute the cosine of A element-wise. A is in degrees. Note that the
     * conversion from radians to degrees is not exact
     *
     * @param A Complex of radians
     *
     * @return a Complex with elements cos(A)
     *************************************************************************/
    public static Numeric cosd ( Numeric A ) {
        double[] realA = A.getData();
        double toRad = 0.017453292519943295;

        Numeric result = new Numeric( A.shape() );
        double[] resultReal = result.getData();

        for ( int i = 0; i < resultReal.length; i++ ) {
            resultReal[i] = Math.cos( realA[i] * toRad );
        }
        return result;
    }

    /**************************************************************************
     * <p>Compute the tangent of A element-wise. A is in radians.
     *
     * @param A Complex of radians
     *
     * @return a Complex with elements tan(A)
     *************************************************************************/
    public static Complex tan ( Complex A ) {
        double[] realA = A.getDataReal();
        double[] imagA = A.getDataImag();
        double a,b,cos2acosh2b;

        Complex result = new Complex( A.shape() );
        double[] resultReal = result.getDataReal();
        double[] resultImag = result.getDataImag();

        for ( int i = 0; i < resultReal.length; i++ ) {
            a = realA[i] * 2;
            b = imagA[i] * 2;
            cos2acosh2b = Math.cos(a) + Math.cosh(b);
            resultReal[i] = Math.sin(a)  / cos2acosh2b;
            resultImag[i] = Math.sinh(b) / cos2acosh2b;
        }
        return result;
    }

    /**************************************************************************
     * <p>Compute the tangent of A element-wise. A is in radians.
     *
     * @param A Complex of radians
     *
     * @return a Complex with elements tan(A)
     *************************************************************************/
    public static Numeric tan ( Numeric A ) {
        double[] realA = A.getData();

        Numeric result = new Numeric( A.shape() );
        double[] resultReal = result.getData();

        for ( int i = 0; i < resultReal.length; i++ )
            resultReal[i] = Math.tan( realA[i] );

        return result;
    }

    /**************************************************************************
     * <p>Compute the tangent of A element-wise. A is in degrees, Not that the
     * conversion to radians is not exact.
     *
     * @param A Complex of radians
     *
     * @return a Complex with elements tan(A)
     *************************************************************************/
    public static Complex tand ( Complex A ) {
        double[] realA = A.getDataReal();
        double[] imagA = A.getDataImag();
        double a,b,cos2acosh2b,toRad = 0.017453292519943295;

        Complex result = new Complex( A.shape() );
        double[] resultReal = result.getDataReal();
        double[] resultImag = result.getDataImag();

        for ( int i = 0; i < resultReal.length; i++ ) {
            a = realA[i] * 2 * toRad;
            b = imagA[i] * 2 * toRad;
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
     * @param A Complex of radians
     *
     * @return a Complex with elements tan(A)
     *************************************************************************/
    public static Numeric tand ( Numeric A ) {
        double[] realA = A.getData();
        double toRad = 0.017453292519943295;

        Numeric result = new Numeric( A.shape() );
        double[] resultReal = result.getData();
        double[] resultImag = result.getData();

        for ( int i = 0; i < resultReal.length; i++ )
            resultReal[i] = Math.tan( realA[i] * toRad );

        return result;
    }

    /**************************************************************************
     * <p>Compute the hyperbolic cosine of A element-wise. A is in radians.
     *
     * @param A Complex of radians
     *
     * @return a Complex with elements cosh(A)
     *************************************************************************/
    public static Complex cosh ( Complex A ) {
        double[] realA = A.getDataReal();
        double[] imagA = A.getDataImag();
        double a,b;

        Complex result = new Complex( A.shape() );
        double[] resultReal = result.getDataReal();
        double[] resultImag = result.getDataImag();

        for ( int i = 0; i < resultReal.length; i++ ) {
            a = realA[i];
            b = imagA[i];
            resultReal[i] = Math.cosh(a) * Math.cos(b);
            resultImag[i] = Math.sinh(a) * Math.sin(b);
        }
        return result;
    }

    /**************************************************************************
     * <p>Compute the hyperbolic cosine of A element-wise. A is in radians.
     *
     * @param A Complex of radians
     *
     * @return a Complex with elements cosh(A)
     *************************************************************************/
    public static Numeric cosh ( Numeric A ) {
        double[] realA = A.getData();

        Numeric result = new Numeric( A.shape() );
        double[] resultReal = result.getData();

        for ( int i = 0; i < resultReal.length; i++ )
            resultReal[i] = Math.cosh( realA[i] );

        return result;
    }

    /**************************************************************************
     * <p>Compute the hyperbolic sine of A element-wise. A is in radians.
     *
     * @param A Complex of radians
     *
     * @return a Complex with elements sinh(A)
     *************************************************************************/
    public static Complex sinh ( Complex A ) {
        double[] realA = A.getDataReal();
        double[] imagA = A.getDataImag();
        double a,b;

        Complex result = new Complex( A.shape() );
        double[] resultReal = result.getDataReal();
        double[] resultImag = result.getDataImag();

        for ( int i = 0; i < resultReal.length; i++ ) {
            a = realA[i];
            b = imagA[i];
            resultReal[i] = Math.sinh(a) * Math.cos(b);
            resultImag[i] = Math.cosh(a) * Math.sin(b);
        }
        return result;
    }

    /**************************************************************************
     * <p>Compute the hyperbolic sine of A element-wise. A is in radians.
     *
     * @param A Complex of radians
     *
     * @return a Complex with elements sinh(A)
     *************************************************************************/
    public static Numeric sinh ( Numeric A ) {
        double[] realA = A.getData();

        Numeric result = new Numeric( A.shape() );
        double[] resultReal = result.getData();

        for ( int i = 0; i < resultReal.length; i++ ) {
            resultReal[i] = Math.sinh( realA[i] );
        }
        return result;
    }

    /**************************************************************************
     * <p>Compute the hyperbolic tangent of A element-wise. A is in radians.
     *
     * @param A Complex of radians
     *
     * @return a Complex with elements tanh(A)
     *************************************************************************/
    public static Complex tanh ( Complex A ) {
        double[] realA = A.getDataReal();
        double[] imagA = A.getDataImag();
        double a,b,cosh2acos2b;

        Complex result = new Complex( A.shape() );
        double[] resultReal = result.getDataReal();
        double[] resultImag = result.getDataImag();

        for ( int i = 0; i < resultReal.length; i++ ) {
            a = realA[i] * 2;
            b = imagA[i] * 2;
            cosh2acos2b = Math.cosh(a) + Math.cos(b);
            resultReal[i] = Math.sinh(a) / cosh2acos2b;
            resultImag[i] = Math.sin(b)  / cosh2acos2b;
        }
        return result;
    }

    /**************************************************************************
     * <p>Compute the hyperbolic tangent of A element-wise. A is in radians.
     *
     * @param A Complex of radians
     *
     * @return a Complex with elements tanh(A)
     *************************************************************************/
    public static Numeric tanh ( Numeric A ) {
        double[] realA = A.getData();

        Numeric result = new Numeric( A.shape() );
        double[] resultReal = result.getData();

        for ( int i = 0; i < resultReal.length; i++ )
            resultReal[i] = Math.tanh( realA[i] );

        return result;
    }

    /**************************************************************************
     * <p>Compute the inverse sine of A element-wise.
     *
     * @param A Numeric
     *
     * @return a Numeric with elements arcsin(A)
     *************************************************************************/
    public static Numeric asin ( Numeric A ) {
        double[] data = A.getData();

        Numeric result = new Numeric( A.shape() );
        double[] resultData = result.getData();

        for ( int i = 0; i < resultData.length; i++ ) {
            resultData[i] = Math.asin( data[i] );
        }
        return result;
    }

    /**************************************************************************
     * <p>Compute the inverse sine of A element-wise.
     *
     * @param A Numeric
     *
     * @return a Numeric with elements arcsin(A)
     *************************************************************************/
    public static Numeric acos ( Numeric A ) {
        double[] data = A.getData();

        Numeric result = new Numeric( A.shape() );
        double[] resultData = result.getData();

        for ( int i = 0; i < resultData.length; i++ ) {
            resultData[i] = Math.acos( data[i] );
        }
        return result;
    }

    /**************************************************************************
     * <p>Compute the inverse sine of A element-wise.
     *
     * @param A Numeric
     *
     * @return a Numeric with elements arcsin(A)
     *************************************************************************/
    public static Numeric atan ( Numeric A ) {
        double[] data = A.getData();

        Numeric result = new Numeric( A.shape() );
        double[] resultData = result.getData();

        for ( int i = 0; i < resultData.length; i++ ) {
            resultData[i] = Math.atan( data[i] );
        }
        return result;
    }

    /**************************************************************************
     * <p>Compare A and B element wise. If one of the values is scalar, compare
     * that value to the other Complex. Only compares real components as
     * comparison is not well defined for complex values.
     *
     * @param A the first Complex
     * @param B the second Complex
     *
     * @return a Logical that indexes A > B
     *************************************************************************/
    public static Logical gre ( Complex A, Complex B ) {
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
     * that value to the other Complex. Only compares real components as
     * comparison is not well defined for complex values.
     *
     * @param A the first Complex
     * @param B the second Complex
     *
     * @return a Logical that indexes A >= B
     *************************************************************************/
    public static Logical greq ( Complex A, Complex B ) {
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
     * that value to the other Complex. Only compares real components as
     * comparison is not well defined for complex values.
     *
     * @param A the first Complex
     * @param B the second Complex
     *
     * @return a Logical that indexes A < B
     *************************************************************************/
    public static Logical less ( Complex A, Complex B ) {
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
     * that value to the other Complex. Only compares real components as
     * comparison is not well defined for complex values.
     *
     * @param A the first Complex
     * @param B the second Complex
     *
     * @return a Logical that indexes A <= B
     *************************************************************************/
    public static Logical leq ( Complex A, Complex B ) {
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
     * @param A the first Complex
     * @param B the second Complex
     *
     * @return a Logical that indexes A == B
     *************************************************************************/
    public static Logical equal ( Complex A, Complex B, double tolerance ) {
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
     * @param A the first Complex
     * @param B the second Complex
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
     * @param A the first Complex
     * @param B the second Complex
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


    /*Ensure that either one Complex is scalar or both have the same shape*/
    private static void validateDimensionsFatal ( NDArray N1, NDArray N2 ) {
        if ( NDArray.dimensionsMatch(N1,N2) )
            return;
        if ( N1.isScalar() || N2.isScalar() )
            return;
        if ( N1.isVector() && N2.isVector() && (N1.getSize() == N2.getSize()) )
            return;

        throw new IllegalDimensionException(
            "Element wise operation: dimensions must be equal or one NDArray must be scalar"
        );
    }




}
