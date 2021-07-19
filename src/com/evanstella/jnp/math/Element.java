package com.evanstella.jnp.math;

import com.evanstella.jnp.core.Numeric;

public final class Element {

    // no instances for you
    private Element () {}


    /**************************************************************************
     * Calculates each element in N to the power of the inputted exponent. If
     * the element is both positive and real, just use built in Math.pow;
     * otherwise use DeMoivre's Theorem: z^n = r^n(cos(n*theta)+i*sin(n*theta))
     *
     * @param N         the Numeric to calculate pow for
     * @param exponent  the exponent to raise each element to
     *
     * @return a Numeric whose elements are equal to N^exponent
     *************************************************************************/
    public static Numeric pow ( Numeric N, double exponent ) {
        Numeric result = N.copy();
        double[] realData = result.getDataReal();
        double[] imagData = result.getDataImag();
        boolean isComplex = imagData != null;

        double x,y,r,theta;
        for ( int i = 0; i < realData.length; i++ ) {
            x = realData[i];
            if ( !isComplex ) y = 0;
            else y = imagData[i];
            // only do all the extra work if we have to
            if (x >= 0 && y == 0) {
                realData[i] = Math.pow(x, exponent);
                continue;
            }
            if ( !isComplex ) {
                isComplex = true;
                imagData = result.initializeDataImag();
            }
            r = Math.sqrt( x*x + y*y );
            r = Math.pow( r, exponent );
            theta = Math.atan2(y,x);
            theta = exponent * theta;
            realData[i] = r * Math.cos( theta );
            imagData[i] = r * Math.sin( theta );
        }

        return result;
    }

    /**************************************************************************
     * Calculates N^exponent for each element in N. The following is used:
     *
     * <p>For w = a+bi, z = c+di: w^z = e^(z*ln(r)+i*theta) for r = sqrt(a^2+b^2)
     * and theta = atan2(b,a) ...
     * <p>w^z = r^c * e^(-d*theta) *
     *          [ cos(d*ln(r) + c*theta) + i*sin(d*ln(r) + c*theta) ]
     *
     * @param expReal   the real component of the exponent
     * @param expImag   the imaginary component of the element
     * @param N         the Numeric exponent
     *
     * @return a Numeric with the same size of N whose elements = base^N
     *************************************************************************/
    public static Numeric pow ( Numeric N, double expReal, double expImag ) {
        Numeric result = N.copy();
        double[] realData = result.getDataReal();
        double[] imagData = result.getDataImag();

        // we are now assuming complex output
        if ( imagData == null )
            imagData = result.initializeDataImag();

        double a,b,r,rc,theta,Omega;
        for ( int i = 0; i < realData.length; i++ ) {
            a = realData[i];
            b = imagData[i];
            r = Math.sqrt( a*a + b*b );
            theta = Math.atan2( b, a );
            rc = Math.pow( r, expReal ) * Math.exp( -expImag * theta );
            Omega = ( expImag * Math.log(r) ) + ( expReal * theta );
            realData[i] = rc * Math.cos(Omega);
            imagData[i] = rc * Math.sin(Omega);
        }

        return result;
    }

    /**************************************************************************
     * Calculates base^N for each element in N. If the base is positive and
     * real and the exponent is real, the function just computes base^N to
     * avoid extra computation. Otherwise the following is used:
     *
     * <p>For w = a+bi, z = c+di: w^z = e^(z*ln(r)+i*theta) for r = sqrt(a^2+b^2)
     * and theta = atan2(b,a) ...
     * <p>w^z = r^c * e^(-d*theta) *
     *          [ cos(d*ln(r) + c*theta) + i*sin(d*ln(r) + c*theta) ]
     *
     * @param base      the base be raised to each element
     * @param N         the Numeric exponent
     *
     * @return a Numeric with the same size of N whose elements = base^N
     *************************************************************************/
    public static Numeric pow ( double base, Numeric N ) {
        Numeric result = N.copy();
        double[] realData = result.getDataReal();
        double[] imagData = result.getDataImag();
        boolean isComplex = imagData != null;
        boolean negBase = base < 0.0;

        double c,d,r,rc,theta,Omega;
        // base has no complex component, r = sqrt(a^2) = abs(a) and theta = 0 or pi;
        r = (base <= 0.0) ? 0.0 - base : base;
        theta = (base <= 0.0) ? Math.PI : 0;
        for ( int i = 0; i < realData.length; i++ ) {
            c = realData[i];
            if ( !isComplex ) d = 0;
            else d = imagData[i];
            if ( !negBase && d == 0 ) {
                realData[i] = Math.pow( base, c );
                continue;
            }
            if ( !isComplex ) {
                isComplex = true;
                imagData = result.initializeDataImag();
            }
            rc = Math.pow( r, c ) * Math.exp( -d * theta );
            Omega = d * Math.log(r) + c * theta;
            realData[i] = rc * Math.cos(Omega);
            imagData[i] = rc * Math.sin(Omega);
        }

        return result;
    }

    /**************************************************************************
     * Calculates base^N for each element in N. The following is used:
     *
     * <p>For w = a+bi, z = c+di: w^z = e^(z*ln(r)+i*theta) for r = sqrt(a^2+b^2)
     * and theta = atan2(b,a) ...
     * <p>w^z = r^c * e^(-d*theta) *
     *          [ cos(d*ln(r) + c*theta) + i*sin(d*ln(r) + c*theta) ]
     *
     * @param baseReal  the real component of the base be raised to each element
     * @param baseImag  the imaginary component of the base
     * @param N         the Numeric exponent
     *
     * @return a Numeric with the same size of N whose elements = base^N
     *************************************************************************/
    public static Numeric pow ( double baseReal, double baseImag, Numeric N ) {
        Numeric result = N.copy();
        double[] realData = result.getDataReal();
        double[] imagData = result.getDataImag();

        // we are now assuming complex output
        if ( imagData == null )
            imagData = result.initializeDataImag();

        double c,d,r,rc,theta,Omega;
        r = Math.sqrt( baseReal*baseReal + baseImag*baseImag );
        theta = Math.atan2(baseImag, baseReal);

        for ( int i = 0; i < realData.length; i++ ) {
            c = realData[i];
            d = imagData[i];
            rc = Math.pow( r, c ) * Math.exp( -d * theta );
            Omega = ( d * Math.log(r) ) + ( c * theta );
            realData[i] = rc * Math.cos(Omega);
            imagData[i] = rc * Math.sin(Omega);
        }

        return result;
    }

}
