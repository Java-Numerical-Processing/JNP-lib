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

package com.evanstella.jnp.analysis;

import com.evanstella.jnp.core.*;
import com.evanstella.jnp.math.ExecutionInterruptedException;

import java.util.concurrent.CountDownLatch;

/******************************************************************************
 * <p>IntegralExecutor is a handler for multithreaded integration computations.
 *
 * @author Evan Stella
 *****************************************************************************/
public class IntegralExecutor extends ParallelExecutor {

    /**************************************************************************
     * <p>Constructor.
     *
     * @param numThreads    The number of threads to create for the
     *                      handler.
     *************************************************************************/
    public IntegralExecutor ( int numThreads ) { super( numThreads ); }

    /**************************************************************************
     * <p>Constructor. Create IntegralExecutor from ParallelExecutor so the
     * two are bound: shutting down one will shut down the other.
     *
     * @param P             The ParallelExecutor to bind to this new Executor
     *************************************************************************/
    public IntegralExecutor ( ParallelExecutor P ) { P.bind( this ); }

    /**************************************************************************
     * <p>Computes the integral of inputted the numeric data using trapezoidal
     * numerical integration. Y is the output of some numeric function on X.
     * The data in X and Y need not be evenly spaced but they must be vectors
     * of the same size.
     *
     * @param X        The x (independent) data
     * @param Y        The y (dependent) data
     *
     * @return A scalar with the computed value.
     *************************************************************************/
    public Numeric integral ( Numeric X, Numeric Y ) {
        validateIntegral1DInputs( X, Y );
        double[] x = X.getData();
        double[] y = Y.getData();
        double[] sumParts = new double[threadCount];

        //create worker
        CountDownLatch count = new CountDownLatch( threadCount );
        class worker implements Runnable {
            final int ind, startIdx, endIdx;
            public worker ( int thread, int start, int end ) {
                ind = thread; startIdx = start; endIdx = end;
            }
            public void run ( ) {
                try {
                    sumParts[ind] = 0.0;
                    for ( int i = startIdx+1; i < endIdx; i++ )
                        sumParts[ind] += ( y[i-1] + y[i] ) * ( x[i] - x[i-1] );
                }
                finally { count.countDown(); }
            }
        }

        int start = 0, increment = x.length / threadCount;
        for ( int i = 0; i < threadCount-1; i++ ) {
            executorService.execute( new worker( i, start, start+increment ) );
            start += increment;
        }
        executorService.execute( new worker( threadCount-1, start, x.length ) );
        this.await( count );

        double sum = 0.0;
        for ( double d : sumParts )
            sum += d;
        return Numeric.Scalar( sum/2 );
    }

    /**************************************************************************
     * <p>Computes the integral of the inputted Numeric expression over the
     * inputted bounds using Simpsons 1/3 rule for numerical integration using
     * quadratics.
     *
     * @param expr          The numeric expression to evaluate
     * @param start         The starting bound (inclusive)
     * @param end           The ending bound (inclusive)
     *
     * @return A scalar with the computed value.
     *************************************************************************/
    public Numeric integral ( NumericExpression1D expr, double start, double end ) {
        if ( start >= end )
            throw new IllegalArgumentException (
                "Start value must be less than ending value.");
        double[] sumParts = new double[threadCount];
        int numPts = 10000;
        // number of points per thread, must be even for Simpsons rule
        numPts /= threadCount;
        numPts = ( numPts % 2 == 0 ) ? numPts : numPts + 1;
        // 2 extra for first and last evaluation
        double h = ( end - start ) / ( 2 + numPts * threadCount - 1 );

        //create worker
        CountDownLatch count = new CountDownLatch( threadCount );
        final int increment = numPts;
        class worker implements Runnable {
            final int ind;
            final double startVal;
            public worker ( int thread, double start ) { ind = thread; startVal = start; }
            public void run ( ) {
                try {
                    sumParts[ind] = 0.0;
                    double coef = 2.0, x = startVal;
                    for ( int i = 0; i < increment; i++ ) {
                        coef = ( coef == 2.0 ) ? 4.0 : 2.0;
                        sumParts[ind] += expr.evaluate( x ) * coef;
                        x += h;
                    }
                } catch ( RuntimeException ex ) {
                    throw new ExecutionInterruptedException(
                        "Error executing numeric expression:\n" + ex );
                } finally { count.countDown(); }
            }
        }

        start = start + h;
        for ( int i = 0; i < threadCount; i++ ) {
            executorService.execute( new worker( i, start ) );
            start += increment * h;
        }
        this.await( count );

        // the first and last value with a coefficient of 1
        double sum = expr.evaluate( start ) + expr.evaluate( end );
        for ( double d : sumParts ) sum += d;

        return Numeric.Scalar( sum * h / 3 );
    }


                        /* Internal Methods */
    private void validateIntegral1DInputs ( NDArray A, NDArray B ) {
        if ( A.isVector() && B.isVector() )
            return;
        if ( A.getSize() == B.getSize() )
            return;

        throw new IllegalDimensionException (
            "Integral inputs must be vectors of the same size."
        );

    }
}
