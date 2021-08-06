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

import com.evanstella.jnp.core.Complex;
import com.evanstella.jnp.core.IllegalDimensionException;
import com.evanstella.jnp.core.Numeric;

import java.util.concurrent.CountDownLatch;

/******************************************************************************
 * <p>MatrixExecutor a is a handler for multithreaded matrix operations. Once
 * an instance of the class has been created, it can be used to perform most
 * of the same operations in Matrix using multithreading.
 *
 * @author Evan Stella
 *****************************************************************************/
public final class MatrixExecutor extends ParallelExecutor {

    /**************************************************************************
     * <p>Constructor.
     *
     * @param threadCount   The number of threads to create for the
     *                      handler.
     *************************************************************************/
    public MatrixExecutor ( int threadCount ) {
        super(threadCount);
    }

    /**************************************************************************
     * <p>Transpose the inputted matrix
     *
     * @param matrix    The matrix to transpose
     *
     * @return a matrix equal to the transpose of the input
     *************************************************************************/
    public Numeric transpose ( Numeric matrix ) {
        Matrix.validateMatrixFatal( matrix );

        int[] shape = matrix.shape();
        double[] data = matrix.getData();
        final int c = shape[1], r = shape[0];

        Numeric transposed = new Numeric( shape[1], shape[0] );
        double[] newData = transposed.getData();

        //create worker
        CountDownLatch count = new CountDownLatch(threadCount);
        class worker implements Runnable {
            final int startIdx, endIdx;
            public worker (int start, int end) {startIdx = start; endIdx = end;}
            public void run ( ) {
                try {
                    int curR, curC;
                    for (int i = startIdx; i < endIdx; i++) {
                        curR = i / c;
                        curC = i % c;
                        newData[ curC*r+curR ] = data[i];
                    }
                }
                finally { count.countDown(); }
            }
        }

        int start = 0, increment = data.length / threadCount;
        for ( int i = 0; i < threadCount-1; i++ ) {
            executorService.execute( new worker( start, start+increment ) );
            start = start+increment;
        }
        executorService.execute( new worker( start, data.length ) );
        this.await( count );
        return transposed;
    }

    /**************************************************************************
     * <p>Transpose the inputted matrix
     *
     * @param matrix    The matrix to transpose
     *
     * @return a matrix equal to the transpose of the input
     **************************************************************************/
    public Complex transpose ( Complex matrix ) {
        Matrix.validateMatrixFatal( matrix );

        int[] shape = matrix.shape();
        double[] dataR = matrix.getDataReal();
        double[] dataI = matrix.getDataImag();
        final int c = shape[1], r = shape[0];

        Complex transposed = new Complex( shape[1], shape[0] );
        double[] newDataR = transposed.getDataReal();
        double[] newDataI = transposed.getDataImag();

        //create worker
        CountDownLatch count = new CountDownLatch( threadCount );
        class worker implements Runnable {
            final int startIdx, endIdx;
            public worker (int start, int end) {startIdx = start; endIdx = end;}
            public void run ( ) {
                try {
                    int curR, curC;
                    for (int i = startIdx; i < endIdx; i++) {
                        curR = i / c;
                        curC = i % c;
                        newDataR[ curC*r+curR ] = dataR[i];
                        newDataI[ curC*r+curR ] = dataI[i];
                    }
                }
                finally { count.countDown(); }
            }
        }

        int start = 0, increment = dataR.length / threadCount;
        for ( int i = 0; i < threadCount-1; i++ ) {
            executorService.execute( new worker( start, start+increment ) );
            start = start+increment;
        }
        executorService.execute( new worker( start, dataR.length ) );
        this.await( count );
        return transposed;
    }

    /**************************************************************************
     * <p>Conjugate transpose the inputted matrix. Not that the complex
     * transpose equal to the transpose of a matrix except the conjugate of
     * each element is taken
     *
     * @param matrix    The matrix to transpose
     *
     * @return a matrix equal to the transpose of the input
     *************************************************************************/
    public Complex ctranspose ( Complex matrix ) {
        Matrix.validateMatrixFatal( matrix );

        int[] shape = matrix.shape();
        double[] dataR = matrix.getDataReal();
        double[] dataI = matrix.getDataImag();
        final int c = shape[1], r = shape[0];

        Complex transposed = new Complex( shape[1], shape[0] );
        double[] newDataR = transposed.getDataReal();
        double[] newDataI = transposed.getDataImag();

        //create worker
        CountDownLatch count = new CountDownLatch( threadCount );
        class worker implements Runnable {
            final int startIdx, endIdx;
            public worker (int start, int end) {startIdx = start; endIdx = end;}
            public void run ( ) {
                try {
                    int curR, curC;
                    for (int i = startIdx; i < endIdx; i++) {
                        curR = i / c;
                        curC = i % c;
                        newDataR[ curC*r+curR ] =  dataR[i];
                        newDataI[ curC*r+curR ] = -dataI[i];
                    }
                }
                finally { count.countDown(); }
            }
        }

        int start = 0, increment = dataR.length / threadCount;
        for ( int i = 0; i < threadCount-1; i++ ) {
            executorService.execute( new worker( start, start+increment ) );
            start = start+increment;
        }
        executorService.execute( new worker( start, dataR.length ) );
        this.await( count );
        return transposed;
    }


    /**************************************************************************
     * <p>Matrix multiplication of m1 and m2. Number of columns of m1 must be
     * equal to the number of rows of m2.
     *
     * @param m1    The first matrix
     * @param m2    The second matrix
     *
     * @return a matrix m1*m2.
     **************************************************************************/
    public Numeric mul ( Numeric m1, Numeric m2 ) {
        Matrix.validateMatrixFatal( m1 );
        Matrix.validateMatrixFatal( m2 );
        if ( m1.shape()[1] != m2.shape()[0] )
            throw new IllegalDimensionException(
                "Matrix multiplication: invalid dimensions."
            );

        final int c1 = m1.shape()[1], c2 = m2.shape()[1];
        final int r = m1.shape()[0], c = m2.shape()[1];
        double[] X = m1.getData(); double[] Y = m2.getData();
        Numeric result = new Numeric( r, c );
        double[] Z = result.getData();

        //create worker
        CountDownLatch count = new CountDownLatch( threadCount );
        class worker implements Runnable {
            final int i,j;
            public worker ( int i, int j ) { this.i = i; this.j = j; }
            public void run ( ) {
                try {
                    int ind = i*c+j;
                    Z[ind] = 0;
                    for( int k = 0; k < c1; k++ )
                        Z[ind] += X[i*c1+k] * Y[k*c2+j];
                }
                finally { count.countDown(); }
            }
        }

        for ( int i = 0; i < r; i++ )
            for( int j = 0; j < c; j++ )
                executorService.execute( new worker(i, j) );

        this.await( count );
        return result;
    }




}
