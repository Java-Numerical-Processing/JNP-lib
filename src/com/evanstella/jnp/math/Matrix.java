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
import com.evanstella.jnp.core.NDArray;
import com.evanstella.jnp.core.Numeric;

/******************************************************************************
 * <p>Matrix encapsulates all of the Matrix operations that can be done on
 * NDArrays. These operations are implemented as static methods and are
 * executed in series. For multithreaded operations use MatrixExecutor.
 *
 * @author Evan Stella
 *****************************************************************************/
public final class Matrix {

    // no instances for you
    private Matrix ( ) {}

    /**************************************************************************
     * <p>Matrix multiplication of m1 and m2. Number of columns of m1 must be
     * equal to the number of rows of m2.
     *
     * @param m1    The first matrix
     * @param m2    The second matrix
     *
     * @return a matrix m1*m2.
     *************************************************************************/
    public static Numeric mul ( Numeric m1, Numeric m2 ) {
        Matrix.validateMatrixFatal( m1 );
        Matrix.validateMatrixFatal( m2 );
        int c1 = m1.shape()[1], c2 = m2.shape()[1];

        if ( m1.shape()[1] != m2.shape()[0] )
            throw new IllegalDimensionException(
                "Matrix multiplication: invalid dimensions."
            );

        int r = m1.shape()[0], c = m2.shape()[1];
        double[] X = m1.getData();
        double[] Y = m2.getData();
        Numeric result = new Numeric( r, c );
        double[] Z = result.getData();

        for ( int ind, i = 0; i < r ;i++ ) {
            for( int j = 0; j < c; j++ ) {
                ind = i*c+j;
                Z[ind] = 0;
                for( int k = 0; k < c1; k++ ) {
                    Z[ind] += X[i*c1+k] * Y[k*c2+j];
                }
            }
        }

        return result;
    }

    /**************************************************************************
     * <p>Compute the trace of the input matrix
     *
     * @param matrix    the matrix
     *
     * @return the trace of the matrix.
     **************************************************************************/
    public static Numeric trace ( Numeric matrix ) {
        Matrix.validateSquareMatrix( matrix );
        double[] data = matrix.getData();

        double tr = 0;
        int cols = matrix.shape()[0];
        for ( int i = 0; i < cols; i++ ) {
            tr += data[ i * cols + i ];
        }

        return Numeric.Scalar( tr );
    }

    /**************************************************************************
     * <p>Compute the trace of the input matrix
     *
     * @param matrix    the matrix
     *
     * @return the trace of the matrix.
     **************************************************************************/
    public static Complex trace ( Complex matrix ) {
        Matrix.validateSquareMatrix( matrix );
        double[] dataR = matrix.getDataReal();
        double[] dataI = matrix.getDataImag();

        double trR = 0;
        double trI = 0;
        int cols = matrix.shape()[0];
        for ( int ind, i = 0; i < cols; i++ ) {
            ind = i*cols+i;
            trR += dataR[ind];
            trI += dataI[ind];
        }

        return Complex.Scalar( trR, trI );
    }


    /*package private*/ static void validateMatrixFatal ( NDArray N ) {
        if ( N.shape().length == 2 )
            return;

        throw new IllegalDimensionException(
            "Matrix operation: Inputs must be matrices (2 dimensions)."
        );
    }

    /*package private*/ static void validateSquareMatrix( NDArray N ) {
        validateMatrixFatal( N );
        int[] shape = N.shape();
        if ( shape[0] == shape[1] )
            return;

        throw new IllegalDimensionException(
            "Input matrix must be square (NxN)."
        );
    }
}
