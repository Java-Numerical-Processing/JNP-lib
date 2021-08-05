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

import com.evanstella.jnp.core.IllegalDimensionException;
import com.evanstella.jnp.core.NDArray;
import com.evanstella.jnp.core.Numeric;

/******************************************************************************
 * Matrix encapsulates all of the Matrix operations that can be done on
 * NDArrays. TODO
 *
 * @author Evan Stella
 *****************************************************************************/
public final class Matrix {

    // no instances for you
    private Matrix ( ) {}

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
        Numeric result = new Numeric( r, c);
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

    public static void validateMatrixFatal ( NDArray N ) {
        if ( N.shape().length == 2 )
            return;

        throw new IllegalDimensionException(
            "Matrix operation: Inputs must be matrices (2 dimensions)."
        );
    }

    private static void validateSquareMatrix( NDArray N ) {
        validateMatrixFatal( N );
        int[] shape = N.shape();
        if ( shape[0] == shape[1] )
            return;

        throw new IllegalDimensionException(
            "Input matrix must be square (NxN)."
        );
    }
}
