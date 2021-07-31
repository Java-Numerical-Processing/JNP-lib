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

    public static void validateMatrixFatal ( Numeric N ) {
        if ( N.shape().length == 2 )
            return;

        throw new IllegalDimensionException(
            "Matrix operation: Inputs must be matrices (2 dimensions)."
        );
    }

    private static void validateSquareMatrix( Numeric N ) {
        validateMatrixFatal( N );
        int[] shape = N.shape();
        if ( shape[0] == shape[1] )
            return;

        throw new IllegalDimensionException(
            "Input matrix must be square (NxN)."
        );
    }
}
