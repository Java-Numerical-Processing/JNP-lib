package com.evanstella.jnp.math;

import com.evanstella.jnp.core.IllegalDimensionException;
import com.evanstella.jnp.core.Numeric;

public final class Matrix {

    // no instances for you
    private Matrix ( ) {}







    private static void validateMatrixFatal ( Numeric N ) {
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
