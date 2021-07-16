package com.evanstella.jnp.core;

public abstract class NDArray {

    protected int[] shape;

    public abstract NDArray get ( );

    public abstract NDArray reshape ( );

    public abstract NDArray diagonal ( );

    public abstract NDArray transpose ( );

    public abstract NDArray flatten ( );

    public abstract NDArray copy ( );

    public int sub2ind ( int ...sub ) {
        if ( sub.length != shape.length )
            throw new IllegalDimensionException(
                "Number of subscript dimensions must match data dimensions."
            );

        int ind = 0, prod = 1;
        for ( int i = sub.length-1; i >= 0; --i ) {
            try {
                ind += prod * sub[i];
                prod *= shape[i];
            } catch ( ArrayIndexOutOfBoundsException E ) {
                throw new IllegalDimensionException(
                    "Subscript " + sub[i] + " out of bounds for dimension " +
                        i + " of length " + shape[i]
                );
            }
        }
        return ind;
    }

    public int[] ind2sub ( int ind ) {
        return null;
    }

    public int[] shape ( ) {
        int[] ret = new int[shape.length];
        System.arraycopy( shape, 0, ret, 0, shape.length );
        return ret;
    }

    public String toString ( ) {
        StringBuilder s = new StringBuilder();
        for ( int i = 0; i < shape.length-1; i++ ) {
            s.append(shape[i]).append("x");
        }
        s.append(shape[shape.length - 1]).append(" Array");
        return s.toString();
    }

}
