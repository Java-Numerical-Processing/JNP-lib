package com.evanstella.jnp.core;

public abstract class NDArray {

    protected int[] shape;


    public abstract NDArray get ( );

    public abstract NDArray reshape ( int ...dimensions );

    public abstract NDArray transpose ( );

    public abstract NDArray flatten ( );

    public abstract NDArray copy ( );

    public abstract NDArray squeeze ( );

    public abstract void resize ( int ...dimensions );

    public int sub2ind ( int ...sub ) {

        if ( sub.length != shape.length )
            throw new IllegalDimensionException(
                "Number of subscripts must match data dimensions."
            );

        int ind = 0, prod = 1;
        for ( int i = sub.length-1; i >= 0; --i ) {
            ind += prod * sub[i];
            prod *= shape[i];
        }
        return ind;
    }

    public int[] ind2sub ( int ind ) {
        if ( ind < 0 )
            throw new IllegalDimensionException(
                "Index can not be less than 0."
            );

        int[] sub = new int[shape.length];
        int s, prod = 1;
        for ( int i = shape.length-1; i >= 0; --i ) {
            s = ind % shape[i];
            ind -= s;
            ind /= shape[i];
            sub[i] = s;
        }

        return sub;
    };

    public static int sub2ind ( int[] dimensions, int[] sub ) {

        if ( sub.length != dimensions.length )
            throw new IllegalDimensionException(
                    "Number of subscripts must match data dimensions."
            );

        int ind = 0, prod = 1;
        for ( int i = sub.length-1; i >= 0; --i ) {
            try {
                ind += prod * sub[i];
                prod *= dimensions[i];
            } catch ( ArrayIndexOutOfBoundsException E ) {
                throw new IllegalDimensionException(
                        "Subscript " + sub[i] + " out of bounds for dimension " +
                                i + " of length " + dimensions[i]
                );
            }
        }
        return ind;
    }

    public static int[] ind2sub ( int[] dimensions, int ind ) {
        if ( ind < 0 )
            throw new IllegalDimensionException(
                    "Index can not be less than 0."
            );

        int[] sub = new int[dimensions.length];
        int s, prod = 1;
        for ( int i = dimensions.length-1; i >= 0; --i ) {
            s = ind % dimensions[i];
            ind -= s;
            ind /= dimensions[i];
            sub[i] = s;
        }

        return sub;
    };

    public int[] shape ( ) {
        int[] ret = new int[shape.length];
        System.arraycopy( shape, 0, ret, 0, shape.length );
        return ret;
    }

    public String toString ( ) {
        StringBuilder s = new StringBuilder();
        if ( shape.length == 1 )
            s.append("1x");
        for ( int i = 0; i < shape.length-1; i++ ) {
            s.append(shape[i]).append("x");
        }
        s.append(shape[shape.length - 1]).append(" Array");
        return s.toString();
    }

}
