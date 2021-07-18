package com.evanstella.jnp.core;

public abstract class NDArray {

    protected int[] shape;


    /**************************************************************************
     * Reshapes the NDArray by simply changing the dimensions, not the data.
     * This means the output will still be in row-major order and the
     * dimensions must satisfy the requirement that the number of elements
     * must not change.
     *************************************************************************/
    public abstract NDArray reshape ( int ...dimensions );

    /**************************************************************************
     * Transposes the NDArray. Only defined for Matrices (<3 dimensions); any
     * input of a higher dimension than 2 will cause an exception
     *************************************************************************/
    public abstract NDArray transpose ( );

    /**************************************************************************
     * Slices the data into a sub NDArray
     *
     * @param dimensions    An int[] for each dimension, e.i. if the data is 3
     *                      dimensional, inputs would be int[],int[],int[]
     *                      (or int[3][]). Each int[] should have either one
     *                      value (int[]{n}) specifying the nth index of the
     *                      dimension, or two values specifying a range of
     *                      values (int[]{0,4}: inclusive, exclusive).
     *
     *                      Example: for a 5x5 matrix:
     *                      slice(new int{0,5}, new int{2}) returns the 3rd
     *                      column of the matrix.
     *
     * @return a reference to the sliced NDArray.
     *************************************************************************/
    public abstract NDArray slice ( int[] ...dimensions );

    /**************************************************************************
     * flattens the NDArray to one dimension without changing any of the data.
     *************************************************************************/
    public abstract NDArray flatten ( );

    /**************************************************************************
     * Creates a deep copy of this NDArray object and returns a reference to
     * the copy.
     *************************************************************************/
    public abstract NDArray copy ( );

    /**************************************************************************
     * Returns an NDArray with all dimensions of length 1 removed. If the
     * resulting data is one dimensional, the resulting object will be a row
     * vector (1xN) if the first dimension (row) was length 1, and a column
     * vector (Nx1) otherwise.
     *
     * @return a reference to the squeezed NDArray.
     *************************************************************************/
    public abstract NDArray squeeze ( );

    /**************************************************************************
     * Resizes the array containing the data while keeping data in its
     * current subscripted position.
     *
     * Should be done infrequently as this requires copying each element to a
     * new array. There is also a little overhead associated with calculating
     * the new element indices. Elements in newly allocated space are
     * initialized to 0
     *************************************************************************/
    public abstract void resize ( int ...dimensions );

    /**************************************************************************
     * Converts the inputted subscript to a linear index for this NDArray
     *
     * @param sub the subscript to calculate the index for
     *
     * @return the calculated linear index
     *************************************************************************/
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

    /**************************************************************************
     * Converts the inputted linear index to a subscript for this NDArray
     *
     * @param ind the index to calculate the subscript for
     *
     * @return the calculated subscript
     *************************************************************************/
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

    /**************************************************************************
     * Converts the inputted subscript to a linear index an NDArray with the
     * inputted dimensions
     *
     * @param dimensions    the dimensions to calculate for
     * @param sub           the subscript to calculate the index for
     *
     * @return the calculated index
     *************************************************************************/
    public static int sub2ind ( int[] dimensions, int[] sub ) {

        if ( sub.length != dimensions.length )
            throw new IllegalDimensionException(
                "Number of subscripts must match data dimensions."
            );

        int ind = 0, prod = 1;
        for ( int i = sub.length-1; i >= 0; --i ) {
            ind += prod * sub[i];
            prod *= dimensions[i];
        }
        return ind;
    }

    /**************************************************************************
     * Converts the inputted linear index to a subscript for an NDArray with
     * the inputted dimensions
     *
     * @param dimensions    the dimensions to calculate for
     * @param ind           the index to calculate the subscript for
     *
     * @return the calculated subscript
     *************************************************************************/
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

    /**************************************************************************
     * Returns the shape of the NDArray
     *
     * @return a copy of this.shape
     *************************************************************************/
    public int[] shape ( ) {
        int[] ret = new int[shape.length];
        System.arraycopy( shape, 0, ret, 0, shape.length );
        return ret;
    }

    /**************************************************************************
     * Overrides Object.toString()
     *************************************************************************/
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


    /**************************************************************************
     *                          Internal Methods
     *************************************************************************/

    /* sub2ind with no dimension checking */
    protected static int sub2indNoCheck ( int[] dimensions, int[] sub ) {
        int ind = 0, prod = 1;
        for ( int i = sub.length-1; i >= 0; --i ) {
            ind += prod * sub[i];
            prod *= dimensions[i];
        }
        return ind;
    }

    /* ind2sub with no dimension checking */
    protected static int[] ind2subNoCheck ( int[] dimensions, int ind ) {
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

}
