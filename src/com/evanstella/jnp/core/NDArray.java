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

package com.evanstella.jnp.core;

/******************************************************************************
 * <p>An abstract class encapsulating the logic and behavior of an N-dimensional
 * array. The reason this is done instead of using generics is primarily speed
 * and memory overhead: an array of primitive doubles is significantly more
 * performant than that of Java's class Double, and Java does not support
 * primitive generics for obvious reasons.
 *
 * @author Evan Stella
 *****************************************************************************/
public abstract class NDArray {

    protected int[] shape;


    /**************************************************************************
     * <p>Reshapes the NDArray by simply changing the dimensions, not the data.
     * This means the output will still be in row-major order and the
     * dimensions must satisfy the requirement that the number of elements
     * must not change.
     *************************************************************************/
    public abstract NDArray reshape ( int ...dimensions );

    /**************************************************************************
     * <p>Transposes the NDArray. Only defined for Matrices (<3 dimensions); any
     * input of a higher dimension than 2 will cause an exception
     *************************************************************************/
    public abstract NDArray transpose ( );

    /**************************************************************************
     * <p>Slices the data into a sub NDArray
     *************************************************************************/
    public abstract NDArray slice ( int[] ...dimensions );

    /**************************************************************************
     * <p>flattens the NDArray to one dimension without changing any of the data.
     *************************************************************************/
    public abstract NDArray flatten ( );

    /**************************************************************************
     * <p>Creates a deep copy of this NDArray object and returns a reference to
     * the copy.
     *************************************************************************/
    public abstract NDArray copy ( );

    /**************************************************************************
     * <p>Returns an NDArray with all dimensions of length 1 removed. If the
     * resulting data is one dimensional, the resulting object will be a row
     * vector (1xN) if the first dimension (row) was length 1, and a column
     * vector (Nx1) otherwise.
     *************************************************************************/
    public abstract NDArray squeeze ( );

    /**************************************************************************
     * <p>Resizes the array containing the data while keeping data in its
     * current subscripted position.
     *
     * Should be done infrequently as this requires copying each element to a
     * new array. There is also a little overhead associated with calculating
     * the new element indices. Elements in newly allocated space are
     * initialized to 0
     *************************************************************************/
    public abstract void resize ( int ...dimensions );

    /**************************************************************************
     * <p>Determine if NDArray is scalar (1x1)
     *************************************************************************/
    public abstract boolean isScalar ( );

    /**************************************************************************
     * <p>Determine if NDArray is a vector (1xN), (Nx1)
     *************************************************************************/
    public abstract boolean isVector ( );

    /**************************************************************************
     * <p>Returns the size of the NDArray in values.
     *************************************************************************/
    public int getSize ( ) {
        int size = 1;
        for ( int n : shape )
            size *= n;
        return size;
    }

    /**************************************************************************
     * <p>Converts the inputted subscript to a linear index for this NDArray
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
     * <p>Converts the inputted linear index to a subscript for this NDArray
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
     * <p>Converts the inputted subscript to a linear index an NDArray with the
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
     * <p>Converts the inputted linear index to a subscript for an NDArray with
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
     * <p>Determine if the dimensions of A1 and A2 are the same.
     *
     * @param A1    The first NDArray to compare
     * @param A2    The second NDArray to compare
     *
     * @return whether the dimensions of A1 = those of A2
     *************************************************************************/
    public static boolean dimensionsMatch ( NDArray A1, NDArray A2 ) {
        if ( A1.shape.length != A2.shape.length )
            return false;

        for ( int i = 0; i < A1.shape.length; i++ ) {
            if ( A1.shape[i] != A2.shape[i] )
                return false;
        }

        return true;
    }

    /**************************************************************************
     * <p>Returns the shape of the NDArray
     *
     * @return a copy of this.shape
     *************************************************************************/
    public int[] shape ( ) {
        int[] ret = new int[shape.length];
        System.arraycopy( shape, 0, ret, 0, shape.length );
        return ret;
    }

    /**************************************************************************
     * <p>Overrides Object.toString()
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
