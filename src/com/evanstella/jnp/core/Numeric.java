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

import com.evanstella.jnp.math.Element;

/******************************************************************************
 * An NDArray of numeric values. Data is stored internally as a 1 dimensional
 * array of doubles.
 *
 * @author Evan Stella
 *****************************************************************************/
public class Numeric extends NDArray {

    protected double[] data;


    /**************************************************************************
     * <p>Class constructor. Initializes a Numeric with parameterized dimensions.
     * The dimensions are deep copied to limit their write access to the object.
     *
     * @param dimensions The arbitrary dimensions for the N-D Numeric.
     *************************************************************************/
    public Numeric ( int ...dimensions ) {
        shape = new int[dimensions.length];
        System.arraycopy(dimensions, 0, shape,0, dimensions.length);
        int size = 1;
        for ( int n : dimensions )
            size *= n;
        if ( size < 1 )
            throw new IllegalDimensionException(
                "Dimension lengths must be positive and non-zero."
            );
        data = new double[size];
    }


    /**************************************************************************
     * <p>Class constructor. Initialize a Numeric from a double[]. Only a
     * shallow copy of the array is made
     *
     * @param indata The array to initialize from
     *************************************************************************/
    public Numeric ( double[] indata ) {
        data = indata;
        shape = new int[]{ indata.length };
    }

    /**************************************************************************
     * <p>Class constructor. Initialize a Numeric from a double[][]. Makes a
     * a deep copy in order to aggregate the double[][] into a double[].
     *
     * @param indata The array to initialize from
     *************************************************************************/
    public Numeric ( double[][] indata ) {
        data = new double[indata.length * indata[0].length];
        int ind = 0;
        for ( double[] row : indata ) {
            if ( row.length != indata[0].length )
                throw new IllegalDimensionException(
                        "Input dimensions must be consistent."
                );
            for ( double d : row )
                data[ind++] = d;
        }
        shape = new int[]{ indata.length, indata[0].length };
    }

    /**************************************************************************
     * <p>Initializes a Logical with all true values.
     *
     * @param dimensions The dimensions for the N-D Logical.
     *
     * @return a reference the initialized Logical
     *************************************************************************/
    public static Numeric Ones ( int ...dimensions ) {
        Numeric N = new Numeric( dimensions );
        java.util.Arrays.fill(N.data, 1.0);
        return N;
    }

    /**************************************************************************
     * <p>Initializes a Numeric with all zeros. Literally the same as
     * calling the constructor but here for completeness' sake.
     *
     * @param dimensions the dimensions for the N-D Numeric.
     *
     * @return a reference the initialized Numeric
     *************************************************************************/
    public static Numeric Zeros ( int ...dimensions ) {
        return new Numeric( dimensions );
    }

    /**************************************************************************
     * <p>Initializes a 1x1 scalar Numeric with the inputted real data.
     *
     * @param real   the real scalar
     *
     * @return a reference the initialized Numeric
     *************************************************************************/
    public static Numeric Scalar ( double real ) {
        return new Numeric( new double[]{real} );
    }


    /**************************************************************************
     * <p>Initializes a Numeric with random real values using java.util.Random.
     *
     * @param dimensions the dimensions for the N-D Numeric.
     *
     * @return a reference the initialized Numeric
     *************************************************************************/
    public static Numeric Rand ( int ...dimensions ) {
        java.util.Random R = new java.util.Random( );
        Numeric N = new Numeric( dimensions );

        for ( int i = 0; i < N.data.length; i++ )
            N.data[i] = R.nextDouble();

        return N;
    }

    /**************************************************************************
     * <p>Initializes a Numeric with random real values with a seed using
     * java.util.Random.
     *
     * @param seed       the PRNG seed.
     * @param dimensions the dimensions for the N-D Numeric.
     *
     * @return a reference the initialized Numeric
     *************************************************************************/
    public static Numeric Rand ( long seed,  int ...dimensions ) {
        java.util.Random R = new java.util.Random( seed );
        Numeric N = new Numeric( dimensions );

        for ( int i = 0; i < N.data.length; i++ )
            N.data[i] = R.nextDouble();

        return N;
    }

    /**************************************************************************
     * <p>Initializes a row vector of linearly spaced elements.
     *
     * @param start     the value of the first element
     * @param end       the value of the last element
     * @param numPts    the number of points to sample
     *
     * @return A Numeric of linearly space elements.
     *************************************************************************/
    public static Numeric LinSpace ( double start, double end, int numPts ) {
        Numeric linspace = new Numeric( 1, numPts );
        double[] data = linspace.data;

        double step = ( end - start ) / (numPts-1);
        double val = start;
        int ind = 0;
        while ( ind < numPts ) {
            data[ind++] = val;
            val += step;
        }
        return linspace;
    }

    /**************************************************************************
     * <p>Initializes a row vector of logarithmically (base 10) spaced elements.
     *
     * @param start     the value of the first element
     * @param end       the value of the last element
     * @param numPts    the number of points to sample
     *
     * @return A Numeric of logarithmically space elements.
     *************************************************************************/
    public static Numeric LogSpace ( double start, double end, int numPts ) {
        return Element.pow( Numeric.Scalar(10), LinSpace( start, end, numPts) );
    }

    /**************************************************************************
     * <p>Initializes a row vector of logarithmically spaced elements.
     *
     * @param b         the base of the log scale
     * @param s         the value of the first element
     * @param end       the value of the last element
     * @param pts       the number of points to sample
     *
     * @return A Numeric of logarithmically space elements.
     *************************************************************************/
    public static Numeric LogSpace ( double b, double s, double end, int pts ) {
        return Element.pow( Numeric.Scalar(b), LinSpace( s, end, pts) );
    }

    /**************************************************************************
     * <p>Initializes an identity matrix with the inputted number of rows and
     * columns
     *
     * @param r   The number of rows for the matrix
     * @param c   The number of columns for the matrix
     *
     * @return a reference to the initialized matrix
     *************************************************************************/
    public static Numeric Eye ( int r, int c ) {
        Numeric eye = new Numeric( r, c );
        for ( int i = 0; i < r && i < c; i++ ) {
            eye.data[i * c + i] = 1;
        }
        return eye;
    }

    /**************************************************************************
     * <p>Initializes a square matrix with the diagonal equal to the inputted
     * vector. All off-diagonal values are initialized to zero.
     *
     * @param N     A Numeric vector ( 1 dimensional row or column vector )
     *
     * @return a reference to the initialized matrix
     *************************************************************************/
    public static Numeric Diagonal ( Numeric N ) {
        if ( N.shape.length > 2 )
            throw new IllegalDimensionException(
            "Input must be a row or column vector (shape = N or 1xN or Nx1)."
        );
        if ( N.shape.length != 1 && N.shape[0] != 1 && N.shape[1] != 1 )
            throw new IllegalDimensionException(
            "Input must be a row or column vector (shape = N or 1xN or Nx1)."
        );

        int size = N.data.length;
        Numeric result = new Numeric( size, size );
        for ( int ind, i = 0; i < size; i++ ) {
            ind = i * size + i;
            result.data[ind] = N.data[i];
        }
        return result;
    }

    /**************************************************************************
     * <p>Check whether this Numeric is a scalar (size of 1)
     *
     * @return true if this Numeric is scalar
     *************************************************************************/
    public boolean isScalar ( ) {
        return data.length == 1;
    }

    /**************************************************************************
     * <p>Check whether this Numeric is a scalar (size of 1)
     *
     * @return true if this Numeric is scalar
     *************************************************************************/
    public boolean isVector ( ) {
        if ( shape.length == 1 )
            return true;
        return shape.length == 2 && (shape[0] == 1 || shape[1] == 1);
    }

    /**************************************************************************
     * <p>If this numeric is a scalar, return the scalar real value
     *
     * @return true if this Numeric is scalar
     *************************************************************************/
    public double value ( ) {
        if ( !isScalar() )
            throw new IllegalDimensionException(
                "Error using value(): Numeric must be scalar."
            );
        return data[0];
    }

    /**************************************************************************
     * <p>Gets a reference to the raw real data contained in the Numeric. Only
     * recommended if you need fast access to the data in the object.
     *
     * @return a reference to the Logical data.
     *************************************************************************/
    public double[] getData ( ) {
        return data;
    }

    /**************************************************************************
     * <p>Sets the value of the data at the subscript with the inputted value
     * WITHOUT resizing if the subscript is out of bounds of the data
     * dimensions. If this is the case, an IllegalDimensionException will be
     * thrown.
     *
     * @param indata    the data to set (real component)
     * @param sub           the subscript at which to set the data
     *************************************************************************/
    public void set ( double indata, int[] sub ) {
        int ind = sub2ind( shape, sub );
        if ( ind >= data.length ) {
            throw new IllegalDimensionException(
                "Subscript out of bounds for dimensions"
            );
        }
        data[ind] = indata;
    }

    /**************************************************************************
     * <p>Sets the value of the data at the subscript with the inputted value
     * WITHOUT resizing if the subscript is out of bounds of the data
     * dimensions.
     *
     * @param indata    the data to set (real component)
     * @param ind           the index at which to set the data
     *************************************************************************/
    public void set ( double indata, int ind ) {
        data[ind] = indata;
    }

    /**************************************************************************
     * <p>Sets the value of the data at the subscript with the inputted value. If
     * the subscript is out of bounds of the data dimensions, the data will be
     * resized. This can be slower on larger data sets as this
     * this operation is O(n) for data of size n.
     *
     * @param indata    the data to set (real component)
     * @param sub           the subscript at which to set the data
     *************************************************************************/
    public void setAt ( double indata, int[] sub ) {
        int ind = sub2ind( sub );
        if ( ind >= data.length ) {
            int[] newDims = new int[shape.length];
            for ( int i = 0; i < shape.length; i++ ) {
                if ( shape[i] > sub[i] )
                    newDims[i] = shape[i];
                else
                    newDims[i] = sub[i] + 1;
            }
            resizeNoCheck( newDims );
            ind = sub2ind( sub );
        }
        data[ind] = indata;
    }

    /**************************************************************************
     * <p>Indexes the Numeric using the inputted Logical as a mask and sets the
     * indexed values to val.
     *
     * @param val       The real value to set the indexed elements to
     * @param Inds      Logical to be used as an index into the data
     *************************************************************************/
    public void set ( double val, Logical Inds ) {
        if ( Inds.shape.length != shape.length)
            throw new IllegalDimensionException(
                "Index dimensions must be equal to data dimensions."
        );
        for ( int i = 0; i < shape.length; i++ ) {
            if ( Inds.shape[i] != shape[i] )
                throw new IllegalDimensionException(
                    "Index dimensions must be equal to data dimensions."
            );
        }

        for ( int i = 0; i < data.length; i++ ) {
            if ( Inds.data[i] ) {
                data[i] = val;
            }
        }
    }

    /**************************************************************************
     * <p>Indexes the Numeric using the inputted Logical as a mask and sets the
     * indexed values to val + current value.
     *
     * @param val   The real component of the value to set the indexed
     *                  elements to
     * @param Inds      Logical to be used as an index into the data
     *************************************************************************/
    public void setAdd ( double val, Logical Inds ) {
        if ( !NDArray.dimensionsMatch( this, Inds ) )
            throw new IllegalDimensionException(
                    "Index dimensions must be equal to data dimensions."
            );

        for ( int i = 0; i < data.length; i++ ) {
            if ( Inds.data[i] ) {
                data[i] += val;
            }
        }
    }

    /**************************************************************************
     * <p>Indexes the Numeric using the inputted Logical as a mask and sets the
     * indexed values to val * current value.
     *
     * @param val The real component of the value to set the indexed
     *                  elements to
     * @param Inds      Logical to be used as an index into the data
     *************************************************************************/
    public void setMul ( double val, Logical Inds ) {
        if ( !NDArray.dimensionsMatch( this, Inds ) )
            throw new IllegalDimensionException(
                    "Index dimensions must be equal to data dimensions."
            );

        for ( int i = 0; i < data.length; i++ ) {
            if ( Inds.data[i] ) {
                data[i] = data[i] * val;
            }
        }
    }

    /**************************************************************************
     * <p>Indexes the Numeric using the inputted Logical as a mask. Returns a row
     * vector with the indexed values.
     *
     * @param Inds      Logical to be used as an index into the data
     *
     * @return a reference to the initialized vector.
     *************************************************************************/
    public Numeric get ( Logical Inds ) {

        if ( !NDArray.dimensionsMatch( this, Inds ) )
            throw new IllegalDimensionException(
                "Index dimensions must be equal to data dimensions."
        );
        // get number of elements
        int numElements = 0;
        for ( int i = 0; i < data.length; i++ )
            if ( Inds.data[i] ) numElements++;

        Numeric result = new Numeric( 1, numElements );
        for ( int ind = 0, i = 0; i < data.length; i++ ) {
            if ( Inds.data[i] ) {
                result.data[ind] = data[i];
            }
        }
        return result;
    }

    /**************************************************************************
     * <p>Returns a scalar numeric with the value at the inputted subscript.
     *
     * @param sub   the subscript to retrieve data at
     *
     * @return a reference to the initialized scalar.
     *************************************************************************/
    public Numeric get ( int... sub ) {
        int ind = sub2ind( shape, sub );
        if ( ind >= data.length )
            throw new IllegalDimensionException(
                "Subscript out of range of data dimensions");

        Numeric ret = new Numeric( 1,1 );
        ret.data[0] = data[ind];
        return ret;
    }

    /**************************************************************************
     * <p>Returns a scalar numeric with the value at the inputted linear
     * index.
     *
     * @param ind   the linear index to retrieve data at
     *
     * @return a reference to the initialized scalar.
     *************************************************************************/
    public Numeric get ( int ind ) {
        if ( ind >= data.length )
            throw new IllegalDimensionException(
                "Subscript out of range of data dimensions");

        Numeric ret = new Numeric( 1,1 );
        ret.data[0] = data[ind];
        return ret;
    }

    /**************************************************************************
     * <p>Find the index of the max value in the data. Only compares the real
     * component of the data since comparisons between complex numbers are not
     * well defined.
     *
     * @return the linear index of the max value in the data
     *************************************************************************/
    public int findMax ( ) {
        double max = data[0];
        int ind = 0;
        for ( int i = 0; i < data.length; i++ ) {
            if ( data[i] > max) {
                max = data[i];
                ind = i;
            }
        }
        return ind;
    }

    /**************************************************************************
     * <p>Find the index of the absolute max value in the data. Only compares the
     * real component of the data since comparisons between complex numbers are
     * not well defined.
     *
     * @return the linear index of the max value in the data
     *************************************************************************/
    public int findAbsMax ( ) {
        double max = data[0];
        double absVal, tmp;
        int ind = 0;
        for ( int i = 0; i < data.length; i++ ) {
            tmp = data[i];
            absVal = (tmp <= 0.0) ? 0.0 - tmp : tmp;
            if ( absVal > max) {
                max = absVal;
                ind = i;
            }
        }
        return ind;
    }

    /**************************************************************************
     * <p>Find the index of the min value in the data. Only compares the real
     * component of the data since comparisons between complex numbers are not
     * well defined.
     *
     * @return the linear index of the min value in the data
     *************************************************************************/
    public int findMin ( ) {
        double min = data[0];
        int ind = 0;
        for ( int i = 0; i < data.length; i++ ) {
            if ( data[i] < min) {
                min = data[i];
                ind = i;
            }
        }
        return ind;
    }

    /**************************************************************************
     * <p>Find the index of the absolute min value in the data. Only compares the
     * real component of the data since comparisons between complex numbers are
     * not well defined.
     *
     * @return the linear index of the min value in the data
     *************************************************************************/
    public int findAbsMin ( ) {
        double min = data[0];
        double absVal, tmp;
        int ind = 0;
        for ( int i = 0; i < data.length; i++ ) {
            tmp = data[i];
            absVal = (tmp <= 0.0) ? 0.0 - tmp : tmp;
            if ( absVal < min) {
                min = absVal;
                ind = i;
            }
        }
        return ind;
    }

    /**************************************************************************
     * <p>Reshapes the Numeric by simply changing the dimensions, not the data.
     * This means the output will still be in row-major order and the dimensions
     * must satisfy the requirement that the number of elements must not change.
     *
     * @param dimensions    The dimensions to shape the data to
     *
     * @return a reference to the reshaped Numeric
     *************************************************************************/
    public Numeric reshape ( int ...dimensions ) {
        int size = 1;
        for ( int n : dimensions )
            size *= n;
        if ( size != data.length )
            throw new IllegalDimensionException(
                "Invalid Dimensions. Number of elements must be the same"
            );
        Numeric reshaped = new Numeric( dimensions );
        copyData( reshaped );
        return reshaped;
    }

    /**************************************************************************
     * <p>Transposes the Numeric. Only defined for Matrices (<3 dimensions); any
     * input of a higher dimension than 2 will cause an exception
     *
     * @return a reference to the transposed data.
     *************************************************************************/
    public Numeric transpose ( ) {
        if ( shape.length > 2 )
            throw new IllegalDimensionException(
                "Transpose operation is not defined for data with more than 2 dimensions"
            );

        Numeric transposed;
        if ( shape.length == 1 )
            transposed = new Numeric( shape[0], 1 );
        else
            transposed = new Numeric( shape[1], shape[0] );

        int r = transposed.shape[0];
        int c = transposed.shape[1];
        for ( int i = 0; i < r; i++ )
            for (int j = 0; j < c; j++)
                transposed.data[i*c+j] = data[j*r+i];

        return transposed;
    }

    /**************************************************************************
     * <p>flattens the Numeric to one dimension without changing any of the data.
     *
     * @return a reference to the flattened Numeric
     *************************************************************************/
    public Numeric flatten ( ) {
        Numeric flat = new Numeric( 1, data.length );
        copyData( flat );
        return flat;
    }

    /**************************************************************************
     * <p>Creates a deep copy of this Numeric object and returns a reference to
     * the copy.
     *
     * @return a reference to the copied object
     *************************************************************************/
    public Numeric copy ( ) {
        Numeric copy = new Numeric( this.shape );
        copyData( copy );
        return copy;
    }

    /**************************************************************************
     * <p>Returns a Numeric with all dimensions of length 1 removed. If the
     * resulting data is one dimensional, the resulting object will be a row
     * vector (1xN) if the first dimension (row) was length 1, and a column
     * vector (Nx1) otherwise.
     *
     * @return a reference to the resulting Numeric.
     *************************************************************************/
    public Numeric squeeze ( ) {
        int newSize = 0;
        for ( int n : shape ) {
            if (n > 1) newSize++;
        }
        boolean makeColumn = false;
        if ( shape[0] != 1 && newSize == 1 ) {
            makeColumn = true;
            newSize++;
        }
        int[] newShape = new int[newSize];
        int ind = 0;
        for ( int n : shape ) {
            if (n > 1) newShape[ind++] = n;
        }
        if ( makeColumn )
            newShape[ind] = 1;
        Numeric squeezed = new Numeric( newShape );
        copyData( squeezed );
        return squeezed;
    }

    /**************************************************************************
     * <p>Resizes the array containing the Numeric data while keeping data in its
     * current subscripted position.
     *
     * <p>Should be done infrequently as this requires copying each element to a
     * new array. There is also a little overhead associated with calculating
     * the new element indices. Elements in newly allocated space are
     * initialized to 0.0
     *
     * @param dimensions The new dimensions for the data
     *************************************************************************/
    public void resize ( int ...dimensions ) {
        for ( int i = 0; i < shape.length; i++ )
            if ( shape[i] > dimensions[i] )
                throw new IllegalDimensionException(
                        "Resized dimensions must be greater than current dimensions"
                );
        int newSize = 1;
        for ( int n : dimensions )
            newSize *= n;
        double[] newdata = new double[newSize];
        int[] sub;
        int newInd;

        for ( int i = 0; i < data.length; i++ ) {
            sub = ind2sub( shape, i );
            newInd = sub2ind( dimensions, sub );
            newdata[newInd] = data[i];
        }
        data = newdata;
        shape = dimensions;
    }

    /**************************************************************************
     * <p>Slices the data into a sub Numeric
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
     * @return a reference to the sliced Numeric.
     *************************************************************************/
    public Numeric slice ( int[] ...dimensions ) {
        if ( dimensions.length != shape.length )
            throw new IllegalDimensionException(
                "Number of slice dimensions must match number of data dimensions");
        int[] newShape = new int[shape.length];
        int size = 1;
        for ( int i = 0; i < shape.length; i++ ) {
            int[] dim = dimensions[i];
            if ( dim.length == 1 )
                newShape[i] = 1;
            else if ( dim.length == 2 ) {
                if ( dim[1] > shape[i] )
                    throw new IllegalDimensionException(
                        "Slice index " + dim[1] + " out of bounds of dimension " + i);
                if ( dim[0] >= dim[1] )
                    throw new IllegalDimensionException(
                        "Slice index 1 must be greater than index 2 for dimension " + i);
                newShape[i] = dim[1] - dim[0];
                size *= (dim[1] - dim[0]);
            }
            else
                throw new IllegalDimensionException(
                    "Slice indices must be a single index or two values");
            if ( newShape[i] < 1 )
                throw new IllegalDimensionException(
                    "Slice indices must be positive");
        }

        double[] newData = sliceData( dimensions, newShape, size );

        Numeric sliced = new Numeric( newShape );
        sliced.data = newData;
        return sliced;
    }

    /**************************************************************************
     * <p>Concatenate the inputted Numeric to the end of this Numeric. Note that
     * the number of dimensions and the dimension lengths, except for the
     * dimension being concatenated along, must be equal.
     *
     * @param dimension the dimension to concatenate along
     * @param N         the Numeric to concatenate
     *
     * @return a reference to the new Numeric
     *************************************************************************/
    public Numeric concat ( int dimension, Numeric N ) {
        if ( dimension >= shape.length )
            throw new IllegalDimensionException(
                "Logical does not have " + (dimension+1) + " dimensions.");
        if ( shape.length != N.shape.length )
            throw new IllegalDimensionException(
                "Both objects must have the same number of dimensions.");
        for ( int i = 0; i < shape.length; i++ ) {
            if ( i != dimension && shape[i] != N.shape[i] )
                throw new IllegalDimensionException(
                "All dimensions but the concat dimension must be equal in length.");
        }

        Numeric catted = this.copy( );
        int[] newShape = catted.shape( );
        newShape[dimension] = shape[dimension] + N.shape[dimension];
        catted.resizeNoCheck( newShape );

        int[] sub;
        int ind;
        for ( int i = 0; i < N.data.length; i++ ) {
            sub = ind2subNoCheck( N.shape, i );
            sub[dimension] += shape[dimension];
            ind = sub2indNoCheck( newShape, sub );
            catted.data[ind] = N.data[i];
        }

        return catted;
    }

    /**************************************************************************
     * <p>Creates a string representation of the Numeric for printing. Will only
     * show actual data for 1 and 2 dimensional data as higher dimensional
     * data is difficult to display well in a string.
     *
     * @return a string representation of the Logical
     *************************************************************************/
    public String toString ( ) {
        int maxInd = findAbsMax( );
        double max = data[maxInd];
        int L = String.format( "%.4f", max ).length();

        StringBuilder s = new StringBuilder(super.toString() + " <Numeric>\n");
        if ( shape.length < 3 ) {
            int row, col;
            if ( shape.length == 2 ) {
                row = shape[0];
                col = shape[1];
            } else {
                row = 1;
                col = shape[0];
            }
            for ( int i = 0; i < row; i++ ) {
                if ( data[i * col] < 0 )
                    s.append("[ ");
                else
                    s.append("[  ");
                for ( int j = 0; j < col; j++ ) {
                    String tmp;
                    double real = data[i * col + j];
                    if ( real < 0 )
                        tmp = String.format("%" + (L+1) + ".4f", data[i * col + j]);
                    else
                        tmp = String.format("%" + L + ".4f", data[i * col + j]);

                    if ( j < col-1 && data[i * col + j+1] < 0 )
                        s.append(tmp).append(" ");
                    else
                        s.append(tmp).append("  ");
                }
                s.append("]\n");
            }
        }

        return s.toString();
    }

    /**************************************************************************
     * <p>Compares two Numeric objects to check if they are equal in both
     * dimension and data. NOT RECOMMENDED when dealing with arithmetic because
     * of the possibility of precision issues. Use equalTolerance() instead.
     *
     * @param N the other Numeric to compare this one to
     *
     * @return  true if the two Numeric objects are equal, false otherwise.
     *************************************************************************/
    public boolean equals ( Numeric N ) {
        if ( this == N )
            return true;
        if ( !NDArray.dimensionsMatch( this, N ) )
            return false;

        for ( int i = 0; i < data.length; i++ ) {
            if ( this.data[i] != N.data[i] )
                return false;
        }

        return true;
    }

    /**************************************************************************
     * <p>Compares two Numeric objects to check if they are equal with a
     * tolerance
     *
     * @param N         the other Numeric to compare to
     * @param tolerance the tolerance for the difference between two
     *                  corresponding data in the two Numerics
     *
     * @return true if the two Numerics are equal to the specified tolerance
     *************************************************************************/
    public boolean equalsTolerance ( Numeric N, double tolerance ) {
        if ( this == N )
            return true;
        if ( !NDArray.dimensionsMatch( this, N ) )
            return false;

        double diff;
        for ( int i = 0; i < data.length; i++ ) {
            diff = this.data[i] - N.data[i];
            diff = (diff <= 0.0) ? 0.0 - diff : diff; // get abs value
            if ( diff > tolerance )
                return false;
        }

        return true;
    }


    /**************************************************************************
     *                          Internal Methods
     *************************************************************************/

    /* Resize data without checking dimensions. */
    protected void resizeNoCheck ( int ...dimensions ) {
        int newSize = 1;
        for ( int n : dimensions )
            newSize *= n;

        double[] newdata = new double[newSize];
        int[] sub;
        int newInd;
        for ( int i = 0; i < data.length; i++ ) {
            sub = ind2sub( shape, i );
            newInd = sub2ind( dimensions, sub );
            newdata[newInd] = data[i];
        }
        data = newdata;
        shape = dimensions;
    }

    /* Moved this to its own function bc slice was getting a little long */
    protected double[] sliceData ( int[][] dims, int[] newShape, int size ) {
        // get first subscript in the slice
        int[] sub = new int[shape.length];
        int ind = 0;
        for ( int[] dim : dims )
            sub[ind++] = dim[0];

        // get the last subscript in the slice
        int[] subLast = new int[shape.length];
        ind = 0;
        for ( int[] dim : dims ) {
            if ( dim.length == 1 )
                subLast[ind++] = dim[0];
            else
                subLast[ind++] = dim[1]-1;
        }

        // start at first subscript, get all elements that fall in the slice
        double[] newdata = new double[size];
        int startInd = sub2indNoCheck( shape, sub );
        int endInd = sub2indNoCheck( shape, subLast );
        ind = 0;
        newdata[ind++] = data[startInd];
        for ( int i = startInd+1; i < endInd; i++ ) {
            sub = ind2subNoCheck( shape, i );
            boolean add = true;
            for ( int j = 0; j < shape.length; j++ ) {
                if ( dims[j].length == 1 && sub[j] != dims[j][0] ) {
                    add = false;
                    break;
                } else if ( dims[j].length == 2 &&
                        (sub[j] < dims[j][0] || sub[j] >= dims[j][1]) ) {
                    add = false;
                    break;
                }
            }
            if ( add )
                newdata[ind++] = data[i];
        }
        if ( startInd != endInd )
            newdata[ind] = data[endInd];

        return newdata;
    }

    /* copy this objects data into dest */
    protected void copyData ( Numeric dest ) {
        System.arraycopy(data, 0, dest.data, 0, data.length );
    }

}
