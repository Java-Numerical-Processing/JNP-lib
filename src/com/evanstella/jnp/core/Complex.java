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
 * An NDArray of Complex values. Complex numbers are stored as two doubles,
 * one storing the real component and the other the imaginary component. This
 * means a Complex NDArray is twice as large in memory as a regular Numeric
 * NDArray holding the same number of values.
 *
 * @author Evan Stella
 *****************************************************************************/
public class Complex extends NDArray {

    protected double[] dataReal;
    protected double[] dataImag;


    /**************************************************************************
     * <p>Class constructor. Initializes a Complex with parameterized dimensions.
     * The dimensions are deep copied to limit their write access to the object.
     *
     * @param dimensions The arbitrary dimensions for the N-D Complex.
     *************************************************************************/
    public Complex ( int ...dimensions ) {
        shape = new int[dimensions.length];
        System.arraycopy(dimensions, 0, shape,0, dimensions.length);
        int size = 1;
        for ( int n : dimensions )
            size *= n;
        if ( size < 1 )
            throw new IllegalDimensionException(
                "Dimension lengths must be positive and non-zero."
            );
        dataReal = new double[size];
        dataImag = new double[size];
    }

    /**************************************************************************
     * <p>Class constructor. Initializes a Complex from a Numeric.
     *
     * @param N The Numeric to initialize from.
     *************************************************************************/
    public Complex ( Numeric N ) {
        shape = new int[N.shape.length];
        System.arraycopy(N.shape, 0, shape,0, N.shape.length);
        int size = N.data.length;
        if ( size < 1 )
            throw new IllegalDimensionException(
                "Dimension lengths must be positive and non-zero."
            );
        dataReal = new double[size];
        dataImag = new double[size];
        System.arraycopy(N.data, 0, dataReal,0, N.data.length);
    }


    /**************************************************************************
     * <p>Class constructor. Initialize a Complex from a double[]. Only a
     * shallow copy of the array is made
     *
     * @param indata The array to initialize from
     *************************************************************************/
    public Complex ( double[] indata ) {
        dataReal = indata;
        dataImag = new double[indata.length];
        shape = new int[]{ indata.length };
    }

    /**************************************************************************
     * <p>Class constructor. Initialize a Complex from a double[][]. Makes a
     * a deep copy in order to aggregate the double[][] into a double[].
     *
     * @param indata The array to initialize from
     *************************************************************************/
    public Complex ( double[][] indata ) {
        dataReal = new double[indata.length * indata[0].length];
        dataImag = new double[indata.length * indata[0].length];
        int ind = 0;
        for ( double[] row : indata ) {
            if ( row.length != indata[0].length )
                throw new IllegalDimensionException(
                        "Input dimensions must be consistent."
                );
            for ( double d : row )
                dataReal[ind++] = d;
        }
        shape = new int[]{ indata.length, indata[0].length };
    }

    /**************************************************************************
     * <p>Class constructor. Initialize a Complex from a double[]. Only a
     * shallow copy of the array is made
     *
     * @param inDataReal The array to initialize the real data component from
     * @param inDataImag The array to initialize the imaginary data component
     *                   from
     *************************************************************************/
    public Complex ( double[] inDataReal, double[] inDataImag ) {
        if ( inDataImag.length != inDataReal.length )
            throw new IllegalDimensionException(
                    "Real and imaginary data must be the same size"
            );

        dataReal = inDataReal;
        dataImag = inDataImag;
        shape = new int[]{ inDataReal.length };
    }

    /**************************************************************************
     * <p>Class constructor. Initialize a Complex from a double[][]. Makes a
     * a deep copy in order to aggregate the double[][] into a double[].
     *
     * @param inDataReal The array to initialize the real data component from
     * @param inDataImag The array to initialize the imaginary data component
     *                   from
     *************************************************************************/
    public Complex ( double[][] inDataReal, double[][] inDataImag ) {
        if ( inDataImag.length != inDataReal.length ||
                inDataImag[0].length != inDataReal[0].length)
            throw new IllegalDimensionException(
                    "Real and imaginary data must be the same size"
            );
        dataReal = new double[inDataReal.length * inDataReal[0].length];
        dataImag = new double[inDataImag.length * inDataImag[0].length];
        int ind = 0;
        for ( int i = 0; i < inDataReal.length; i++ ) {
            if ( inDataReal[i].length != inDataReal[0].length ||
                    inDataImag[i].length != inDataImag[0].length )
                throw new IllegalDimensionException(
                        "Input dimensions must be consistent."
                );
            for ( int j = 0; j < inDataReal[i].length; i++ ) {
                dataReal[ind] = inDataReal[i][j];
                dataImag[ind] = inDataImag[i][j];
                ind++;
            }
        }
        shape = new int[]{ inDataReal.length, inDataReal[0].length };
    }

    /**************************************************************************
     * <p>Initializes a Complex with all zeros. Literally the same as
     * calling the constructor but here for completeness' sake.
     *
     * @param dimensions the dimensions for the N-D Complex.
     *
     * @return a reference the initialized Complex
     *************************************************************************/
    public static Complex Zeros ( int ...dimensions ) {
        return new Complex( dimensions );
    }

    /**************************************************************************
     * <p>Initializes a 1x1 scalar Complex with the inputted real data.
     *
     * @param real   the real scalar
     *
     * @return a reference the initialized Complex
     *************************************************************************/
    public static Complex Scalar ( double real ) {
        return new Complex( new double[]{real} );
    }

    /**************************************************************************
     * <p>Initializes a 1x1 scalar Complex with the inputted complex data.
     *
     * @param real   the real component of the scalar
     * @param imag   the imaginary component of the scalar
     *
     * @return a reference the initialized Complex
     *************************************************************************/
    public static Complex Scalar ( double real, double imag ) {
        return new Complex( new double[]{real}, new double[]{imag} );
    }


    /**************************************************************************
     * <p>Initializes a Complex with random real values using java.util.Random.
     *
     * @param dimensions the dimensions for the N-D Complex.
     *
     * @return a reference the initialized Complex
     *************************************************************************/
    public static Complex Rand ( int ...dimensions ) {
        java.util.Random R = new java.util.Random( );
        Complex C = new Complex( dimensions );
        C.dataImag = new double[C.dataReal.length];

        for ( int i = 0; i < C.dataReal.length; i++ ) {
            C.dataReal[i] = R.nextDouble();
            C.dataImag[i] = R.nextDouble();
        }

        return C;
    }

    /**************************************************************************
     * <p>Initializes a Complex with random real values with a seed using
     * java.util.Random.
     *
     * @param seed       the PRNG seed.
     * @param dimensions the dimensions for the N-D Complex.
     *
     * @return a reference the initialized Complex
     *************************************************************************/
    public static Complex Rand ( long seed,  int ...dimensions ) {
        java.util.Random R = new java.util.Random( seed );
        Complex N = new Complex( dimensions );
        N.dataImag = new double[N.dataReal.length];

        for ( int i = 0; i < N.dataReal.length; i++ ) {
            N.dataReal[i] = R.nextDouble();
            N.dataImag[i] = R.nextDouble();
        }

        return N;
    }

    /**************************************************************************
     * <p>Initializes a row vector of linearly spaced elements.
     *
     * @param start     the value of the first element
     * @param end       the value of the last element
     * @param numPts    the number of points to sample
     *
     * @return A Complex of linearly space elements.
     *************************************************************************/
    public static Complex LinSpace ( double start, double end, int numPts ) {
        Complex linspace = new Complex( 1, numPts );
        double[] data = linspace.dataReal;

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
     * @return A Complex of logarithmically space elements.
     *************************************************************************/
    public static Complex LogSpace ( double start, double end, int numPts ) {
        return Element.pow( Complex.Scalar(10), LinSpace( start, end, numPts) );
    }

    /**************************************************************************
     * <p>Initializes a row vector of logarithmically spaced elements.
     *
     * @param b         the base of the log scale
     * @param s         the value of the first element
     * @param end       the value of the last element
     * @param pts       the number of points to sample
     *
     * @return A Complex of logarithmically space elements.
     *************************************************************************/
    public static Complex LogSpace ( double b, double s, double end, int pts ) {
        return Element.pow( Complex.Scalar(b), LinSpace( s, end, pts) );
    }

    /**************************************************************************
     * <p>Initializes a square matrix with the diagonal equal to the inputted
     * vector. All off-diagonal values are initialized to zero.
     *
     * @param N     A Complex vector ( 1 dimensional row or column vector )
     *
     * @return a reference to the initialized matrix
     *************************************************************************/
    public static Complex Diagonal ( Complex N ) {
        if ( N.shape.length > 2 )
            throw new IllegalDimensionException(
                "Input must be a row or column vector (shape = N or 1xN or Nx1)."
            );
        if ( N.shape.length != 1 && N.shape[0] != 1 && N.shape[1] != 1 )
            throw new IllegalDimensionException(
                "Input must be a row or column vector (shape = N or 1xN or Nx1)."
            );

        int size = N.dataReal.length;
        Complex result = new Complex( size, size );
        for ( int ind, i = 0; i < size; i++ ) {
            ind = i * size + i;
            result.dataReal[ind] = N.dataReal[i];
            result.dataImag[ind] = N.dataImag[i];
        }
        return result;
    }

    /**************************************************************************
     * <p>Check whether this Complex is a scalar (size of 1)
     *
     * @return true if this Complex is scalar
     *************************************************************************/
    public boolean isScalar ( ) {
        return dataReal.length == 1;
    }

    /**************************************************************************
     * <p>Check whether this Complex is a scalar (size of 1)
     *
     * @return true if this Complex is scalar
     *************************************************************************/
    public boolean isVector ( ) {
        if ( shape.length == 1 )
            return true;
        return shape.length == 2 && (shape[0] == 1 || shape[1] == 1);
    }

    /**************************************************************************
     * <p>If this Complex is a scalar, return the scalar real value
     *
     * @return the real value of the complex scalar
     *************************************************************************/
    public double valueReal ( ) {
        if ( !isScalar() )
            throw new IllegalDimensionException(
                "Error using value(): Complex must be scalar."
            );
        return dataReal[0];
    }

    /**************************************************************************
     * <p>If this Complex is a scalar, return the scalar imaginary value
     *
     * @return the imaginary value of the complex scalar
     *************************************************************************/
    public double valueImag ( ) {
        if ( !isScalar() )
            throw new IllegalDimensionException(
                "Error using value(): Complex must be scalar."
            );
        return dataImag[0];
    }

    /**************************************************************************
     * <p>Gets a reference to the raw real data contained in the Complex. Only
     * recommended if you need fast access to the data in the object.
     *
     * @return a reference to the Logical data.
     *************************************************************************/
    public double[] getDataReal ( ) {
        return dataReal;
    }

    /**************************************************************************
     * <p>Gets a reference to the raw imaginary data contained in the Logical. Only
     * recommended if you need fast access to the data in the object. If the
     * Complex contains only real numbers, the imaginary component is null.
     *
     * @return a reference to the Logical data.
     *************************************************************************/
    public double[] getDataImag ( ) {
        return dataImag;
    }

    /**************************************************************************
     * <p>Sets the value of the data at the subscript with the inputted value
     * WITHOUT resizing if the subscript is out of bounds of the data
     * dimensions. If this is the case, an IllegalDimensionException will be
     * thrown.
     *
     * @param inDataReal    the data to set (real component)
     * @param inDataImag    the data to set (imaginary component)
     * @param sub           the subscript at which to set the data
     *************************************************************************/
    public void set ( double inDataReal, double inDataImag, int[] sub ) {
        int ind = sub2ind( shape, sub );
        if ( ind >= dataReal.length ) {
            throw new IllegalDimensionException(
                "Subscript out of bounds for dimensions"
            );
        }
        dataReal[ind] = inDataReal;
        dataImag[ind] = inDataImag;
    }

    /**************************************************************************
     * <p>Sets the value of the data at the subscript with the inputted value
     * WITHOUT resizing if the subscript is out of bounds of the data
     * dimensions.
     *
     * @param inDataReal    the data to set (real component)
     * @param inDataImag    the data to set (imaginary component)
     * @param ind           the index at which to set the data
     *************************************************************************/
    public void set ( double inDataReal, double inDataImag, int ind ) {
        dataReal[ind] = inDataReal;
        dataImag[ind] = inDataImag;
    }

    /**************************************************************************
     * <p>Sets the value of the data at the subscript with the inputted value. If
     * the subscript is out of bounds of the data dimensions, the data will be
     * resized. This can be slower on larger data sets as this
     * this operation is O(n) for data of size n.
     *
     * @param inDataReal    the data to set (real component)
     * @param inDataImag    the data to set (imaginary component)
     * @param sub           the subscript at which to set the data
     *************************************************************************/
    public void setAt ( double inDataReal, double inDataImag, int[] sub ) {
        int ind = sub2ind( sub );
        if ( ind >= dataReal.length ) {
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
        dataReal[ind] = inDataReal;
        dataImag[ind] = inDataImag;
    }

    /**************************************************************************
     * <p>Indexes the Complex using the inputted Logical as a mask and sets the
     * indexed values to val.
     *
     * @param valReal   The real component of the value to set the indexed
     *                  elements to
     * @param valImag   The imaginary component of the value
     * @param Inds      Logical to be used as an index into the data
     *************************************************************************/
    public void set ( double valReal, double valImag, Logical Inds ) {
        if ( !NDArray.dimensionsMatch( this, Inds ) )
            throw new IllegalDimensionException(
                "Index dimensions must be equal to data dimensions."
            );

        for ( int i = 0; i < dataReal.length; i++ ) {
            if ( Inds.data[i] ) {
                dataReal[i] = valReal;
                dataImag[i] = valImag;
            }
        }
    }

    /**************************************************************************
     * <p>Indexes the Complex using the inputted Logical as a mask and sets the
     * indexed values to val + current value.
     *
     * @param valReal   The real component of the value to set the indexed
     *                  elements to
     * @param valImag   The imaginary component of the value
     * @param Inds      Logical to be used as an index into the data
     *************************************************************************/
    public void setAdd ( double valReal, double valImag, Logical Inds ) {
        if ( !NDArray.dimensionsMatch( this, Inds ) )
            throw new IllegalDimensionException(
                    "Index dimensions must be equal to data dimensions."
            );

        for ( int i = 0; i < dataReal.length; i++ ) {
            if ( Inds.data[i] ) {
                dataReal[i] += valReal;
                dataImag[i] += valImag;
            }
        }
    }

    /**************************************************************************
     * <p>Indexes the Complex using the inputted Logical as a mask and sets the
     * indexed values to val * current value.
     *
     * @param valReal   The real component of the value to set the indexed
     *                  elements to
     * @param valImag   The imaginary component of the value
     * @param Inds      Logical to be used as an index into the data
     *************************************************************************/
    public void setMul ( double valReal, double valImag, Logical Inds ) {
        if ( !NDArray.dimensionsMatch( this, Inds ) )
            throw new IllegalDimensionException(
                    "Index dimensions must be equal to data dimensions."
            );

        double x,y,u,v;
        for ( int i = 0; i < dataReal.length; i++ ) {
            if ( Inds.data[i] ) {
                x = dataReal[i];
                y = dataImag[i];
                u = valReal;
                v = valImag;
                dataReal[i] = x*u - y*v;
                dataImag[i] = x*v + y*u;
            }
        }
    }

    /**************************************************************************
     * <p>Indexes the Complex using the inputted Logical as a mask. Returns a row
     * vector with the indexed values.
     *
     * @param Inds      Logical to be used as an index into the data
     *
     * @return a reference to the initialized vector.
     *************************************************************************/
    public Complex get ( Logical Inds ) {

        if ( !NDArray.dimensionsMatch( this, Inds ) )
            throw new IllegalDimensionException(
                "Index dimensions must be equal to data dimensions."
            );
        // get number of elements
        int numElements = 0;
        for ( int i = 0; i < dataReal.length; i++ )
            if ( Inds.data[i] ) numElements++;

        Complex result = new Complex( 1, numElements );
        for ( int ind = 0, i = 0; i < dataReal.length; i++ ) {
            if ( Inds.data[i] ) {
                result.dataReal[ind] = dataReal[i];
                result.dataImag[ind] = dataImag[i];
                ind++;
            }
        }
        return result;
    }

    /**************************************************************************
     * <p>Returns a scalar Complex with the value at the inputted subscript.
     *
     * @param sub   the subscript to retrieve data at
     *
     * @return a reference to the initialized scalar.
     *************************************************************************/
    public Complex get ( int... sub ) {
        int ind = sub2ind( shape, sub );
        if ( ind >= dataReal.length )
            throw new IllegalDimensionException(
                "Subscript out of range of data dimensions");

        Complex ret = new Complex( 1,1 );
        ret.dataReal[0] = dataReal[ind];
        ret.dataImag[0] = dataImag[ind] ;
        return ret;
    }

    /**************************************************************************
     * <p>Returns a scalar Complex with the value at the inputted linear
     * index.
     *
     * @param ind   the linear index to retrieve data at
     *
     * @return a reference to the initialized scalar.
     *************************************************************************/
    public Complex get ( int ind ) {
        if ( ind >= dataReal.length )
            throw new IllegalDimensionException(
                    "Subscript out of range of data dimensions");

        Complex ret = new Complex( 1,1 );
        ret.dataReal[0] = dataReal[ind];
        ret.dataImag[0] = dataImag[ind];
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
        double max = dataReal[0];
        int ind = 0;
        for ( int i = 0; i < dataReal.length; i++ ) {
            if ( dataReal[i] > max) {
                max = dataReal[i];
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
        double max = dataReal[0];
        double absVal, tmp;
        int ind = 0;
        for ( int i = 0; i < dataReal.length; i++ ) {
            tmp = dataReal[i];
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
        double min = dataReal[0];
        int ind = 0;
        for ( int i = 0; i < dataReal.length; i++ ) {
            if ( dataReal[i] < min) {
                min = dataReal[i];
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
        double min = dataReal[0];
        double absVal, tmp;
        int ind = 0;
        for ( int i = 0; i < dataReal.length; i++ ) {
            tmp = dataReal[i];
            absVal = (tmp <= 0.0) ? 0.0 - tmp : tmp;
            if ( absVal < min) {
                min = absVal;
                ind = i;
            }
        }
        return ind;
    }

    /**************************************************************************
     * <p>Reshapes the Complex by simply changing the dimensions, not the data.
     * This means the output will still be in row-major order and the dimensions
     * must satisfy the requirement that the number of elements must not change.
     *
     * @param dimensions    The dimensions to shape the data to
     *
     * @return a reference to the reshaped Complex
     *************************************************************************/
    public Complex reshape ( int ...dimensions ) {
        int size = 1;
        for ( int n : dimensions )
            size *= n;
        if ( size != dataReal.length )
            throw new IllegalDimensionException(
                "Invalid Dimensions. Number of elements must be the same"
            );
        Complex reshaped = new Complex( dimensions );
        copyData( reshaped );
        return reshaped;
    }

    /**************************************************************************
     * <p>Transposes the Complex. Only defined for Matrices (<3 dimensions); any
     * input of a higher dimension than 2 will cause an exception
     *
     * @return a reference to the transposed data.
     *************************************************************************/
    public Complex transpose ( ) {
        if ( shape.length > 2 )
            throw new IllegalDimensionException(
                    "Transpose operation is not defined for data with more than 2 dimensions"
            );

        Complex transposed;
        if ( shape.length == 1 )
            transposed = new Complex( shape[0], 1 );
        else
            transposed = new Complex( shape[1], shape[0] );

        int r = transposed.shape[0];
        int c = transposed.shape[1];
        for ( int i = 0; i < r; i++ )
            for (int j = 0; j < c; j++) {
                transposed.dataReal[i*c+j] = dataReal[j*r+i];
                transposed.dataImag[i*c+j] = dataImag[j*r+i];
            }

        return transposed;
    }

    /**************************************************************************
     * <p>flattens the Complex to one dimension without changing any of the data.
     *
     * @return a reference to the flattened Complex
     *************************************************************************/
    public Complex flatten ( ) {
        Complex flat = new Complex( 1, dataReal.length );
        copyData( flat );
        return flat;
    }

    /**************************************************************************
     * <p>Creates a deep copy of this Complex object and returns a reference to
     * the copy.
     *
     * @return a reference to the copied object
     *************************************************************************/
    public Complex copy ( ) {
        Complex copy = new Complex( this.shape );
        copyData( copy );
        return copy;
    }

    /**************************************************************************
     * <p>Returns a Complex with all dimensions of length 1 removed. If the
     * resulting data is one dimensional, the resulting object will be a row
     * vector (1xN) if the first dimension (row) was length 1, and a column
     * vector (Nx1) otherwise.
     *
     * @return a reference to the resulting Complex.
     *************************************************************************/
    public Complex squeeze ( ) {
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
        Complex squeezed = new Complex( newShape );
        copyData( squeezed );
        return squeezed;
    }

    /**************************************************************************
     * <p>Resizes the array containing the Complex data while keeping data in its
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

        boolean resizeImag = dataImag != null;
        double[] newDataReal = new double[newSize];
        double[] newDataImag = new double[newSize];

        int[] sub;
        int newInd;
        for ( int i = 0; i < dataReal.length; i++ ) {
            sub = ind2sub( shape, i );
            newInd = sub2ind( dimensions, sub );
            newDataReal[newInd] = dataReal[i];
            newDataImag[newInd] = dataImag[i];
        }
        dataReal = newDataReal;
        dataImag = newDataImag;
        shape = dimensions;
    }

    /**************************************************************************
     * <p>Slices the data into a sub Complex
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
     * @return a reference to the sliced Complex.
     *************************************************************************/
    public Complex slice ( int[] ...dimensions ) {
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

        double[][] newData = sliceData( dimensions, newShape, size );

        Complex sliced = new Complex( newShape );
        sliced.dataReal = newData[0];
        sliced.dataImag = newData[1];
        return sliced;
    }

    /**************************************************************************
     * <p>Concatenate the inputted Complex to the end of this Complex. Note that
     * the number of dimensions and the dimension lengths, except for the
     * dimension being concatenated along, must be equal.
     *
     * @param dimension the dimension to concatenate along
     * @param N         the Complex to concatenate
     *
     * @return a reference to the new Complex
     *************************************************************************/
    public Complex concat ( int dimension, Complex N ) {
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

        Complex catted = this.copy( );
        int[] newShape = catted.shape( );
        newShape[dimension] = shape[dimension] + N.shape[dimension];
        catted.resizeNoCheck( newShape );

        int[] sub;
        int ind;
        for ( int i = 0; i < N.dataReal.length; i++ ) {
            sub = ind2subNoCheck( N.shape, i );
            sub[dimension] += shape[dimension];
            ind = sub2indNoCheck( newShape, sub );
            catted.dataReal[ind] = N.dataReal[i];
            catted.dataImag[ind] = N.dataImag[i];
        }

        return catted;
    }

    /**************************************************************************
     * <p>Creates a string representation of the Complex for printing. Will only
     * show actual data for 1 and 2 dimensional data as higher dimensional
     * data is difficult to display well in a string.
     *
     * @return a string representation of the Logical
     *************************************************************************/
    public String toString ( ) {
        boolean isComplex = dataImag != null;
        int maxInd = findAbsMax( );
        double max = dataReal[maxInd];
        int L = String.format( "%.4f", max ).length();

        StringBuilder s = new StringBuilder(super.toString() + " <Complex>\n");
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
                if ( dataReal[i * col] < 0 )
                    if ( isComplex )
                        s.append("[  ");
                    else
                        s.append("[ ");
                else
                if ( isComplex )
                    s.append("[   ");
                else
                    s.append("[  ");
                for ( int j = 0; j < col; j++ ) {
                    String tmp;
                    double real = dataReal[i * col + j];
                    if ( real < 0 )
                        tmp = String.format("%" + (L+1) + ".4f", dataReal[i * col + j]);
                    else
                        tmp = String.format("%" + L + ".4f", dataReal[i * col + j]);
                    if ( isComplex ) {
                        double img = dataImag[i * col + j];
                        if ( img < 0 )
                            tmp = tmp + " -" + String.format("%" + L + ".4fi", -img);
                        else
                            tmp = tmp + " +" + String.format("%" + L + ".4fi", img);
                    }
                    if ( j < col-1 && dataReal[i * col + j+1] < 0 )
                        if ( isComplex )
                            s.append(tmp).append("  ");
                        else
                            s.append(tmp).append(" ");
                    else
                    if ( isComplex )
                        s.append(tmp).append("   ");
                    else
                        s.append(tmp).append("  ");
                }
                s.append("]\n");
            }
        }

        return s.toString();
    }

    /**************************************************************************
     * <p>Compares two Complex objects to check if they are equal in both
     * dimension and data. NOT RECOMMENDED when dealing with arithmetic because
     * of the possibility of precision issues. Use equalTolerance() instead.
     *
     * @param N the other Complex to compare this one to
     *
     * @return  true if the two Complex objects are equal, false otherwise.
     *************************************************************************/
    public boolean equals ( Complex N ) {
        if ( this == N )
            return true;
        if ( !NDArray.dimensionsMatch( this, N ) )
            return false;

        for ( int i = 0; i < dataReal.length; i++ ) {
            if ( this.dataReal[i] != N.dataReal[i] )
                return false;
            if ( this.dataImag[i] != N.dataImag[i] )
                return false;
        }

        return true;
    }

    /**************************************************************************
     * <p>Compares two Complex objects to check if they are equal with a
     * tolerance
     *
     * @param N         the other Complex to compare to
     * @param tolerance the tolerance for the difference between two
     *                  corresponding data in the two Complexes
     *
     * @return true if the two Complexes are equal to the specified tolerance
     *************************************************************************/
    public boolean equalsTolerance ( Complex N, double tolerance ) {
        if ( this == N )
            return true;
        if ( !NDArray.dimensionsMatch( this, N ) )
            return false;


        double diff;
        for ( int i = 0; i < dataReal.length; i++ ) {
            diff = this.dataReal[i] - N.dataReal[i];
            diff = (diff <= 0.0) ? 0.0 - diff : diff; // get abs value
            if ( diff > tolerance )
                return false;
            diff = Math.abs(this.dataImag[i] - N.dataImag[i]);
            diff = (diff <= 0.0) ? 0.0 - diff : diff;
            if (diff > tolerance)
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

        boolean resizeImag = dataImag != null;
        double[] newDataReal = new double[newSize];
        double[] newDataImag = new double[newSize];

        int[] sub;
        int newInd;
        for ( int i = 0; i < dataReal.length; i++ ) {
            sub = ind2sub( shape, i );
            newInd = sub2ind( dimensions, sub );
            newDataReal[newInd] = dataReal[i];
            newDataImag[newInd] = dataImag[i];
        }
        dataReal = newDataReal;
        dataImag = newDataImag;
        shape = dimensions;
    }

    /* Moved this to its own function bc slice was getting a little long */
    protected double[][] sliceData ( int[][] dims, int[] newShape, int size ) {
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
        double[] newDataReal = new double[size];
        double[] newDataImag = new double[size];
        int startInd = sub2indNoCheck( shape, sub );
        int endInd = sub2indNoCheck( shape, subLast );

        ind = 0;
        newDataReal[ind] = dataReal[startInd];
        newDataImag[ind] = dataImag[startInd];
        ind++;

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
            if ( add ) {
                newDataReal[ind] = dataReal[i];
                newDataImag[ind] = dataImag[i];
                ind++;
            }
        }

        if ( startInd != endInd ) {
            newDataReal[ind] = dataReal[endInd];
            newDataImag[ind] = dataImag[endInd];
        }

        return new double[][]{ newDataReal, newDataImag };
    }

    /* copy this objects data into dest */
    protected void copyData ( Complex dest ) {
        System.arraycopy( dataReal, 0, dest.dataReal, 0, dataReal.length );
        System.arraycopy( dataImag, 0, dest.dataImag, 0, dataImag.length );
    }

}

