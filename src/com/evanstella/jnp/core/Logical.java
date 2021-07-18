package com.evanstella.jnp.core;


public class Logical extends NDArray {

    protected boolean[] data;


    /**************************************************************************
     * Class constructor. Initializes a Logical with parameterized dimensions.
     * The dimensions are deep copied to limit their write access to the object.
     *
     * @param dimensions The arbitrary dimensions for the N-D Logical.
     *************************************************************************/
    public Logical ( int ...dimensions ) {
        shape = new int[dimensions.length];
        System.arraycopy(dimensions, 0, shape,0, dimensions.length);
        int size = 1;
        for ( int n : dimensions )
            size *= n;
        data = new boolean[size];
    }

    /**************************************************************************
     * Class constructor. Initialize a Logical from a boolean[]. Only a
     * shallow copy of the array is made
     *
     * @param indata The array to initialize from
     *************************************************************************/
    public Logical ( boolean[] indata ) {
        data = indata;
        shape = new int[]{ indata.length };
    }

    /**************************************************************************
     * Class constructor. Initialize a Logical from a boolean[][]. Makes a
     * a deep copy in order to aggregate the boolean[][] into a boolean[].
     *
     * @param indata The array to initialize from
     *************************************************************************/
    public Logical ( boolean[][] indata ) {
        data = new boolean[indata.length * indata[0].length];
        int ind = 0;
        for ( boolean[] r : indata ) {
            if ( r.length != indata[0].length )
                throw new IllegalDimensionException(
                    "Input dimensions must be consistent."
                );
            for ( boolean b : r )
                data[ind++] = b;
        }
        shape = new int[]{ indata.length, indata[0].length };
    }

    /**************************************************************************
     * Initializes a Logical with all true values.
     *
     * @param dimensions The dimensions for the N-D Logical.
     *
     * @return a reference the initialized Logical
     *************************************************************************/
    public static Logical True ( int ...dimensions ) {
        Logical L = new Logical( dimensions );
        java.util.Arrays.fill(L.data, true);
        return L;
    }

    /**************************************************************************
     * Initializes a Logical with all false values. Literally the same as
     * calling the constructor but here for completeness' sake.
     *
     * @param dimensions the dimensions for the N-D Logical.
     *
     * @return a reference the initialized Logical
     *************************************************************************/
    public static Logical False ( int ...dimensions ) {
        return new Logical( dimensions );
    }

    /**************************************************************************
     * Initializes a Logical with random true/false values using
     * java.util.Random.
     *
     * @param dimensions the dimensions for the N-D Logical.
     *
     * @return a reference the initialized Logical
     *************************************************************************/
    public static Logical Rand ( int ...dimensions ) {
        java.util.Random R = new java.util.Random( );
        Logical L = new Logical( dimensions );

        for ( int i = 0; i < L.data.length; i++ )
            L.data[i] = R.nextBoolean();

        return L;
    }

    /**************************************************************************
     * Initializes a Logical with random true/false values using with a seed
     * java.util.Random.
     *
     * @param seed       the PRNG seed.
     * @param dimensions the dimensions for the N-D Logical.
     *
     * @return a reference the initialized Logical
     *************************************************************************/
    public static Logical Rand ( long seed,  int ...dimensions ) {
        java.util.Random R = new java.util.Random( seed );
        Logical L = new Logical( dimensions );

        for ( int i = 0; i < L.data.length; i++ )
            L.data[i] = R.nextBoolean();

        return L;
    }

    /**************************************************************************
     * Gets a reference to the raw data contained in the Logical. Only
     * recommended if you need fast access to the data in the object.
     *
     * @return a reference to the Logical data.
     *************************************************************************/
    public boolean[] getData ( ) {
        return data;
    }

    /**************************************************************************
     * Sets the value of the data at the subscript with the inputted value
     * WITHOUT resizing if the subscript is out of bounds of the data
     * dimensions. If this is the case, an IllegalDimensionException will be
     * thrown.
     *
     * @param indata    the data to set
     * @param sub       the subscript at which to set the data
     *************************************************************************/
    public void set ( boolean indata, int ...sub ) {
        int ind = sub2ind( sub );
        if ( ind >= data.length ) {
            throw new IllegalDimensionException(
                "Subscript out of bounds for dimensions"
            );
        }
        data[ind] = indata;
    }

    /**************************************************************************
     * Sets the value of the data at the subscript with the inputted value
     * WITHOUT resizing if the subscript is out of bounds of the data
     * dimensions.
     *
     * @param indata    the data to set
     * @param ind       the index at which to set the data
     *************************************************************************/
    public void set ( boolean indata, int ind ) {
        data[ind] = indata;
    }

    /**************************************************************************
     * Sets the value of the data at the subscript with the inputted value. If
     * the subscript is out of bounds of the data dimensions, the data will be
     * resized. This can be slower on larger data sets as this
     * this operation is O(n) for data of size n.
     *
     * @param indata    the data to set
     * @param sub       the subscript at which to set the data
     *************************************************************************/
    public void setAt ( boolean indata, int ...sub ) {
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
     * Returns the data at the inputted subscript.
     *
     * @param sub   the subscript to retrieve data at
     *
     * @return the value at the subscript
     *************************************************************************/
    public boolean get ( int... sub ) {
        int ind = sub2ind( shape, sub );
        if ( ind >= data.length )
            throw new IllegalDimensionException(
                "Subscript out of range of data dimensions");
        return data[ind];
    }

    /**************************************************************************
     * Returns the data at the inputted linear index.
     *
     * @param ind   the linear index to retrieve data at
     *
     * @return the value at the linear index
     *************************************************************************/
    public boolean get ( int ind ) {
        if ( ind >= data.length )
            throw new IllegalDimensionException(
                "Index out of range of data size");
        return data[ind];
    }

    /**************************************************************************
     * Reshapes the Logical by simply changing the dimensions, not the data.
     * This means the output will still be in row-major order and the dimensions
     * must satisfy the requirement that the number of elements must not change.
     *
     * @param dimensions    The dimensions to shape the data to
     *
     * @return a reference to the reshaped Logical
     *************************************************************************/
    public Logical reshape ( int ...dimensions ) {
        int size = 1;
        for ( int n : dimensions )
            size *= n;
        if ( size != data.length )
            throw new IllegalDimensionException(
                "Invalid Dimensions. Number of elements must be the same"
            );
        Logical reshaped = new Logical( dimensions );
        System.arraycopy( data, 0, reshaped.data, 0, data.length );
        return reshaped;
    }

    /**************************************************************************
     * Transposes the Logical. Only defined for Matrices (<3 dimensions); any
     * input of a higher dimension than 2 will cause an exception
     *
     * @return a reference to the transposed data.
     *************************************************************************/
    public Logical transpose ( ) {
        if ( shape.length > 2 )
            throw new IllegalDimensionException(
                "Transpose operation is not defined for data with more than 2 dimensions"
            );

        Logical transposed;
        if ( shape.length == 1 )
            transposed = new Logical( shape[0], 1 );
        else
            transposed = new Logical( shape[1], shape[0] );

        int r = transposed.shape[0];
        int c = transposed.shape[1];
        for ( int i = 0; i < r; i++ )
            for (int j = 0; j < c; j++)
                transposed.data[i*c+j] = data[j*r+i];

        return transposed;
    }

    /**************************************************************************
     * flattens the Logical to one dimension without changing any of the data.
     *
     * @return a reference to the flattened Logical
     *************************************************************************/
    public Logical flatten ( ) {
        Logical flat = new Logical( 1, data.length );
        System.arraycopy( data, 0, flat.data, 0, data.length );
        return flat;
    }

    /**************************************************************************
     * Creates a deep copy of this Logical object and returns a reference to
     * the copy.
     *
     * @return a reference to the copied object
     *************************************************************************/
    public Logical copy ( ) {
        Logical copy = new Logical( this.shape );
        System.arraycopy( data, 0, copy.data, 0, data.length );
        return copy;
    }

    /**************************************************************************
     * Returns a Logical with all dimensions of length 1 removed. If the
     * resulting data is one dimensional, the resulting object will be a row
     * vector (1xN) if the first dimension (row) was length 1, and a column
     * vector (Nx1) otherwise.
     *
     * @return a reference to the squeezed Logical.
     *************************************************************************/
    public Logical squeeze ( ) {
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
        Logical squeezed = new Logical( newShape );
        System.arraycopy( data, 0, squeezed.data, 0, data.length );
        return squeezed;
    }

    /**************************************************************************
     * Resizes the array containing the Logical data while keeping data in its
     * current subscribed position.
     *
     * Should be done infrequently as this requires copying each element to a
     * new array. There is also a little overhead associated with calculating
     * the new element indices. Elements in newly allocated space are
     * initialized to 0 (false)
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
        boolean[] newData = new boolean[newSize];

        int[] sub;
        int newInd;
        for ( int i = 0; i < data.length; i++ ) {
            sub = ind2sub( shape, i );
            newInd = sub2ind( dimensions, sub );
            newData[newInd] = data[i];
        }

        data = newData;
        shape = dimensions;
    }

    /**************************************************************************
     * Slices the data into a sub Logical
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
     * @return a reference to the sliced Logical.
     *************************************************************************/
    public Logical slice ( int[] ...dimensions ) {
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

        boolean[] newData = sliceData( dimensions, newShape, size );

        Logical sliced = new Logical( newShape );
        sliced.data = newData;
        return sliced;
    }

    /**************************************************************************
     * Concatenate the inputted Logical to the end of this Logical. Note that
     * the number of dimensions and the dimensions lengths, except for the
     * dimension being concatenated along, must be equal.
     *
     * @param dimension the dimension to concatenate along
     * @param L         the Logical to concatenate
     *
     * @return a reference to the new Logical
     *************************************************************************/
    public Logical concat ( int dimension, Logical L ) {
        if ( dimension >= shape.length )
            throw new IllegalDimensionException(
                "Logical does not have " + (dimension+1) + " dimensions.");
        if ( shape.length != L.shape.length )
            throw new IllegalDimensionException(
                "Both objects must have the same number of dimensions.");
        for ( int i = 0; i < shape.length; i++ ) {
            if ( i != dimension && shape[i] != L.shape[i] )
                throw new IllegalDimensionException(
                "All dimensions but the concat dimension must be equal in length.");
        }

        Logical catted = this.copy( );
        int[] newShape = catted.shape( );
        newShape[dimension] = shape[dimension] + L.shape[dimension];
        catted.resizeNoCheck( newShape );
        int[] sub;
        int ind;
        for ( int i = 0; i < L.data.length; i++ ) {
            sub = ind2subNoCheck( L.shape, i );
            sub[dimension] += shape[dimension];
            ind = sub2indNoCheck( newShape, sub );
            catted.data[ind] = L.data[i];
        }

        return catted;
    }

    /**************************************************************************
     * Creates a string representation of the Logical for printing. Will only
     * show actual data for 1 and 2 dimensional data as higher dimensional
     *
     *
     * data is difficult to display well in a string.
     *
     * @return a string representation of the Logical
     *************************************************************************/
    public String toString ( ) {
        StringBuilder s = new StringBuilder(super.toString() + " <Logical>\n");
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
                s.append("[ ");
                for ( int j = 0; j < col; j++ ) {
                    String tmp = data[i*col+j] ? "1 ": "0 ";
                    s.append(tmp);
                }
                s.append("]\n");
            }
        }

        return s.toString();
    }

    /**************************************************************************
     * Compares two Logical objects to check if they are equal in both
     * dimension and data.
     *
     * @param L the other Logical to compare this one to
     *
     * @return  true if the two Logical objects are equal, false otherwise.
     *************************************************************************/
    public boolean equals ( Logical L ) {
        if ( this == L )
            return true;
        if ( L.shape.length != this.shape.length )
            return false;
        for ( int i = 0; i < shape.length; i++ ) {
            if ( this.shape[i] != L.shape[i] )
                return false;
        }
        for ( int i = 0; i < data.length; i++ ) {
            if ( this.data[i] != L.data[i] )
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
        boolean[] newData = new boolean[newSize];

        int[] sub;
        int newInd;
        for ( int i = 0; i < data.length; i++ ) {
            sub = ind2sub( shape, i );
            newInd = sub2ind( dimensions, sub );
            newData[newInd] = data[i];
        }

        data = newData;
        shape = dimensions;
    }

    /* Moved this to its own function bc slice was getting a little long */
    protected boolean[] sliceData ( int[][] dims, int[] newShape, int size ) {
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
        boolean[] newData = new boolean[size];
        int startInd = sub2indNoCheck( shape, sub );
        int endInd = sub2indNoCheck( shape, subLast );
        ind = 0;
        newData[ind++] = data[startInd];
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
                newData[ind++] = data[i];
        }
        if ( startInd != endInd )
            newData[ind] = data[endInd];
        return newData;
    }

}
