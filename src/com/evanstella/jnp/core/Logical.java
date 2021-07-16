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
     * TODO
     *************************************************************************/
    public Logical ( boolean[] indata ) {}

    /**************************************************************************
     * TODO
     *************************************************************************/
    public Logical ( boolean[][] indata ) {}

    /**************************************************************************
     * Initializes a Logical with all true values.
     *
     * @param dimensions The arbitrary dimensions for the N-D Logical.
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
     * @param dimensions the arbitrary dimensions for the N-D Logical.
     *
     * @return a reference the initialized Logical
     *************************************************************************/
    public static Logical False ( int ...dimensions ) {
        return new Logical( dimensions );
    }

    /**************************************************************************
     * Gets a reference to the raw data contained in the logical. Only
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
    public void setFast ( boolean indata, int ...sub ) {
        int ind = sub2ind( sub );
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
    public void setFast ( boolean indata, int ind ) {
        data[ind] = indata;
    }

    /**************************************************************************
     * TODO
     * Sets the value of the data at the subscript with the inputted value. If
     * the subscript is out of bounds of the data dimensions, the data will be
     * resized, this can be undesirable on larger data sets as this
     * this operation is O(n) for data of size n.
     *
     * @param indata    the data to set
     * @param sub       the subscript at which to set the data
     *************************************************************************/
    public void set ( boolean indata, int ...sub ) {
        int ind;
        try {
            ind = sub2ind( sub );
        } catch( IllegalDimensionException E ) {
            if ( sub.length != shape.length )
                throw new IllegalDimensionException(
                        "Number of subscript dimensions must match data dimensions."
                );
            /*TODO*/
            throw new IllegalDimensionException("TODO");
        }
        data[ind] = indata;
    }

    /**************************************************************************
     * TODO
     * Sets the value of the data at the index with the inputted value. If
     * the index is out of bounds of the data dimensions, the data will be
     * resized; this can be undesirable on larger data sets as this
     * this operation is O(n) for data of size n.
     *
     * @param indata    the data to set
     * @param ind       the index at which to set the data
     *************************************************************************/
    public void set ( boolean indata, int ind ) {
        try {
            data[ind] = indata;
        } catch( ArrayIndexOutOfBoundsException E ) {
            /*TODO*/
            throw new IllegalDimensionException("TODO");
        }
    }

    /**************************************************************************
     * TODO
     *************************************************************************/
    public Logical get ( ) {
        return null;
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
     * input of higher dimension than 2 will cause an exception
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
     * flattens the Logical to one dimension, not changing any of the data.
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
     * TODO
     *
     * @param dimension the dimension to concatenate along
     * @param L         the Logical to concatenate
     *
     * @return a reference to the new Logical
     *************************************************************************/
    public Logical concat ( int dimension, Logical L ) {return null;}

    /**************************************************************************
     * Creates a string representation of the Logical for printing. Will only
     * show actual data for 1 and 2 dimensional data as higher dimensional
     * data is difficult to display well in a string.
     *
     * @return a string representation of the Logical
     *************************************************************************/
    public String toString ( ) {
        StringBuilder s = new StringBuilder(super.toString() + " <Logical>\n");
        if ( shape.length < 3 ) {
            int row = shape[0];
            int col = 1;
            if ( shape.length == 2 )
                col = shape[1];
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
     * Private Methods
     *************************************************************************/

    // TODO
    private void resize ( int ...newDims ) {
        int newSize = 1;
        for ( int n : newDims )
            newSize *= n;
        boolean[] newData = new boolean[newSize];
        /*TODO*/
        data = newData;
    }
}
