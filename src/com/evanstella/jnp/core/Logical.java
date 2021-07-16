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
     * Sets the value of the data at the subscript with the inputted value If
     * the subscript is out of bounds of the data dimensions, the data will be
     * resized, but this can be undesirable on larger objects as this
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

    public void set ( boolean indata, int ind ) {
        /*TODO*/
    }

    public Logical get() {
        return null;
    }

    public Logical reshape() {
        return null;
    }

    public Logical diagonal() {
        return null;
    }

    public Logical transpose() {
        return null;
    }

    public Logical flatten() {
        return null;
    }

    public Logical copy() {
        return null;
    }

    /**************************************************************************
     * Gets a reference to the raw data contained in the logical. Only
     * recommended if you need fast access to the data in the object.
     *
     * @return a reference to the Logical data.
     *************************************************************************/
    public boolean[] getDataRaw ( ) {
        return data;
    }

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

    private void resize ( int ...newDims ) {
        int newSize = 1;
        for ( int n : newDims )
            newSize *= n;
        boolean[] newData = new boolean[newSize];
        /*TODO*/
        data = newData;
    }
}
