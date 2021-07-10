package com.evanstella.jmp.core;

public class Mat2 {

    protected double[] data;
    protected int numRow;
    protected int numCol;


    public Mat2 ( int r, int c ) {
        if ( r < 1 || c < 1 )
            throw new IllegalMatrixDimensionException( "Matrix dimensions must be positive" );
        numRow = r;
        numCol = c;
        data = new double[numRow*numCol];
    }

    public void fill ( double start, double step ) {}

}
