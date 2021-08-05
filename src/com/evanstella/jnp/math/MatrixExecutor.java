package com.evanstella.jnp.math;

import com.evanstella.jnp.core.Complex;
import com.evanstella.jnp.core.IllegalDimensionException;
import com.evanstella.jnp.core.Numeric;

import java.util.concurrent.CountDownLatch;

public final class MatrixExecutor extends ParallelExecutor {

    public MatrixExecutor ( int threadCount ) {
        super(threadCount);
    }


    public Numeric transpose ( Numeric matrix ) {
        Matrix.validateMatrixFatal( matrix );

        int[] shape = matrix.shape();
        double[] data = matrix.getData();
        final int c = shape[1], r = shape[0];

        Numeric transposed = new Numeric( shape[1], shape[0] );
        double[] newData = transposed.getData();

        //create worker
        CountDownLatch count = new CountDownLatch(threadCount);
        class worker implements Runnable {
            final int startIdx, endIdx;
            public worker ( int start, int end ) { startIdx = start; endIdx = end; }
            public void run ( ) {
                try {
                    int curR, curC;
                    for (int i = startIdx; i < endIdx; i++) {
                        curR = i / c;
                        curC = i % c;
                        newData[ curC*r+curR ] = data[i];
                    }
                }
                finally { count.countDown(); }
            }
        }

        int start = 0, increment = data.length / threadCount;
        for ( int i = 0; i < threadCount-1; i++ ) {
            executorService.execute( new worker( start, start+increment ) );
            start = start+increment;
        }
        executorService.execute( new worker( start, data.length ) );
        this.await( count );
        return transposed;
    }

    public Complex ctranspose (Complex matrix ) {
        Matrix.validateMatrixFatal( matrix );

        int[] shape = matrix.shape();
        double[] dataR = matrix.getDataReal();
        double[] dataI = matrix.getDataImag();
        final int c = shape[1], r = shape[0];

        Complex transposed = new Complex( shape[1], shape[0] );
        double[] newDataR = transposed.getDataReal();
        double[] newDataI = transposed.getDataImag();

        //create worker
        CountDownLatch count = new CountDownLatch( threadCount );
        class worker implements Runnable {
            final int startIdx, endIdx;
            public worker ( int start, int end ) { startIdx = start; endIdx = end; }
            public void run ( ) {
                try {
                    int curR, curC;
                    for (int i = startIdx; i < endIdx; i++) {
                        curR = i / c;
                        curC = i % c;
                        newDataR[ curC*r+curR ] =  dataR[i];
                        newDataI[ curC*r+curR ] = -dataI[i];
                    }
                }
                finally { count.countDown(); }
            }
        }

        int start = 0, increment = dataR.length / threadCount;
        for ( int i = 0; i < threadCount-1; i++ ) {
            executorService.execute( new worker( start, start+increment ) );
            start = start+increment;
        }
        executorService.execute( new worker( start, dataR.length ) );
        this.await( count );
        return transposed;
    }


    public Numeric mul ( Numeric m1, Numeric m2 ) {
        Matrix.validateMatrixFatal( m1 );
        Matrix.validateMatrixFatal( m2 );
        if ( m1.shape()[1] != m2.shape()[0] )
            throw new IllegalDimensionException(
                "Matrix multiplication: invalid dimensions."
            );

        final int c1 = m1.shape()[1], c2 = m2.shape()[1];
        final int r = m1.shape()[0], c = m2.shape()[1];
        double[] X = m1.getData(); double[] Y = m2.getData();
        Numeric result = new Numeric( r, c );
        double[] Z = result.getData();

        //create worker
        CountDownLatch count = new CountDownLatch( threadCount );
        class worker implements Runnable {
            final int i,j;
            public worker ( int i, int j ) { this.i = i; this.j = j; }
            public void run ( ) {
                try {
                    int ind = i*c+j;
                    Z[ind] = 0;
                    for( int k = 0; k < c1; k++ )
                        Z[ind] += X[i*c1+k] * Y[k*c2+j];
                }
                finally { count.countDown(); }
            }
        }

        for ( int i = 0; i < r; i++ )
            for( int j = 0; j < c; j++ )
                executorService.execute( new worker(i, j) );

        this.await( count );
        return result;
    }




}
