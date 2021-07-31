package com.evanstella.jnp.math;

import java.util.concurrent.CountDownLatch;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

public class ParallelExecutor {

    public int threadCount;
    protected final ExecutorService executorService;

    // no instances for you
    public ParallelExecutor( int threadCount ) {
        this.threadCount = threadCount;
        executorService = Executors.newFixedThreadPool(threadCount);
    }

    public void shutdown ( ) {
        executorService.shutdown();
    }

    //TODO
    protected static class executionWorker implements Runnable {
        final int startIdx, endIdx;
        public executionWorker ( int start, int end ) {
            startIdx = start;
            endIdx = end;
        }
        public void run ( ) {}
    }

    protected void await ( CountDownLatch count ) {
        try {
            count.await();
        } catch ( InterruptedException e ) {
            throw new ExecutionInterruptedException(
                "Parallel execution was interrupted"
            );
        }
    }

}
