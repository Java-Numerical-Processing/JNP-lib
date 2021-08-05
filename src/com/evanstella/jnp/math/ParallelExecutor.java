package com.evanstella.jnp.math;

import java.util.concurrent.CountDownLatch;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

public abstract class ParallelExecutor {

    public int threadCount;
    protected final ExecutorService executorService;

    public ParallelExecutor( int threadCount ) {
        this.threadCount = threadCount;
        executorService = Executors.newFixedThreadPool(threadCount);
    }

    public void shutdown ( ) {
        executorService.shutdown();
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
