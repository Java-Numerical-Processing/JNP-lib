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

import com.evanstella.jnp.math.ExecutionInterruptedException;

import java.util.concurrent.CountDownLatch;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

/******************************************************************************
 * <p>ParallelExecutor is a parallel process handler. Other executors can be
 * can be created from an instance of this class so that they use the same
 * executorService.
 *
 * @author Evan Stella
 *****************************************************************************/
public class ParallelExecutor {

    protected int threadCount;
    protected ExecutorService executorService;

    protected ParallelExecutor ( ) { };

    /**************************************************************************
     * <p>Constructor. Start a thread pool with the inputted thread count.
     *
     * @param threadCount   The number of threads.
     *************************************************************************/
    public ParallelExecutor ( int threadCount ) {
        this.threadCount = threadCount;
        executorService = Executors.newFixedThreadPool(threadCount);
    }

    /**************************************************************************
     * <p>Bind another Executor to this object so they use the same executor
     * service.
     *
     * @param E     The other executor to bind to; must inherit from this
     *              class
     *************************************************************************/
    public <Executor extends ParallelExecutor> void bind ( Executor E ) {
        E.executorService = executorService;
        E.threadCount = threadCount;
    }

    /**************************************************************************
     * <p>Shutdown the thread pool.
     *************************************************************************/
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
