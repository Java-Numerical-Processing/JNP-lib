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

/******************************************************************************
 * <p>An interface for numerical functions of 2 variables ( z = f(x,y) )
 *
 * @author Evan Stella
 *****************************************************************************/
public interface ComplexExpression2D {

    /**************************************************************************
     * <p>This method needs to be implemented with the logic of the intended
     * single variable function to return the real component of the result.
     *
     * @param reX    the real component of the first input
     * @param imX    the imaginary component of the first input
     * @param reY    the real component of the second input
     * @param imY    the imaginary component of the second input
     *
     * @return the real component of the result of the function
     **************************************************************************/
    double evaluateReal ( double reX, double imX, double reY, double imY );

    /**************************************************************************
     * <p>This method needs to be implemented with the logic of the intended
     * single variable function to return the imaginary component of the result.
     *
     * @param reX    the real component of the first input
     * @param imX    the imaginary component of the first input
     * @param reY    the real component of the second input
     * @param imY    the imaginary component of the second input
     *
     * @return the imaginary component of the result of the function
     **************************************************************************/
    double evaluateImag ( double reX, double imX, double reY, double imY );
}