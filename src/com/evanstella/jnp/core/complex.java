package com.evanstella.jnp.core;

public class complex {

    public double re;
    public double im;

    public complex ( double real, double im ) {
        this.re = real;
        this.im = im;
    }

    public String toString ( ) {
        String s = re + " ";
        s += ( im < 0 ) ? "-" : "+";
        s += " " + Math.abs(im);
        return s;
    }

}
