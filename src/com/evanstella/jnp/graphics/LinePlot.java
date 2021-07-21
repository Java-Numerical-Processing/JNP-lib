package com.evanstella.jnp.graphics;

import com.evanstella.jnp.core.IllegalDimensionException;
import com.evanstella.jnp.core.Numeric;

import javax.swing.JComponent;
import java.awt.*;
import java.awt.event.*;

public class LinePlot extends JComponent {

    private Numeric XData;
    private Numeric YData;

    private double scaleFactor;
    private double scale;
    private int originX;
    private int originY;

    public LinePlot ( Numeric XData, Numeric YData, int h, int w ) {
        super();
        this.XData = XData;
        this.YData = YData;
        originX = h/2;
        originY = w/2;
        scale = 10;
        scaleFactor = 1.1;

        MouseHandler mouseHandler = new MouseHandler( this );
    }

    public void paintComponent (Graphics g) {
        update( g );
    }


    public void update ( Graphics g ) {
        //w is x, and h is y (as in x/y values in a graph)
        int w = originX;
        int h = originY;

        Graphics2D g1 = (Graphics2D) g;
        g1.setRenderingHint(
                RenderingHints.KEY_ANTIALIASING,
                RenderingHints.VALUE_ANTIALIAS_ON);
        g1.setStroke(new BasicStroke(2));
        g1.setColor(Color.black);
        g1.drawLine(0,h,getWidth(),h);
        g1.drawLine(w,0,w,getHeight());
        g1.drawString("0", w - 7, h + 13);

        Graphics2D g2 = (Graphics2D) g;
        g2.setRenderingHint(
                RenderingHints.KEY_ANTIALIASING,
                RenderingHints.VALUE_ANTIALIAS_ON);
        g2.setStroke(new BasicStroke((float)2));
        g2.setColor(Color.red);
        //line1
        Polygon p = new Polygon();
        double[] XPoints = XData.getDataReal();
        double[] YPoints = YData.getDataReal();

        for (int x = 0; x < XPoints.length; x++) {
            p.addPoint((int) (w+scale*XPoints[x]), (int) (h - scale* YPoints[x]));
        }
        g2.drawPolyline(p.xpoints, p.ypoints, p.npoints);
    }

    public void setXData ( Numeric XData ) {
        if ( !XData.isVector() && !XData.isScalar() ) {
            throw new IllegalDimensionException(
                "Line Plot: XData must be a vector or scalar quantity."
            );
        }
        this.XData = XData;
    }

    public Numeric getXData ( ) {
        return XData;
    }

    public void setYData ( Numeric YData ) {
        if ( !YData.isVector() && !YData.isScalar() ) {
            throw new IllegalDimensionException(
                    "Line Plot: YData must be a vector or scalar quantity."
            );
        }
        this.YData = YData;
    }

    public Numeric getYData ( ) {
        return YData;
    }


    private class MouseHandler implements MouseMotionListener, MouseWheelListener {

        public MouseHandler ( JComponent plot ) {
            plot.addMouseMotionListener(this);
            plot.addMouseWheelListener(this);
        }

        @Override
        public void mouseDragged(MouseEvent e) {
            originX = e.getX();
            originY = e.getY();
            repaint();
        }

        @Override
        public void mouseMoved(MouseEvent e) {}

        @Override
        public void mouseWheelMoved(MouseWheelEvent e) {
            int notches = e.getWheelRotation();
            if ( notches < 0 ) scale /= scaleFactor;
            else scale *= scaleFactor;

            if (scale < 1) scale = 1;
            if (scale > 10000) scale = 10000;
            repaint();
        }
    }

}