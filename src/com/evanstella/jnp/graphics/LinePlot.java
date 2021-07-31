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

package com.evanstella.jnp.graphics;

import com.evanstella.jnp.core.IllegalDimensionException;
import com.evanstella.jnp.core.Numeric;

import javax.swing.JComponent;
import java.awt.*;
import java.awt.event.*;

/******************************************************************************
 * A JComponent for a line plot graph. TODO.
 *
 * @author Evan Stella
 *****************************************************************************/
public class LinePlot extends JComponent {

    private Numeric XData;
    private Numeric YData;

    private double zoom;
    private double scale;
    private int originX;
    private int originY;
    private int dragAnchorX;
    private int dragAnchorY;

    public LinePlot ( Numeric XData, Numeric YData, int w, int h ) {
        super();
        this.XData = XData;
        this.YData = YData;
        originX = 0;
        originY = 0;
        scale = 10;
        zoom = 1;

        MouseHandler mouseHandler = new MouseHandler( this );
    }

    public void paintComponent (Graphics g) {
        update( g );
    }


    public void update ( Graphics g ) {
        //w is x, and h is y (as in x/y values in a graph)
        int w = 400;
        int h = 300;

        Graphics2D canvas = (Graphics2D) g;
        setTranslation( canvas );
        setZoom( canvas );

        Graphics2D g2d = (Graphics2D) g;
        g2d.setRenderingHint(
            RenderingHints.KEY_ANTIALIASING,
            RenderingHints.VALUE_ANTIALIAS_ON);
        g2d.setStroke(new BasicStroke((float)2));
        g2d.setColor(Color.blue);

        int maxX = 0, maxY = 0, minX = 0, minY = 0;

        Polygon p = new Polygon();
        double[] XPoints = XData.getData();
        double[] YPoints = YData.getData();

        for ( int px, py, x = 0; x < XPoints.length; x++ ) {
            px = (int) (w+scale * XPoints[x]);
            if ( px > maxX ) maxX = px;
            if ( px < minX ) minX = px;

            py = (int) (h-scale * YPoints[x]);
            if ( py > maxY ) maxY = py;
            if ( py < minY ) minY = py;

            p.addPoint( px, py );
        }

        g2d.drawPolyline( p.xpoints, p.ypoints, p.npoints );

        g2d.setStroke(new BasicStroke(2));
        g2d.setColor(Color.black);
        g2d.drawLine(minX-50,h,maxX+50,h);
        g2d.drawLine(w,minY-50,w,maxY+50);
        g2d.drawString("(0,0)", w + 2, h + 15);

    }

    private void setTranslation ( Graphics2D canvas ) {
        canvas.translate(originX,originY);
    }

    private void setZoom ( Graphics2D canvas ) {
        //canvas.translate(originX,originY);
        canvas.scale(zoom,zoom);
        canvas.translate(-originX,-originY);
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
            dragAnchorX = originX;
            dragAnchorY = originY;

            plot.addMouseMotionListener(this);
            plot.addMouseWheelListener(this);
        }

        @Override
        public void mouseDragged ( MouseEvent e ) {
            int deltaX = dragAnchorX - e.getX();
            int deltaY = dragAnchorY - e.getY();

            originX -= deltaX * Math.abs(zoom);
            originY -= deltaY * Math.abs(zoom);

            repaint();

            dragAnchorX = e.getX();
            dragAnchorY = e.getY();
        }

        @Override
        public void mouseMoved(MouseEvent e) {
            dragAnchorX = e.getX();
            dragAnchorY = e.getY();
        }

        @Override
        public void mouseWheelMoved(MouseWheelEvent e) {

            zoom += (.1 * e.getWheelRotation());
            zoom = Math.max(0.1, zoom);
            zoom = Math.min(zoom, 10);
            originX = e.getX();
            originY = e.getY();

            repaint();
        }
    }

}
