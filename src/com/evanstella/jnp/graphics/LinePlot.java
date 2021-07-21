package com.evanstella.jnp.graphics;

import com.evanstella.jnp.core.IllegalDimensionException;
import com.evanstella.jnp.core.Numeric;

import javax.swing.JComponent;
import java.awt.*;

public class LinePlot extends JComponent {

    private Numeric XData;
    private Numeric YData;

    private int scale;

    public LinePlot ( Numeric XData, Numeric YData ) {
        super();
        this.XData = XData;
        this.YData = YData;
    }

    public void paintComponent (Graphics g) {
        update( g );
    }

    public void update ( Graphics g ) {
        //w is x, and h is y (as in x/y values in a graph)
        int w = this.getWidth()/2;
        int h = this.getHeight()/2;


        Graphics2D g1 = (Graphics2D) g;
        g1.setRenderingHint(
                RenderingHints.KEY_ANTIALIASING,
                RenderingHints.VALUE_ANTIALIAS_ON);
        g1.setStroke(new BasicStroke(2));
        g1.setColor(Color.black);
        g1.drawLine(0,h,w*2,h);
        g1.drawLine(w,0,w,h*2);
        g1.drawString("0", w - 7, h + 13);


        Graphics2D g2 = (Graphics2D) g;
        g2.setStroke(new BasicStroke((float)1.5));
        g2.setColor(Color.red);
        //line1
        Polygon p = new Polygon();
        scale = 20;
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

}
