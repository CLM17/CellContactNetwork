/*
 * To the extent possible under law, the ImageJ developers have waived
 * all copyright and related or neighboring rights to this tutorial code.
 *
 * See the CC0 1.0 Universal license for details:
 *     http://creativecommons.org/publicdomain/zero/1.0/
 */

package com.mycompany.imagej;

import java.util.Iterator;
import ij.process.ImageProcessor;
import ij.ImagePlus;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.TreeSet;
import ij.IJ;
import net.imagej.ops.OpService;
import org.scijava.ui.UIService;
import net.imagej.Dataset;
import org.scijava.plugin.Parameter;
import java.io.File;
import org.scijava.plugin.Plugin;
import org.scijava.command.Command;
import net.imglib2.type.numeric.RealType;


/**
 * This example illustrates how to create an ImageJ {@link Command} plugin.
 * <p>
 * The code here is a simple Gaussian blur using ImageJ Ops.
 * </p>
 * <p>
 * You should replace the parameter fields with your own inputs and outputs,
 * and replace the {@link run} method implementation with your own logic.
 * </p>
 */
@Plugin(type = Command.class, menuPath = "Plugins>Network Creator")
public class NetworkCreator<T extends RealType<T>> implements Command
{
    @Parameter(label = "Please select directory for output file", style = "directory")
    private File outputdir;
    @Parameter
    private Dataset currentData;
    @Parameter
    private UIService uiService;
    @Parameter
    private OpService opService;
    
    public void run() {
        final ImagePlus imp = IJ.getImage();
        final TreeSet<LabelPair> list = new TreeSet<LabelPair>();
        final ImageProcessor ip = imp.getChannelProcessor();
        
        // Loop through image, but don't include top and left edges.
        System.out.println("Entering for-loop");
        for (int x = 2; x < ip.getWidth(); ++x) {
            for (int y = 2; y < ip.getHeight(); ++y) {
            	
                final int center = (int)ip.getf(x, y);
                
                if (center > 0) {
                    final int north = (int)ip.getf(x, y - 2);
                    final int west = (int)ip.getf(x - 2, y);
                    final int northwest = (int)ip.getf(x - 2, y - 2);
                    
                    if (north > 0 && center != north) {
                        final LabelPair pairN = new LabelPair(center, north);
                        if (!list.contains(pairN)) {
                            list.add(pairN);
                        }
                    }
                    if (west > 0 && center != west) {
                        final LabelPair pairW = new LabelPair(center, west);
                        if (!list.contains(pairW)) {
                            list.add(pairW);
                        }
                    }
                    if (northwest > 0 && center != northwest) {
                        final LabelPair pairNW = new LabelPair(center, northwest);
                        if (!list.contains(pairNW)) {
                            list.add(pairNW);
                        }
                    }
                }
            }
        }
        
        PrintWriter pw = null;
        try {
            pw = new PrintWriter(new File(this.outputdir + "/edges.csv"));
        }
        catch (FileNotFoundException e) {
            e.printStackTrace();
        }
        for (final LabelPair pair : list) {
            pw.print(String.valueOf(pair.label1) + "," + pair.label2);
            pw.print("\n");
        }
        pw.close();
    }
    
    public final class LabelPair implements Comparable<LabelPair>
    {
        public int label1;
        public int label2;
        
        public LabelPair(final int label1, final int label2) {
            if (label1 < label2) {
                this.label1 = label1;
                this.label2 = label2;
            }
            else {
                this.label1 = label2;
                this.label2 = label1;
            }
        }
        
        @Override
        public int compareTo(final LabelPair pair) {
            if (this.label1 < pair.label1) {
                return -1;
            }
            if (this.label1 > pair.label1) {
                return 1;
            }
            if (this.label2 < pair.label2) {
                return -1;
            }
            if (this.label2 > pair.label2) {
                return 1;
            }
            return 0;
        }
    }
}
