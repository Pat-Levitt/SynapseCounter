import ij.*;
import ij.process.*;
import ij.gui.*;
import java.awt.*;
import java.util.*;
import ij.plugin.*;
import ij.plugin.frame.*;
import ij.plugin.filter.*;
import ij.measure.*;
import ij.gui.Roi;
import ij.gui.ShapeRoi;
import java.util.Hashtable;
import java.util.Enumeration;
import java.io.*;
import java.util.Arrays;



public class SynapseCounter_v1 implements PlugIn {

    	static String title="Example";
    	static int width=512,height=512;
    	static int minPSD95=80,maxPSD95=255,minGFP=60,maxGFP=255,minSynI=120;
    	static double minArea=0.06,maxArea=2;
    
   
	public void run(String arg) {
//Build the dialog box to input the threshods and particle size.//
      		GenericDialog gd = new GenericDialog("Set Thresholds");
      		gd.addNumericField("minthreshold for PSD95 channel: ", minPSD95, 0);
      		gd.addNumericField("minthreshod for GFP channel: ", minGFP,0);
      		gd.addNumericField("minthreshod for SynI channel: ", minSynI,0);
      		gd.addNumericField("minimal particle size: ", minArea,4);
      		gd.addNumericField("maximal particle size: ", maxArea,4);
      		gd.showDialog();
      		if (gd.wasCanceled()) return;
      		minPSD95 = (int)gd.getNextNumber();
      		minGFP = (int)gd.getNextNumber();
      		minSynI = (int)gd.getNextNumber();
      		minArea = gd.getNextNumber();
      		maxArea = gd.getNextNumber();
//Store the thresholds information.
      		IJ.log("Thresholds:");
      	    	IJ.log("minPSD95: "+minPSD95+" maxPSD95: "+maxPSD95+" minGFP: "+minGFP+" maxGFP: "+maxGFP+" minSynI: "+minSynI);
      	    	IJ.log("particle size: "+minArea+" "+maxArea);
//read images, channels have to be splitted before running the plugin, three images representating three channels are read once together.//
		String dir = IJ.getDirectory("Choose a Folder");
		if (dir==null) return;
		String[] list = (new File(dir)).list();
		if (list==null) return;
		Analyzer.setMeasurement(Measurements.LABELS, true);    		
		for (int i=0; i<list.length/3; i++) {
			if (list[i].startsWith(".")) continue;
          		String path1 = dir+list[3*i];
          		String path2 = dir+list[3*i+1];
			String path3 = dir+list[3*i+2];  
			IJ.showProgress(i+1, list.length);
			ImagePlus imp1 = IJ.openImage(path1);
		        ImagePlus imp2 = IJ.openImage(path2);			
			ImagePlus imp3 = IJ.openImage(path3);						
			ImageConverter ic1 = new ImageConverter(imp1);
			ImageConverter ic2 = new ImageConverter(imp2);
			ImageConverter ic3 = new ImageConverter(imp3);
//Images have to be 8-bit//			
		        if (imp1.getType() != ImagePlus.GRAY8) {
            			ic1.convertToGray8(); 
			        IJ.showMessage("Error", "8-Bit Grayscale Image Required");	
			}
		        if (imp2.getType() != ImagePlus.GRAY8) {
            			ic2.convertToGray8(); 
			        IJ.showMessage("Error", "8-Bit Grayscale Image Required");	
			}
		        if (imp3.getType() != ImagePlus.GRAY8) {
            			ic3.convertToGray8(); 
			        IJ.showMessage("Error", "8-Bit Grayscale Image Required");
			}
//analyze the image of PSD95 channel, apply threshold and make selection, add the selection to ROI manager, split the selection//									
			ImageProcessor ip = imp3.getProcessor();
			ThresholdToSelection tts = new ThresholdToSelection();
        		RoiManager roim= new RoiManager();
        		java.awt.List listRoi=roim.getList();
        		Hashtable rois=roim.getROIs();
        		Analyzer analyzer=new Analyzer(imp3);
    			int lutupdate=2;
    			ip.setThreshold(minPSD95, maxPSD95, lutupdate);
        		Roi roi=tts.convert(ip);
        		WindowManager.setTempCurrentImage(imp3);
		  	roim.add(imp3,roi,0);
  			roim.select(0);
       			roim.runCommand("split");
       			roim.runCommand("delete");
//apply the splitted selections from PSD95 channel to GFP channel, if the selection larger than 4 pixel and at least one pixel with intensity larger than// 
//minimal threhold of GFP channel, then keep this selection; otherwise, remove the selection.//       			
			WindowManager.setTempCurrentImage(imp1);
		        roim.select(imp1,-1);
        		roim.runCommand("Measure");
        		ResultsTable measureResult=Analyzer.getResultsTable();        
        		int count = listRoi.getItemCount();
        		int measureCount=measureResult.getCounter();
        		for (int j=count-1;j>=0;j--){
        			if (measureResult.getValue("Area",j)<=minArea||measureResult.getValue("Area",j)>=maxArea||measureResult.getValue("Max",j)<=minGFP){
        			measureResult.deleteRow(j);
        			rois.remove(listRoi.getItem(j));
				listRoi.remove(j);
       				}
		        }
		        measureResult.reset();
		        int count2 = listRoi.getItemCount();
//apply remained selections to SynI channel, if the selection has at least one pixel with intensity larger than minimal threshold of SynI channel, keep this selection
//otherwise, remove the selection//
		        WindowManager.setTempCurrentImage(imp2);
		        roim.select(imp2,-1);
		        roim.runCommand("Measure");
		        for (int k=count2-1;k>=0;k--){
		        	if (measureResult.getValue("Max",k)<=minSynI){
        			measureResult.deleteRow(k);
        			rois.remove(listRoi.getItem(k));
				listRoi.remove(k);
       				}        	
		        }
		        int count3 = listRoi.getItemCount();
//save the measurement and rois
		        try{
		        measureResult.saveAs(path3+".txt");
		        }catch (IOException e)
		        {IJ.showMessage("Error");
		        }
		        measureResult.reset();
		        roim.runCommand("save", path3+".zip");
//measure the area of dendrites, GFP area after thresholding
		 	WindowManager.setTempCurrentImage(imp1);
		        ImageProcessor ipGFP = imp1.getProcessor();
    			ipGFP.setThreshold(minGFP, maxGFP, lutupdate);
    			ThresholdToSelection ttsGFP = new ThresholdToSelection();
    			Roi roiGFP=ttsGFP.convert(ipGFP);
    			roim.add(imp1,roiGFP,0);
  			roim.runCommand("Measure");
    			ResultsTable measureResultGFP=Analyzer.getResultsTable();
    			IJ.log(path3+"R&B Channel:"+count2+" R&G&B Channel: "+count3+"GFPArea: "+measureResultGFP.getValue("Area",count3));
    			measureResult.reset();
    			roim.close(); 
		}
	}
}				      