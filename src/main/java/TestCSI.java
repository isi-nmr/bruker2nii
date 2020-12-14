import bruker_plugin_lib.Bruker;
import org.apache.commons.lang3.ArrayUtils;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.factory.Nd4j;
import org.nd4j.linalg.indexing.INDArrayIndex;
import org.nd4j.linalg.indexing.NDArrayIndex;
import org.nd4j.linalg.inverse.InvertMatrix;

import java.util.Arrays;
import java.util.List;
import java.util.stream.IntStream;

public class TestCSI {
    public static void main(String[] args) {
        Bruker2nii bruker2nii = new Bruker2nii("D:\\DATA SETs\\for test Jbruker\\20200612_094625_lego_phantom_3_1_2\\30\\fid");
        bruker2nii.convert("testCSI.nii.gz");

//
//        Bruker brukerCSI = new Bruker();
//
//        String pathTempcsi = "D:\\DATA SETs\\for test Jbruker\\20200612_094625_lego_phantom_3_1_2\\30\\pdata\\1\\2dseq";
//        brukerCSI.setPath(Paths.get(pathTempcsi));
//        DataBruker dataCSI = brukerCSI.getData();
//        double[][][][] CSI2nii = new double[2][16][16][2048];
//
//
//        NiftiHeader header = new NiftiHeader(CSI2nii.length/2, CSI2nii[0].length, CSI2nii[0][0].length, (int) CSI2nii[0][0][0].length);
//        header.datatype = NiftiHeader.NIFTI_TYPE_COMPLEX64;
//        NiftiVolume volumeW = new NiftiVolume(header);
//        volumeW.header.pixdim[4] = (float) (1/brukerCSI.getJcampdx().getSW(4000));
//        volumeW.header.xyzt_units = NiftiHeader.NIFTI_UNITS_SEC | NiftiHeader.NIFTI_UNITS_METER;
//
//        // Qform
////        INDArray affinMat = getAffineMat(brukerCSI);
//        volumeW.header.qform_code = NiftiHeader.NIFTI_XFORM_ALIGNED_ANAT;
//        float quatern_b = 1;
//        float quatern_c = 1;
//        float quatern_d = 1;
//        volumeW.header.quatern = new float[]{quatern_b, quatern_c, quatern_d};
//        volumeW.header.qfac = 1;
//        //qform end
//
//        volumeW.header.intent_code = 2001; // each voxel has time_series
//
//
//        for (int i = 0; i < 16; i++) {
//            for (int j = 0; j < 16 ; j++) {
//                INDArrayIndex[] indx = { NDArrayIndex.all() , NDArrayIndex.point(i), NDArrayIndex.point(j)};
////                CSI2nii[0][i][j] = dataCSI.real.get(indx).toDoubleVector();
////                CSI2nii[1][i][j] = dataCSI.imag.get(indx).toDoubleVector();
//                for (int t = 0; t < 2048 ; t++) {
//                    volumeW.data.set(j,i,0,t, CSI2nii[0][j][i][t]);
//                    volumeW.data.set(j,i,0,t, CSI2nii[1][j][i][t]);
//                }
//
//            }
//        }
//
//        try {
//            volumeW.write("testCSI.nii.gz");
//        } catch (IOException e) {
//            e.printStackTrace();
//        }
    }

    private static INDArray getAffineMat(Bruker brukerCSI) {
        INDArray visuCoreSize = brukerCSI.getJcampdx().getVisu_pars().getINDArray("VisuCoreSize");
        INDArray visuCoreExtent = brukerCSI.getJcampdx().getVisu_pars().getINDArray("VisuCoreExtent");
        INDArray resolution = visuCoreExtent.div(visuCoreSize);
        double VisuCoreFrameThickness = brukerCSI.getJcampdx().getVisu_pars().getFloat("VisuCoreFrameThickness");
//        INDArray spatialResolution = Nd4j.append(resolution.get(NDArrayIndex.interval(1, 3)), 1, VisuCoreFrameThickness, -1);
        INDArray spatialResolution = Nd4j.append(resolution, 1, VisuCoreFrameThickness, -1);
        INDArray visuCoreOrientation = brukerCSI.getJcampdx().getVisu_pars().getINDArray("VisuCoreOrientation");
        INDArray visuCorePosition = brukerCSI.getJcampdx().getVisu_pars().getINDArray("VisuCorePosition");
        String visuSubjectPosition = brukerCSI.getJcampdx().getVisu_pars().getString("VisuSubjectPosition");
        INDArray pvm_SPackArrGradOrient = brukerCSI.getJcampdx().getMethod().getINDArray("PVM_SPackArrGradOrient");
        INDArray visuCoreUnits = brukerCSI.getJcampdx().getVisu_pars().getINDArray("VisuCoreUnits");
        // to do check mm
        // add frame count to visuCoreSize (just for image)
        Float VisuCoreDataSlope = brukerCSI.getJcampdx().getVisu_pars().getFloat("VisuCoreDataSlope");
        Float VisuCoreDataOffs = brukerCSI.getJcampdx().getVisu_pars().getFloat("VisuCoreDataOffs");
        // slope Correction
        // offset Correction
        String[] subPos = new String[] {"Head_Prone", "Head_Supine"};
        if(Arrays.stream(subPos).anyMatch(e -> e.contains(visuSubjectPosition))){
            System.out.println("Known case ('Head_Prone' or  'Head_Supine' for the parameter 'visu_pars.VisuSubjectPosition.)");
        }
        INDArray visuCoreOrientationReShaped = visuCoreOrientation.get(NDArrayIndex.point(0),NDArrayIndex.all()).reshape('f', new int[]{3, 3});
        INDArray result = Nd4j.eye(4);
        result.put(new INDArrayIndex[] {NDArrayIndex.interval(0,3), NDArrayIndex.interval(0,3)}, visuCoreOrientationReShaped);
        result.put(new INDArrayIndex[] {NDArrayIndex.interval(0,3), NDArrayIndex.point(3)}, visuCorePosition.getRows(0));
//                  Invert the orientation matrix, according to nifti convention and Bruker manual.
        InvertMatrix.invert(result, true);
        // calculate determinant
        float[][] defMat = new float[][]{{1, 0, 0}, {0, 0, 1}, {0, 1, 0}};
        INDArray result_orientation =
                result.get(new INDArrayIndex[]{NDArrayIndex.interval(0, 3), NDArrayIndex.interval(0, 3)}).mmul(Nd4j.createFromArray(defMat));
//        from SAR to ASL
        if (pivot(result_orientation.getColumn(0).toDoubleVector())>0) {
            result_orientation.getColumn(0).muli(-1);
        }
        if (pivot(result_orientation.getColumn(1).toDoubleVector())>0) {
            result_orientation.getColumn(1).muli(-1);
        }
        if (pivot(result_orientation.getColumn(2).toDoubleVector())<0) {
            result_orientation.getColumn(2).muli(-1);
        }
        result_orientation.mmuli(Nd4j.diag(spatialResolution));
        result.put(new INDArrayIndex[] {NDArrayIndex.interval(0,3), NDArrayIndex.interval(0,3)}, result_orientation);
        // consider_subject_position:
        // keep same det
        return result;
    }
    public static Double pivot(double[] arr) {
        Double[] doubleArray = ArrayUtils.toObject(arr);
        List<Double> arrList = Arrays.asList(doubleArray);
        int indx = IntStream.range(0, arr.length).reduce((a,b) -> ((double) Math.abs(arrList.get(a))) < ((double) Math.abs(arrList.get(b))) ? b: a).getAsInt();
        return arr[indx];
    }
    public static INDArray calc_eulerangle(INDArray mati) {
        INDArray mat = mati.get(NDArrayIndex.point(0), NDArrayIndex.all(), NDArrayIndex.all());
        double sy = Math.sqrt(mat.getDouble(0, 0) * mat.getDouble(0, 0) + mat.getDouble(1, 0) * mat.getDouble(1, 0));
        double x, y, z;
        if(sy < 1e-6) {
            x = Math.atan2(-mat.getDouble(1,2), mat.getDouble(1,1));
            y = Math.atan2(-mat.getDouble(2,0), sy);
            z = 0;
        } else {
            x = Math.atan2(mat.getDouble(2,1), mat.getDouble(2,2));
            y = Math.atan2(-mat.getDouble(2,0), sy);
            z = Math.atan2(mat.getDouble(1,0), mat.getDouble(0,0));
        }
        return Nd4j.create(new double[] {x, y, z});
    }
}

