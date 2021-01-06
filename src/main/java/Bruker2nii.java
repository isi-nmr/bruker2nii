import bruker_plugin_lib.Bruker;
import bruker_plugin_lib.DataBruker;
import checkers.units.quals.A;
import com.ericbarnhill.niftijio.NiftiHeader;
import com.ericbarnhill.niftijio.tools.IndexIterator;
import org.apache.commons.lang3.ArrayUtils;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.factory.Nd4j;
import org.nd4j.linalg.indexing.INDArrayIndex;
import org.nd4j.linalg.indexing.NDArrayIndex;
import org.nd4j.linalg.inverse.InvertMatrix;

import java.io.IOException;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.stream.IntStream;

import static org.nd4j.linalg.ops.transforms.Transforms.sqrt;
import static org.nd4j.linalg.ops.transforms.Transforms.*;
public class Bruker2nii {
    String pathBruker;
    String pathNii;
    Bruker bruker;
    DataBruker dataBruker;
    NiftiMRS niftiMRS;

    public Bruker2nii(String pathBruker) {
        this.pathBruker = pathBruker;
        bruker = new Bruker();
        bruker.setPath(Paths.get(pathBruker));
        dataBruker = bruker.getData();


    }
    public void convert(String s, char d) {
        convert(s, d, false, false);
    }
    public void convert(String s) {
        convert(s, 't', false, false);
    }

    public boolean convert(String pathNii, char d, boolean b, boolean b1) {
        this.pathNii = pathNii;
        if (isSVS()) {
            boolean squeez = false;
            int brukerLength = bruker.getDims().length;
            if (bruker.getDims()[brukerLength-1] == 1)
                squeez = true;
            if (squeez)
                brukerLength -=1;

            int[] niftishape = new int[3 + brukerLength];
            int[] shape = Arrays.stream(bruker.getDims()).mapToInt(i -> (int) i).toArray();
            System.arraycopy(new int[]{1, 1, 1}, 0, niftishape, 0, 3);
            if (squeez)
                System.arraycopy(shape, 0, niftishape, 3, bruker.getDims().length-1);
            else System.arraycopy(shape, 0, niftishape, 3, bruker.getDims().length);

            niftiMRS = new NiftiMRS(niftishape);

            INDArray PVM_VoxArrSize = bruker.getJcampdx().getMethod().getINDArray("PVM_VoxArrSize");
            for (int i = 1; i < 4; i++) { // what is 4?
                niftiMRS.getNifti().header.pixdim[i] = PVM_VoxArrSize.getFloat(i - 1);
            }
            // here pixdim 5 to 7 must be set
            niftiMRS.getNifti().header.pixdim[4] = bruker.getJcampdx().getVisu_pars().getINDArray("VisuCoreExtent").getFloat(0)
                    / bruker.getJcampdx().getVisu_pars().getINDArray("VisuCoreSize").getFloat(0);
            ArrayList<int[]> idcs;
            if (squeez)
                idcs = new IndexIterator().iterate(Arrays.copyOfRange(shape, 1, shape.length-1));
            else idcs = new IndexIterator().iterate(Arrays.copyOfRange(shape, 1, shape.length));

            for (int[] idc : idcs) {
                int[] idx = new int[niftishape.length];
                System.arraycopy(idc, 0, idx, 4, idc.length);
                for (int i = 0; i < shape[0]; i++) {
                    idx[3] = i;
                    idx[0] = 0;
                    int[] dataidx = Arrays.copyOfRange(idx, 3, idx.length);
                    niftiMRS.getNifti().data.set(idx, dataBruker.real.getFloat(dataidx));
                    idx[0] = 1;
                    niftiMRS.getNifti().data.set(idx, dataBruker.imag.getFloat(dataidx));
                }
            }
            String NUCLEUS = bruker.getJcampdx().getAcqp().getString("NUCLEUS");
            if (NUCLEUS.contains("1H"))
                niftiMRS.getJson().ResonantNucleus = new String[]{Nucleus.N_1H.toString()};
            Float TrnsFreq = bruker.getJcampdx().getAcqp().getFloat("SFO1");
            niftiMRS.getJson().SpectrometerFrequency = new Double[]{Double.valueOf(TrnsFreq)};

            if (!squeez) { // what is it really!? avg or coil
            niftiMRS.getJson().setDim_5(DIM_KEYS.DIM_MEAS);
                niftiMRS.getJson().setDim_5_info("repetition");
                niftiMRS.getJson().setDim_5_header(TAGS.EchoTime, Collections.singletonList(bruker.getJcampdx().getTE(10000)));
                niftiMRS.getJson().setDim_5_header(TAGS.RepetitionTime, Collections.singletonList(bruker.getJcampdx().getTR(10000)));
            }
            niftiMRS.getNifti().header.descrip = new StringBuffer("Converted by JBruker2nii API");
        }
        else {

            int[] shape_nifti = new int[3 + 1];
            int[] shape_bruker = Arrays.stream(bruker.getDims()).mapToInt(i -> (int) i).toArray();
            shape_nifti[0] = shape_bruker[1];
            shape_nifti[1] = shape_bruker[2];
            if (shape_bruker.length == 4)
                shape_nifti[2] = shape_bruker[3];
            else
                shape_nifti[2] = 1;
            shape_nifti[3] = shape_bruker[0];
            niftiMRS = new NiftiMRS(shape_nifti);
            INDArray PVM_VoxArrSize = bruker.getJcampdx().getMethod().getINDArray("PVM_VoxArrSize");
            INDArray VisuCoreExtent = bruker.getJcampdx().getVisu_pars().getINDArray("VisuCoreExtent");
            INDArray VisuCoreSize = bruker.getJcampdx().getVisu_pars().getINDArray("VisuCoreSize");
            for (int i = 1; i < 3; i++) {
                niftiMRS.getNifti().header.pixdim[i] = VisuCoreExtent.getFloat(i) / VisuCoreSize.getFloat(i);
            }
            niftiMRS.getNifti().header.pixdim[3] = bruker.getJcampdx().getVisu_pars().getFloat("VisuCoreFrameThickness");
            niftiMRS.getNifti().header.pixdim[4] = VisuCoreExtent.getFloat(0) / VisuCoreSize.getFloat(0);
            ArrayList VisuCoreUnits = bruker.getJcampdx().getVisu_pars().getArrayList("VisuCoreUnits");
            double[][][] CSI2nii = new double[(int) bruker.getDims()[1]][(int) bruker.getDims()[2]][(int) bruker.getDims()[0]];
            NiftiHeader header = new NiftiHeader((int) bruker.getDims()[1], (int) bruker.getDims()[1], 1, (int) bruker.getDims()[0]);
            niftiMRS.getNifti().header.t_unit_code = NiftiHeader.NIFTI_UNITS_SEC;
            if (VisuCoreUnits.contains("mm"))
                niftiMRS.getNifti().header.xyz_unit_code = NiftiHeader.NIFTI_UNITS_MM;
            else if (VisuCoreUnits.contains("m"))
                niftiMRS.getNifti().header.xyz_unit_code = NiftiHeader.NIFTI_UNITS_METER;

                for (int i = 0; i < bruker.getDims()[1]; i++) {
                    for (int j = 0; j < bruker.getDims()[2]; j++) {
                        INDArrayIndex[] indx = {NDArrayIndex.all(), NDArrayIndex.point(i), NDArrayIndex.point(j/2)};
                        double[] tempR = dataBruker.real.get(indx).toDoubleVector();
                        double[] tempI = dataBruker.imag.get(indx).toDoubleVector();
                        int[] idx = new int[]{i, j, 0, 0};
                        for (int t = 0; t < bruker.getDims()[0]; t++) {
                            idx[3] = t;
                            idx[0] = 0;
                            niftiMRS.getNifti().data.set(idx, tempR[t]);
                            idx[0] = 1;
                            niftiMRS.getNifti().data.set(idx, tempI[t]);

                        }

                    }
                }

            INDArray affine_mat = getAffineMat();
            // if sample upside down
            //            q mat
            INDArray transpos = affine_mat.get(new INDArrayIndex[]{NDArrayIndex.interval(0, 3), NDArrayIndex.point(3)});
            INDArray RZS = affine_mat.get(new INDArrayIndex[]{NDArrayIndex.interval(0, 3), NDArrayIndex.interval(0, 3)});
            INDArray zoom = sqrt(RZS.mul(RZS).sum(0));
            INDArray rotate = RZS.divRowVector(zoom);
            // qfactor
//            int nRows = rotate.rows();
//            int nColumns = rotate.columns();
//            INDArray S = Nd4j.zeros(1, nRows);
//            INDArray U = Nd4j.zeros(nRows, nRows);
//            INDArray V = Nd4j.zeros(nColumns, nColumns);
//            Nd4j.getBlasWrapper().lapack().gesvd(rotate, S, U, V);
//            INDArray PR = U.mmul(V);
//            // error proning
//            PR.reshape(new int[]{1,9});
//          s mat
            niftiMRS.getNifti().header.srow_x = affine_mat.getRow(0).toFloatVector();
            niftiMRS.getNifti().header.srow_y = affine_mat.getRow(1).toFloatVector();
            niftiMRS.getNifti().header.srow_z = affine_mat.getRow(2).toFloatVector();
            String NUCLEUS = bruker.getJcampdx().getAcqp().getString("NUCLEUS");
            if (NUCLEUS.contains("1H"))
                niftiMRS.getJson().ResonantNucleus = new String[]{Nucleus.N_1H.toString()};
            Float TrnsFreq = bruker.getJcampdx().getAcqp().getFloat("SFO1");
            niftiMRS.getJson().SpectrometerFrequency = new Double[]{Double.valueOf(TrnsFreq)};
            niftiMRS.getNifti().header.descrip = new StringBuffer("Converted by JBruker2nii API");
            }
        try {
            niftiMRS.write(pathNii, b, b1);
        } catch (IOException e) {
            e.printStackTrace();
        }

        return true;






//        if (!bruker.isRaw()) {
//            if (isSVS()) {
////                int sigLen = bruker.getJcampdx().getVisu_pars().getINDArray("VisuAcqSize").getInt(0);
//                if (bruker.isIR()) {
//                    int sigLen = bruker.getJcampdx().getVisu_pars().getINDArray("VisuAcqSize").getInt(0);
//                    double[] data2nii = new double[sigLen];
//                    int[] niftishape = new int[3 + bruker.getDims().length];
//                    int[] shape = Arrays.stream(bruker.getDims()).mapToInt(i -> (int) i).toArray();
//                    System.arraycopy(new int[]{1, 1, 1}, 0, niftishape, 0, 3);
//                    System.arraycopy(shape, 0, niftishape, 3, bruker.getDims().length);
//
//                    niftiMRS = new NiftiMRS(niftishape);
//                    INDArray PVM_VoxArrSize = bruker.getJcampdx().getMethod().getINDArray("PVM_VoxArrSize");
//                    for (int i = 1; i < 4; i++) {
//                        niftiMRS.getNifti().header.pixdim[i] = PVM_VoxArrSize.getFloat(i - 1);
//                    }
//                    // here pixdim 5 to 7 must be set
//                    niftiMRS.getNifti().header.pixdim[4] = bruker.getJcampdx().getVisu_pars().getINDArray("VisuCoreExtent").getFloat(0)
//                            / bruker.getJcampdx().getVisu_pars().getINDArray("VisuCoreSize").getFloat(0);
//                    ArrayList<int[]> idcs = new IndexIterator().iterate(Arrays.copyOfRange(shape, 1, shape.length));
//                    for (int[] idc : idcs) {
//                        int[] idx = new int[niftishape.length];
//                        System.arraycopy(idc, 0, idx, 4, idc.length);
//                        for (int i = 0; i < shape[0]; i++) {
//                            idx[3] = i;
//                            idx[0] = 0;
//                            int[] dataidx = Arrays.copyOfRange(idx, 3, idx.length);
//                            niftiMRS.getNifti().data.set(idx, dataBruker.real.getFloat(dataidx));
//                            idx[0] = 1;
//                            niftiMRS.getNifti().data.set(idx, dataBruker.imag.getFloat(dataidx));
//                        }
//                    }
////                    double[][] data2nii = new double[2][sigLen];
////                    data2nii[0] = dataBruker.real.getColumn(0).toDoubleVector();
////                    data2nii[1] = dataBruker.imag.getColumn(0).toDoubleVector();
////                    // check the possibility of repitted measure!
////                    int[] niftishape = new int[]{1, 1, 1, (int) data2nii[0].length};
////                    niftiMRS = new NiftiMRS(niftishape);
////                    INDArray PVM_VoxArrSize = bruker.getJcampdx().getMethod().getINDArray("PVM_VoxArrSize");
////                    for (int i = 1; i < 4; i++) {
////                        niftiMRS.getNifti().header.pixdim[i] = PVM_VoxArrSize.getFloat(i - 1);
////                    }
////                    if (d == 'f') {
////                        niftiMRS.getNifti().header.pixdim[4] = bruker.getJcampdx().getVisu_pars().getINDArray("VisuCoreExtent").getFloat(0);
////                        niftiMRS.getNifti().header.xyz_unit_code = NiftiHeader.NIFTI_UNITS_MM;
////                        niftiMRS.getNifti().header.t_unit_code = NiftiHeader.NIFTI_UNITS_PPM;
//////                        for (int i = 0; i < data2nii[0].length; i++) {
//////                            for (int j = 0; j < 2; j++) {
////////                                volumeW.data.set(j, 0, 0, i, data2nii[j][i]);
//////                            }
//////                        }
////                        int[] shape = Arrays.stream(bruker.getDims()).mapToInt(i -> (int) i).toArray();
////                        ArrayList<int[]> idcs = new IndexIterator().iterate(Arrays.copyOfRange(shape, 1, shape.length));
////
////                        for (int[] idc : idcs) {
////                            int[] idx = new int[niftishape.length];
////                            System.arraycopy(idc, 0, idx, 4, idc.length);
////                            for (int i = 0; i < shape[0]; i++) {
////                                idx[3] = i;
////                                idx[0] = 0;
////                                int[] dataidx = Arrays.copyOfRange(idx, 3, idx.length);
////                                niftiMRS.getNifti().data.set(idx, dataBruker.real.getFloat(dataidx));
////                                idx[0] = 1;
////                                niftiMRS.getNifti().data.set(idx, dataBruker.imag.getFloat(dataidx));
////                            }
////                        }
//                        String NUCLEUS = bruker.getJcampdx().getAcqp().getString("NUCLEUS");
//                        if (NUCLEUS.contains("1H"))
//                            niftiMRS.getJson().ResonantNucleus = new String[]{Nucleus.N_1H.toString()};
//                        Float TrnsFreq = bruker.getJcampdx().getAcqp().getFloat("SFO1");
//                        niftiMRS.getJson().SpectrometerFrequency = new Double[]{Double.valueOf(TrnsFreq)};
////                volumeW.header.scl_slope = bruker.getJcampdx().getVisu_pars().getINDArray("VisuCoreDataSlope").getFloat(0);
////                volumeW.header.scl_inter = bruker.getJcampdx().getVisu_pars().getINDArray("VisuCoreDataOffs").getFloat(0);
//
////                volumeW.header.cal_max = bruker.getJcampdx().getVisu_pars().getFloat("VisuCoreDataMax");
////                volumeW.header.cal_min = bruker.getJcampdx().getVisu_pars().getFloat("VisuCoreDataMin");
//                        niftiMRS.getJson().setDim_5(DIM_KEYS.DIM_MEAS);
//                        niftiMRS.getJson().setDim_5_info("repetition");
//                        niftiMRS.getJson().setDim_5_header(TAGS.EchoTime, Collections.singletonList(bruker.getJcampdx().getTE(10000)));
//                        niftiMRS.getJson().setDim_5_header(TAGS.RepetitionTime, Collections.singletonList(bruker.getJcampdx().getTR(10000)));
//
//                        niftiMRS.getNifti().header.descrip = new StringBuffer("Converted by JBruker2nii API");
//                    }
////                    else if (d == 't') {
////                        FastFourierTransformer fft = new FastFourierTransformer(DftNormalization.STANDARD);
////                        Complex[] complexdata = TransformUtils.createComplexArray(data2nii);
////                        Complex[] timedata = fft.transform(complexdata, TransformType.INVERSE);
////                        data2nii = TransformUtils.createRealImaginaryArray(timedata);
////                        niftiMRS.getNifti().header.xyz_unit_code = NiftiHeader.NIFTI_UNITS_MM;
////                        niftiMRS.getNifti().header.t_unit_code = NiftiHeader.NIFTI_UNITS_SEC;
////                        niftiMRS.getNifti().header.pixdim[4] = (float) (1 / bruker.getJcampdx().getSW(4000));
////                        for (int i = 0; i < data2nii[0].length; i++) {
////                            for (int j = 0; j < 2; j++) {
////                                niftiMRS.getNifti().data.set(j, 0, 0, i, data2nii[j][i]);
////                            }
////                        }
////                    }
//                }
////                else {
////                    double[] data2nii = new double[sigLen];
////                    data2nii = dataBruker.real.getColumn(0).toDoubleVector();
//////                data2nii[1]  = dataBruker.imag.getColumn(0).toDoubleVector();
////                    niftiHeader = new NiftiHeader(1, 1, 1, (int) data2nii.length);
////                    niftiHeader.datatype = NiftiHeader.NIFTI_TYPE_FLOAT64;
////                    volumeW = new NiftiVolume(niftiHeader);
////                    INDArray PVM_VoxArrSize = bruker.getJcampdx().getMethod().getINDArray("PVM_VoxArrSize");
////                    for (int i = 1; i < 4; i++) {
////                        volumeW.header.pixdim[i] = PVM_VoxArrSize.getFloat(i - 1);
////                    }
////                    volumeW.header.pixdim[4] = bruker.getJcampdx().getVisu_pars().getINDArray("VisuCoreExtent").getFloat(0)
////                            / bruker.getJcampdx().getVisu_pars().getINDArray("VisuCoreSize").getFloat(0);
////                    for (int i = 0; i < data2nii.length; i++) {
//////                        volumeW.data.set(0, 0, 0, i, data2nii[i]);
////                    }
////
////                }
////                volumeW.header.scl_slope = bruker.getJcampdx().getVisu_pars().getINDArray("VisuCoreDataSlope").getFloat(0);
////                volumeW.header.scl_inter = bruker.getJcampdx().getVisu_pars().getINDArray("VisuCoreDataOffs").getFloat(0);
////
////                volumeW.header.cal_max = bruker.getJcampdx().getVisu_pars().getFloat("VisuCoreDataMax");
////                volumeW.header.cal_min = bruker.getJcampdx().getVisu_pars().getFloat("VisuCoreDataMin");
////                volumeW.header.descrip = new StringBuffer("Converted by JBruker2nii API");
//// to do affine matrix + 1r+1d + time
//
////            } else {
////
////                // if ( 3d and 4d )
////                double[][][] CSI2nii = new double[(int) bruker.getDims()[1]][(int) bruker.getDims()[2]][(int) bruker.getDims()[0]];
////                NiftiHeader header = new NiftiHeader((int) bruker.getDims()[1], (int) bruker.getDims()[1], 1, (int) bruker.getDims()[0]);
////                header.datatype = NiftiHeader.NIFTI_TYPE_COMPLEX64;
////                volumeW = new NiftiVolume(header);
////
////                INDArray VisuCoreExtent = bruker.getJcampdx().getVisu_pars().getINDArray("VisuCoreExtent");
////                INDArray VisuCoreSize = bruker.getJcampdx().getVisu_pars().getINDArray("VisuCoreSize");
////
////                for (int i = 1; i < 3; i++) {
////                    volumeW.header.pixdim[i] = VisuCoreExtent.getFloat(i) / VisuCoreSize.getFloat(i);
////                }
////                volumeW.header.pixdim[3] = bruker.getJcampdx().getVisu_pars().getFloat("VisuCoreFrameThickness");
////                volumeW.header.pixdim[4] = VisuCoreExtent.getFloat(0) / VisuCoreSize.getFloat(0);
////                ArrayList VisuCoreUnits = bruker.getJcampdx().getVisu_pars().getArrayList("VisuCoreUnits");
////
////                if (VisuCoreUnits.contains("[ppm]"))
////                    volumeW.header.t_unit_code = NiftiHeader.NIFTI_UNITS_PPM;
////                if (VisuCoreUnits.contains("mm"))
////                    volumeW.header.xyz_unit_code = NiftiHeader.NIFTI_UNITS_MM;
////                else if (VisuCoreUnits.contains("m"))
////                    volumeW.header.xyz_unit_code = NiftiHeader.NIFTI_UNITS_METER;
////
////
////                volumeW.header.intent_code = 2001; // each voxel has time_series
////                double[] tempR;
////                double[] tempI;
////                for (int i = 0; i < bruker.getDims()[1]; i++) {
////                    for (int j = 0; j < bruker.getDims()[2]; j=j+1) {
////                        INDArrayIndex[] indx = {NDArrayIndex.all(), NDArrayIndex.point(i), NDArrayIndex.point(j)};
////                        tempR = dataBruker.real.get(indx).toDoubleVector();
////
////                        for (int t = 0; t < bruker.getDims()[0]; t++) {
//////                            volumeW.data.set(j, i, 0, t, tempR[t]);
////                        }
////
////                    }
////                }
//
//        } else if (bruker.isRaw()) {
//            if (isSVS()) {
//                int sigLen = bruker.getJcampdx().getVisu_pars().getINDArray("VisuAcqSize").getInt(0);
//                double[] data2nii = new double[sigLen];
//                int[] niftishape = new int[3 + bruker.getDims().length];
//                int[] shape = Arrays.stream(bruker.getDims()).mapToInt(i -> (int) i).toArray();
//                System.arraycopy(new int[]{1, 1, 1}, 0, niftishape, 0, 3);
//                System.arraycopy(shape, 0, niftishape, 3, bruker.getDims().length);
//
//                niftiMRS = new NiftiMRS(niftishape);
//                INDArray PVM_VoxArrSize = bruker.getJcampdx().getMethod().getINDArray("PVM_VoxArrSize");
//                for (int i = 1; i < 4; i++) {
//                    niftiMRS.getNifti().header.pixdim[i] = PVM_VoxArrSize.getFloat(i - 1);
//                }
//                // here pixdim 5 to 7 must be set
//                niftiMRS.getNifti().header.pixdim[4] = bruker.getJcampdx().getVisu_pars().getINDArray("VisuCoreExtent").getFloat(0)
//                        / bruker.getJcampdx().getVisu_pars().getINDArray("VisuCoreSize").getFloat(0);
//                ArrayList<int[]> idcs = new IndexIterator().iterate(Arrays.copyOfRange(shape, 1, shape.length));
//                for (int[] idc : idcs) {
//                    int[] idx = new int[niftishape.length];
//                    System.arraycopy(idc, 0, idx, 4, idc.length);
//                    for (int i = 0; i < shape[0]; i++) {
//                        idx[3] = i;
//                        idx[0] = 0;
//                        int[] dataidx = Arrays.copyOfRange(idx, 3, idx.length);
//                        niftiMRS.getNifti().data.set(idx, dataBruker.real.getFloat(dataidx));
//                        idx[0] = 1;
//                        niftiMRS.getNifti().data.set(idx, dataBruker.imag.getFloat(dataidx));
//                    }
//                }
//                String NUCLEUS = bruker.getJcampdx().getAcqp().getString("NUCLEUS");
//                if (NUCLEUS.contains("1H"))
//                    niftiMRS.getJson().ResonantNucleus = new String[]{Nucleus.N_1H.toString()};
//                Float TrnsFreq = bruker.getJcampdx().getAcqp().getFloat("SFO1");
//                niftiMRS.getJson().SpectrometerFrequency = new Double[]{Double.valueOf(TrnsFreq)};
////                volumeW.header.scl_slope = bruker.getJcampdx().getVisu_pars().getINDArray("VisuCoreDataSlope").getFloat(0);
////                volumeW.header.scl_inter = bruker.getJcampdx().getVisu_pars().getINDArray("VisuCoreDataOffs").getFloat(0);
//
////                volumeW.header.cal_max = bruker.getJcampdx().getVisu_pars().getFloat("VisuCoreDataMax");
////                volumeW.header.cal_min = bruker.getJcampdx().getVisu_pars().getFloat("VisuCoreDataMin");
//                niftiMRS.getJson().setDim_5(DIM_KEYS.DIM_MEAS);
//                niftiMRS.getJson().setDim_5_info("repetition");
//                niftiMRS.getJson().setDim_5_header(TAGS.EchoTime, Collections.singletonList(bruker.getJcampdx().getTE(10000)));
//                niftiMRS.getJson().setDim_5_header(TAGS.RepetitionTime, Collections.singletonList(bruker.getJcampdx().getTR(10000)));
//
//                niftiMRS.getNifti().header.descrip = new StringBuffer("Converted by JBruker2nii API");
//// to do affine matrix + 1r+1d + time
//            } else {
//                // if ( 3d and 4d )
////                double[][][] CSI2nii = new double[(int) bruker.getDims()[1]][(int) bruker.getDims()[2]][(int) bruker.getDims()[0]];
////                NiftiHeader header = new NiftiHeader((int) bruker.getDims()[1], (int) bruker.getDims()[1], 1, (int) bruker.getDims()[0]);
////                header.datatype = NiftiHeader.NIFTI_TYPE_COMPLEX64;
////                volumeW = new NiftiVolume(header);
////
////                INDArray VisuCoreExtent = bruker.getJcampdx().getVisu_pars().getINDArray("VisuCoreExtent");
////                INDArray VisuCoreSize = bruker.getJcampdx().getVisu_pars().getINDArray("VisuCoreSize");
////
////                for (int i = 1; i < 3; i++) {
////                    volumeW.header.pixdim[i] = VisuCoreExtent.getFloat(i) / VisuCoreSize.getFloat(i);
////                }
////                volumeW.header.pixdim[3] = bruker.getJcampdx().getVisu_pars().getFloat("VisuCoreFrameThickness");
////                volumeW.header.pixdim[4] = VisuCoreExtent.getFloat(0) / VisuCoreSize.getFloat(0);
////                ArrayList VisuCoreUnits = bruker.getJcampdx().getVisu_pars().getArrayList("VisuCoreUnits");
////
////                if (VisuCoreUnits.contains("[ppm]"))
////                    volumeW.header.t_unit_code = NiftiHeader.NIFTI_UNITS_PPM;
////                if (VisuCoreUnits.contains("mm"))
////                    volumeW.header.xyz_unit_code = NiftiHeader.NIFTI_UNITS_MM;
////                else if (VisuCoreUnits.contains("m"))
////                    volumeW.header.xyz_unit_code = NiftiHeader.NIFTI_UNITS_METER;
////
////
////                volumeW.header.intent_code = 2001; // each voxel has time_series
////                double[] tempR;
////                double[] tempI;
////                for (int i = 0; i < bruker.getDims()[1]; i++) {
////                    for (int j = 0; j < 2*bruker.getDims()[2]; j=j+2) {
////                        INDArrayIndex[] indx = {NDArrayIndex.all(), NDArrayIndex.point(i), NDArrayIndex.point(j/2)};
////                        tempR = dataBruker.real.get(indx).toDoubleVector();
////                        tempI = dataBruker.imag.get(indx).toDoubleVector();
////                        for (int t = 0; t < bruker.getDims()[0]; t++) {
//////                            volumeW.data.set(j, i, 0, t, tempR[t]);
//////                            volumeW.data.set(j+1, i, 0, t, tempI[t]);
////                        }
////
////                    }
////                }}
//            }
//        }


        }

    private boolean isSVS() {
        try {
            if (bruker.getJcampdx().getVisu_pars().getArrayList("VisuCoreDimDesc")
                    .equals(new ArrayList() {{add("spectroscopic");}})) {
                return true;
            }
        } catch (Exception e) {
        }
        return false;
    }

    private INDArray getAffineMat() {
        // if it is CSI

        INDArray visuCoreSize = null;
        INDArray visuCoreExtent = null;
        if (!bruker.isImage()) {
            visuCoreSize = bruker.getJcampdx().getVisu_pars().getINDArray("VisuCoreSize").get(NDArrayIndex.interval(1, 3));
            visuCoreExtent = bruker.getJcampdx().getVisu_pars().getINDArray("VisuCoreExtent").get(NDArrayIndex.interval(1, 3));
        } else if (bruker.isImage()) {// if it is Image
        visuCoreSize = bruker.getJcampdx().getVisu_pars().getINDArray("VisuCoreSize");
        visuCoreExtent = bruker.getJcampdx().getVisu_pars().getINDArray("VisuCoreExtent");
        }
        // null checking
        INDArray resolution = visuCoreExtent.div(visuCoreSize);
        double VisuCoreFrameThickness = bruker.getJcampdx().getVisu_pars().getFloat("VisuCoreFrameThickness");
//        INDArray spatialResolution = Nd4j.append(resolution.get(NDArrayIndex.interval(1, 3)), 1, VisuCoreFrameThickness, -1);
        INDArray spatialResolution = Nd4j.append(resolution, 1, VisuCoreFrameThickness, -1);
        INDArray visuCoreOrientation = bruker.getJcampdx().getVisu_pars().getINDArray("VisuCoreOrientation");
        INDArray visuCorePosition = bruker.getJcampdx().getVisu_pars().getINDArray("VisuCorePosition");
        String visuSubjectPosition = bruker.getJcampdx().getVisu_pars().getString("VisuSubjectPosition");
        INDArray pvm_SPackArrGradOrient = bruker.getJcampdx().getMethod().getINDArray("PVM_SPackArrGradOrient");
        INDArray visuCoreUnits = bruker.getJcampdx().getVisu_pars().getINDArray("VisuCoreUnits");
        // to do check mm
        // add frame count to visuCoreSize (just for image)
        Float VisuCoreDataSlope = bruker.getJcampdx().getVisu_pars().getFloat("VisuCoreDataSlope");
        Float VisuCoreDataOffs = bruker.getJcampdx().getVisu_pars().getFloat("VisuCoreDataOffs");
        // slope Correction
        // offset Correction
        String[] subPos = new String[]{"Head_Prone", "Head_Supine"};
        if (Arrays.stream(subPos).anyMatch(e -> e.contains(visuSubjectPosition))) {
            System.out.println("Known case ('Head_Prone' or  'Head_Supine' for the parameter 'visu_pars.VisuSubjectPosition.)");
        }
        INDArray visuCoreOrientationReShaped = visuCoreOrientation.get(NDArrayIndex.point(0), NDArrayIndex.all()).reshape('f', new int[]{3, 3});
        INDArray result = Nd4j.eye(4);
        result.put(new INDArrayIndex[]{NDArrayIndex.interval(0, 3), NDArrayIndex.interval(0, 3)}, visuCoreOrientationReShaped);
        result.put(new INDArrayIndex[]{NDArrayIndex.interval(0, 3), NDArrayIndex.point(3)}, visuCorePosition.getRows(0));
//                  Invert the orientation matrix, according to nifti convention and Bruker manual.
        InvertMatrix.invert(result, true);
        // calculate determinant
        float[][] defMat = new float[][]{{1, 0, 0}, {0, 0, 1}, {0, 1, 0}};
        INDArray result_orientation =
                result.get(new INDArrayIndex[]{NDArrayIndex.interval(0, 3), NDArrayIndex.interval(0, 3)}).mmul(Nd4j.createFromArray(defMat));
//        from SAR to ASL
        if (pivot(result_orientation.getColumn(0).toDoubleVector()) > 0) {
            result_orientation.getColumn(0).muli(-1);
        }
        if (pivot(result_orientation.getColumn(1).toDoubleVector()) > 0) {
            result_orientation.getColumn(1).muli(-1);
        }
        if (pivot(result_orientation.getColumn(2).toDoubleVector()) < 0) {
            result_orientation.getColumn(2).muli(-1);
        }
        result_orientation.mmuli(Nd4j.diag(spatialResolution));
        result.put(new INDArrayIndex[]{NDArrayIndex.interval(0, 3), NDArrayIndex.interval(0, 3)}, result_orientation);
        // consider_subject_position:
        // keep same det
        return result;
    }

    public static Double pivot(double[] arr) {
        Double[] doubleArray = ArrayUtils.toObject(arr);
        List<Double> arrList = Arrays.asList(doubleArray);
        int indx = IntStream.range(0, arr.length).reduce((a, b) -> ((double) Math.abs(arrList.get(a))) < ((double) Math.abs(arrList.get(b))) ? b : a).getAsInt();
        return arr[indx];
    }

    public static INDArray calc_eulerangle(INDArray mati) {
        INDArray mat = mati.get(NDArrayIndex.point(0), NDArrayIndex.all(), NDArrayIndex.all());
        double sy = Math.sqrt(mat.getDouble(0, 0) * mat.getDouble(0, 0) + mat.getDouble(1, 0) * mat.getDouble(1, 0));
        double x, y, z;
        if (sy < 1e-6) {
            x = Math.atan2(-mat.getDouble(1, 2), mat.getDouble(1, 1));
            y = Math.atan2(-mat.getDouble(2, 0), sy);
            z = 0;
        } else {
            x = Math.atan2(mat.getDouble(2, 1), mat.getDouble(2, 2));
            y = Math.atan2(-mat.getDouble(2, 0), sy);
            z = Math.atan2(mat.getDouble(1, 0), mat.getDouble(0, 0));
        }
        return Nd4j.create(new double[]{x, y, z});
    }

}
